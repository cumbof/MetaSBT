#!/usr/bin/env python3
"""Define taxonomy-specific boundaries as the minimum and maximum ANI similarities
between all the genomes under a specific cluster.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.4"
__date__ = "Oct 10, 2024"

import argparse as ap
import errno
import os
import shutil
import sys
import time
from datetime import datetime
from functools import partial
from logging import Logger
from pathlib import Path
from typing import Dict, List, Optional

import numpy  # type: ignore

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import build_sh, get_boundaries, init_logger, load_manifest, number, println  # type: ignore
except Exception:
    pass

# Define the module name
TOOL_ID = "boundaries"

# Define the list of dependencies
DEPENDENCIES: List[str] = list()

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Database folder path
    "--log",  # Path to the log file
    "--output",  # Output table path
    "--tmp-dir",  # Temporary folder path
]


def read_params():
    """Read and test the input arguments.

    Returns
    -------
    argparse.ArgumentParser
        The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "Define taxonomy-specific boundaries as the minimum and maximum ANI similarities "
            "between all the genomes under a specific cluster"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--cleanup",
        action="store_true",
        default=False,
        help="Remove temporary data at the end of the pipeline",
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required=True,
        dest="db_dir",
        help="This is the database directory with the taxonomically organised sequence bloom trees",
    )
    p.add_argument(
        "--flat-structure",
        action="store_true",
        default=False,
        dest="flat_structure",
        help="Genomes in the database have been organized without a taxonomic structure",
    )
    p.add_argument(
        "--superkingdom",
        type=str,
        help=(
            "Consider genomes whose lineage belongs to a specific superkingdom. "
            "Warning: it does not compute the boundaries of the specified kingdom"
        ),
    )
    p.add_argument(
        "--log",
        type=os.path.abspath,
        help="Path to the log file. Used to keep track of messages and errors printed on the stdout and stderr"
    )
    p.add_argument(
        "--max-nodes",
        type=number(int, minv=3),
        dest="max_nodes",
        help=(
            "Maximum number of nodes per cluster to be considered for computing boundaries. "
            "Nodes are selected randomly in case the size of clusters is greater than this number. "
            "This must always be greater than or equals to --min-nodes. "
            "WARNING: This works at the species level only!"
        ),
    )
    p.add_argument(
        "--min-nodes",
        type=number(int, minv=3),
        default=3,
        dest="min_nodes",
        help="Consider clusters with a minimum number of nodes only",
    )
    p.add_argument(
        "--nproc",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        help="This argument refers to the number of processors used for parallelizing the pipeline when possible",
    )
    p.add_argument(
        "--output",
        type=os.path.abspath,
        required=True,
        help="Output file with kmer boundaries for each of the taxonomic labels in the database",
    )
    p.add_argument(
        "--tmp-dir",
        type=os.path.abspath,
        required=True,
        dest="tmp_dir",
        help="Path to the folder for storing temporary data",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print messages and errors on the stdout"
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def define_boundaries(
    level_dir: os.path.abspath,
    level_id: str,
    tmp_dir: os.path.abspath,
    output: os.path.abspath,
    kmer_len: int,
    filter_size: int,
    max_nodes: Optional[int]=None,
    min_nodes: int=3,
    nproc: int=1,
) -> None:
    """Compute boundaries for the specified taxonomic level.

    Parameters
    ----------
    level_dir : os.path.abspath
        Path to the taxonomic level folder.
    level_id : str
        The ID of the taxonomic folder.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    output : os.path.abspath
        Path to the output table file with boundaries.
    kmer_len : int
        The length of the kmers.
    filter_size : int
        The size of the bloom filters.
    max_nodes : int, optional
        Consider this number of nodes at most for computing boundaries. This works at the species level only!
    min_nodes : int, default 3
        Consider clusters with at least this number of nodes.
    nproc : int, default 1
        Make the process parallel when possible.
    """

    bf_filepaths = list()

    if level_id == "species":
        genome_ids = list()

        if os.path.isfile(os.path.join(level_dir, "references.txt")):
            # Load the list of reference genomes
            genome_ids.extend([line.strip() for line in open(os.path.join(level_dir, "references.txt")).readlines() if line.strip()])

        if os.path.isfile(os.path.join(level_dir, "mags.txt")):
            # Load the list of MAGs
            genome_ids.extend([line.strip() for line in open(os.path.join(level_dir, "mags.txt")).readlines() if line.strip()])

        for genome_id in genome_ids:
            genome_path = os.path.join(level_dir, "strains", "filters", "{}.bf".format(genome_id))

            if not os.path.isfile(genome_id):
                # In this case, the database has been built by using only 3 representatives per species
                genome_path = os.path.join(level_dir, "filters", "{}.bf".format(genome_id))

            bf_filepaths.append(genome_path)

    else:
        levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

        # Define the next taxonomic level
        next_level_id = levels[levels.index(level_id.lower())+1]

        # Get the root nodes of all the SBTs under a specific taxonomic level
        for filepath in Path(os.path.join(level_dir, "index")).glob("{}__*.bf".format(next_level_id[0])):
            bf_filepaths.append(
                os.path.join(
                    level_dir,
                    os.path.basename(filepath).split(".")[0],  # e.g., p__Nucleocytoviricota.detbrief.rrr.bf
                    "{}.bf".format(os.path.basename(filepath).split(".")[0])
                )
            )

    # In case the number of BFs in the current taxonomic level
    # is greater than or equals to the minimum number of nodes specified in input
    if len(bf_filepaths) >= min_nodes:
        nodes_count = len(bf_filepaths)

        if max_nodes and level_id == "species" and len(bf_filepaths) > max_nodes:
            # Pick `max_nodes` random nodes (only at the species level)
            # Always use the same seed for reproducibility
            rng = numpy.random.default_rng(0)

            # Shuffle the list of genome paths and get the first `max_nodes`
            rng.shuffle(bf_filepaths)

            bf_filepaths = bf_filepaths[:max_nodes]

        # Create a temporary folder for the specific taxonomic level
        tmp_level_dir = os.path.join(tmp_dir, "boundaries", level_id, os.path.basename(level_dir))

        os.makedirs(tmp_level_dir, exist_ok=True)

        # Extract boundaries
        min_score, max_score, centroid = get_boundaries(
            bf_filepaths,
            tmp_level_dir,
            kmer_len,
            filter_size=filter_size,
            nproc=nproc,
        )

        # Get the full lineage from the level folder path
        lineage = get_lineage_from_path(level_dir)

        # In case of --flat-structure
        if len(lineage.strip()) == 0:
            lineage = level_dir

        # Dump results to the boundaries table
        with open(output, "a+") as table:
            table.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    lineage,
                    nodes_count,
                    min_score,
                    max_score,
                    centroid,
                )
            )


def get_lineage_from_path(folder_path: os.path.abspath) -> str:
    """Rebuild a lineage from a taxonomic folder path.

    Parameters
    ----------
    folder_path : os.path.abspath
        Path to the taxonomic folder.

    Returns
    -------
    str
        The taxonomic label.

    Examples
    --------
    >>> from boundaries import get_lineage_from_path
    >>> folder_path = "~/MetaSBT-DBs/k__Viruses/p__Nucleocytoviricota/c__Pokkesviricetes/o__Chitovirales/f__Poxviridae/g__Orthopoxvirus/s__Monkeypox_virus"
    >>> lineage = get_lineage_from_path(folder_path)
    >>> print(lineage)
    k__Viruses|p__Nucleocytoviricota|c__Pokkesviricetes|o__Chitovirales|f__Poxviridae|g__Orthopoxvirus|s__Monkeypox_virus

    Define a taxonomic folder as it appears when indexing Viruses with MetaSBT.
    In this example, the folder points to the Monkeypox virus, and its full lineage is rebuilt starting from the folder path.
    """

    lineage_list = list()

    superkingdom_found = False

    for level in folder_path.split(os.sep):
        # Search for the superkingdom level
        if level.startswith("k__"):
            superkingdom_found = True

        # Start rebuilding the lineage form the superkingdom
        if superkingdom_found:
            lineage_list.append(level)

    return "|".join(lineage_list)


def boundaries(
    db_dir: os.path.abspath,
    tmp_dir: os.path.abspath,
    output: os.path.abspath,
    flat_structure: bool = False,
    max_nodes: Optional[int] = None,
    min_nodes: int = 3,
    superkingdom: Optional[str] = None,
    logger: Optional[Logger] = None,
    verbose: bool = False,
    nproc: int = 1,
) -> None:
    """Define boundaries for each of the taxonomic levels in the database. 
    Boundaries are defined as the minimum and maximum ANI similarities among all the genomes under a specific taxonomic level.

    Parameters
    ----------
    db_dir : os.path.abspath
        Path to the database root folder.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    output : os.path.abspath
        Path to the output table file with boundaries.
    flat_structure : bool, default False
        If True, genomes in the database have been organized without a taxonomic structure.
    max_nodes : int, optional
        Consider this number of nodes at most for computing boundaries. This works at the species level only!
    min_nodes : int, default 3
        Consider clusters with at least this number of nodes.
    superkingdom : str, optional
        Retrieve genomes that belong to a specific superkingdom.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.
    nproc : int, default 1
        Make the process parallel when possible.

    Raises
    ------
    FileNotFoundError
        If the "manifest.txt" file with the database info is not available.
    Exception
        If the "manifest.txt" file does not contain information about the kmer length and filter size.
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Start defining the output table file
    with open(output, "w+") as file:
        # Write header lines
        file.write("# {} version {} ({})\n".format(TOOL_ID, __version__, __date__))

        file.write("# timestamp: {}\n".format(datetime.today().strftime("%Y%m%d")))

        file.write("# --db-dir {}\n".format(db_dir))

        if superkingdom:
            file.write("# --superkingdom {}\n".format(superkingdom))

        file.write("# --min-nodes {}\n".format(min_nodes))

        if max_nodes:
            file.write("# --max-nodes {}\n".format(max_nodes))

        file.write(
            "# {}\t{}\t{}\t{}\t{}\n".format(
                "Lineage",  # Taxonomic label
                "Nodes",  # Number of genomes/clusters under a specific taxonomic level
                "Min ANI",  # The left boundary of a specific lineage (see definition of boundaries)
                "Max ANI",  # The right boundary of a specific lineage (see definition of boundaries)
                "Centroid",  # The cluster centroid defined as the genome that maximize the similarity with all the other genomes
            )
        )

    # Check whether the manifest file exists
    manifest_filepath = os.path.join(db_dir, "manifest.txt")

    if not os.path.isfile(manifest_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    # Load the manifest file
    manifest = load_manifest(manifest_filepath)

    if "kmer_len" not in manifest or "filter_size" not in manifest:
        raise Exception(
            "Manifest file does not contain --kmer-len and --filter-size information: {}".format(manifest_filepath)
        )

    # Check whether the genomes folder exists under the database root directory
    if flat_structure:
        # This means that the database has been build with the --flat-structure option
        printline("Defining boundaries")

        # Treat the database as the species level
        define_boundaries(
            db_dir,
            "species",
            tmp_dir,
            output,
            manifest["kmer_len"],
            manifest["filter_size"],
            max_nodes=max_nodes,
            min_nodes=min_nodes,
            nproc=nproc,
        )

    else:
        # Genomes have been taxonomically organized
        target_dir = db_dir if not superkingdom else os.path.join(db_dir, "k__{}".format(superkingdom))

        levels = ["species", "genus", "family", "order", "class", "phylum"]

        if not superkingdom:
            levels.append("kingdom")

        # Iterate over the taxonomic levels
        for level in levels:
            printline("Defining {} boundaries".format(level))

            for level_dir in Path(target_dir).glob("**/{}__*".format(level[0])):
                if os.path.isdir(str(level_dir)):
                    printline("Processing {}".format(get_lineage_from_path(str(level_dir))))

                    # Load the species manifest in case of the species level
                    species_manifest = dict()

                    if level == "species":
                        # Load the manifest file under the strains folder
                        # Use the species-specific bloom filter size
                        species_manifest_filepath = os.path.join(str(level_dir), "strains", "manifest.txt")

                        if os.path.isfile(species_manifest_filepath):
                            species_manifest = load_manifest(species_manifest_filepath)

                    # Define boundaries for the current taxonomic level
                    define_boundaries(
                        str(level_dir),
                        level,
                        tmp_dir,
                        output,
                        manifest["kmer_len"],
                        species_manifest["filter_size"] if species_manifest else manifest["filter_size"],
                        max_nodes=max_nodes,
                        min_nodes=min_nodes,
                        nproc=nproc,
                    )

        if superkingdom:
            # Also define boundaries for the specified superkingdom
            define_boundaries(
                os.path.join(db_dir, superkingdom),
                "kingdom",
                tmp_dir,
                output,
                manifest["kmer_len"],
                manifest["filter_size"],
                max_nodes=max_nodes,
                min_nodes=min_nodes,
                nproc=nproc,
            )

    # Report the path to the output boundaries table file
    printline("Output table: {}".format(output))


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether the database folder exists
    target_dir = args.db_dir if not args.superkingdom else os.path.join(args.db_dir, "k__{}".format(args.superkingdom))

    if not os.path.isdir(target_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), target_dir)

    # Check whether the output boundaries table alrady exists
    if os.path.isfile(args.output):
        raise Exception("The output boundaries table already exists")

    if args.max_nodes and args.max_nodes < args.min_nodes:
        raise ValueError("--max-nodes must always be greater than or equals to --min-nodes")

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Build a sh script with the command line used to launch the boundaries module
    build_sh(sys.argv, TOOL_ID, args.db_dir)

    t0 = time.time()

    boundaries(
        args.db_dir,
        args.tmp_dir,
        args.output,
        flat_structure=args.flat_structure,
        max_nodes=args.max_nodes,
        min_nodes=args.min_nodes,
        superkingdom=args.superkingdom,
        logger=logger,
        verbose=args.verbose,
        nproc=args.nproc,
    )

    if args.cleanup:
        # Remove the temporary folder
        println(
            "Cleaning up temporary space",
            logger=logger,
            verbose=args.verbose,
        )

        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()

    println(
        "Total elapsed time {}s".format(int(t1 - t0)),
        logger=logger,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
