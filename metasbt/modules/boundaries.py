#!/usr/bin/env python3
"""
Define taxonomy-specific boundaries as the minimum and maximum number of kmers in common 
between all the genomes under a specific cluster
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 26, 2023"

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
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "Define taxonomy-specific boundaries as the minimum and maximum number of kmers in common "
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
        "--consider-mags",
        action="store_true",
        default=False,
        dest="consider_mags",
        help="Also consider MAGs while counting genomes",
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
        help="Consider genomes whose lineage belongs to a specific superkingdom",
    )
    p.add_argument(
        "--log",
        type=os.path.abspath,
        help="Path to the log file. Used to keep track of messages and errors printed on the stdout and stderr"
    )
    p.add_argument(
        "--max-genomes",
        type=number(int, minv=3),
        dest="max_genomes",
        help=(
            "Maximum number of genomes per cluster to be considered for computing boundaries. "
            "Genomes are selected randomly in case the size of clusters is greater than this number. "
            "This must always be greater than or equals to --min-genomes"
        ),
    )
    p.add_argument(
        "--min-genomes",
        type=number(int, minv=3),
        default=3,
        dest="min_genomes",
        help="Consider clusters with a minimum number of genomes only",
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
    level_dir: str,
    level_id: str,
    tmp_dir: str,
    output: str,
    kmer_len: int,
    filter_size: int,
    consider_mags: bool = False,
    max_genomes: Optional[int] = None,
    min_genomes: int = 3,
    nproc: int = 1,
) -> None:
    """
    Compute boundaries for the specified taxonomic level

    :param level_dir:       Path to the taxonomic level folder
    :param level_id:        ID of the taxonomic level
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param kmer_len:        Length of the kmers
    :param filter_size:     Size of the bloom filters
    :param consider_mags:   Consider MAGs while counting genomes
    :param max_genomes:     Consider this number of genomes at most for computing boundaries
    :param min_genomes:     Consider clusters with at least this number of genomes
    :param nproc:           Make the process parallel when possible
    """

    # Search and merge all the reference genomes paths under all references.txt files in the current taxonomic level
    samples: Dict[str, List[str]] = dict()

    # Genomes are usually listed in references.txt files
    references_paths = list(Path(level_dir).glob("**/references.txt"))

    # Databases without a taxonomic structure use genomes.txt
    references_paths.extend(list(Path(level_dir).glob("**/genomes.txt")))

    if consider_mags:
        # Extend the set of genomes to the MAGs
        references_paths.extend(list(Path(level_dir).glob("**/mags.txt")))

    # Take track of the overall number of reference genomes
    # under a specific taxonomic level
    references_count = 0

    for references_path in references_paths:
        path_split = str(references_path).split(os.sep)
        # In case the current level_id is species
        next_level = "NA"
        for path_pos, path_level in enumerate(path_split):
            if path_level.strip():
                if ("{}__".format(path_level[0]) == "{}__".format(level_id[0])) and level_id != "species":
                    next_level = path_split[path_pos + 1]
                    break

        if level_id != "species" and next_level not in samples:
            samples[next_level] = list()

        with open(str(references_path)) as references:
            for line in references:
                line = line.strip()
                if line:
                    if level_id == "species":
                        genome_path = os.path.join(
                            os.path.dirname(str(references_path)),
                            "strains",
                            "filters",
                            "{}.bf".format(line),
                        )

                        samples[line] = [genome_path]

                        references_count += 1

                    else:
                        genome_path = os.path.join(
                            os.path.dirname(str(references_path)),
                            "filters",
                            "{}.bf".format(line),
                        )

                        if os.path.isfile(genome_path):
                            samples[next_level].append(genome_path)

                            references_count += 1

    # In case the current taxonomic level is not the species level
    if level_id != "species":
        # Get rid of clusters with not enough genomes according to min_genomes
        for sample_id in list(samples.keys()):
            if len(samples[sample_id]) < min_genomes:
                # Redefine references count
                references_count -= len(samples[sample_id])

                # Remove cluster
                del samples[sample_id]

    # In case the number of genomes in the current taxonomic level
    # is greater than or equals to the minimum number of genomes specified in input
    if len(samples) >= min_genomes:
        # Create a temporary folder for the specific taxonomic level
        tmp_level_dir = os.path.join(tmp_dir, "boundaries", level_id, os.path.basename(level_dir))
        os.makedirs(tmp_level_dir, exist_ok=True)

        # Consider all the genomes under a particular level
        genome_paths = list()
        for tax_level_id in samples:
            genome_paths.extend(samples[tax_level_id])

        if max_genomes:
            # Always use the same seed for reproducibility
            rng = numpy.random.default_rng(0)

            # Shuffle the list of genome paths and get the first "max_genomes"
            rng.shuffle(genome_paths)
            genome_paths = genome_paths[:max_genomes]

        # Extract boundaries
        all_kmers, min_kmers, max_kmers = get_boundaries(
            genome_paths, tmp_level_dir, kmer_len, filter_size=filter_size, nproc=nproc
        )

        # Get the full lineage from the level folder path
        lineage = get_lineage_from_path(level_dir)

        # In case of --flat-structure
        if len(lineage.strip()) == 0:
            lineage = level_dir

        # Dump results to the boundaries table
        with open(output, "a+") as table:
            table.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    lineage,
                    len(samples) if level_id != "species" else 0,
                    references_count,
                    all_kmers,
                    min_kmers,
                    max_kmers,
                    round(min_kmers / all_kmers, 3),
                    round(max_kmers / all_kmers, 3),
                )
            )


def get_lineage_from_path(folder_path: str) -> str:
    """
    Rebuild a lineage from a taxonomic folder path

    :param folder_path: Path to the taxonomic folder
    :return:            Lineage
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
    db_dir: str,
    tmp_dir: str,
    output: str,
    flat_structure: bool = False,
    max_genomes: Optional[int] = None,
    min_genomes: int = 3,
    consider_mags: bool = False,
    superkingdom: Optional[str] = None,
    logger: Optional[Logger] = None,
    verbose: bool = False,
    nproc: int = 1,
) -> None:
    """
    Define boundaries for each of the taxonomic levels in the database
    Boundaries are defined as the minimum and maximum number of common kmers among all the reference genomes under a specific taxonomic level

    :param db_dir:          Path to the database root folder
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param flat_structure:  Genomes in the database have been organized without a taxonomic structure
    :param min_genomes:     Consider clusters with at least this number of genomes
    :param consider_mags:   Consider MAGs while counting genomes
    :param superkingdom:    Retrieve genomes that belong to a specific superkingdom
    :param logger:          Logger object
    :param verbose:         Print messages on screen
    :param nproc:           Make the process parallel when possible
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
        file.write("# --min-genomes {}\n".format(min_genomes))
        if max_genomes:
            file.write("# --max-genomes {}\n".format(max_genomes))
        file.write(
            "# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                "Lineage",  # Taxonomic label
                "Clusters",  # Number of clusters, valid for every taxonomic level except the species one
                "References",  # Number of reference genomes or clustrs under a specific taxonomic level
                "Kmers",  # Total number of kmers
                "Min kmers",  # Minimum number of common kmers among the reference genomes/clusters
                "Max kmers",  # Maximum number of common kmers among the reference genomes/clusters
                "Min score",  # Percentage of min kmers on the total number of genomes/clusters
                "Max score",  # Percentage of max kmers on the total number of genomes/clusters
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
            consider_mags=consider_mags,
            max_genomes=max_genomes,
            min_genomes=min_genomes,
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
                        species_manifest = load_manifest(species_manifest_filepath)

                    # Define boundaries for the current taxonomic level
                    define_boundaries(
                        str(level_dir),
                        level,
                        tmp_dir,
                        output,
                        manifest["kmer_len"],
                        species_manifest["filter_size"] if species_manifest else manifest["filter_size"],
                        consider_mags=consider_mags,
                        max_genomes=max_genomes,
                        min_genomes=min_genomes,
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
                consider_mags=consider_mags,
                max_genomes=max_genomes,
                min_genomes=min_genomes,
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

    # Check whether --max-genomes >= --min-genomes
    if args.max_genomes and args.max_genomes < args.min_genomes:
        raise ValueError("--max-genomes must always be greater than or equals to --min-genomes")

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
        max_genomes=args.max_genomes,
        min_genomes=args.min_genomes,
        consider_mags=args.consider_mags,
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
