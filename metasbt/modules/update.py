#!/usr/bin/env python3
"""Update a specific database with a new set of reference genomes or metagenome-assembled genomes.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.2"
__date__ = "Mar 27, 2024"

import argparse as ap
import errno
import math
import multiprocessing as mp
import os
import shutil
import sys
import time
from collections import Counter
from functools import partial
from itertools import takewhile
from logging import Logger
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import tqdm  # type: ignore

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import (  # type: ignore  # isort: skip
        bfaction,
        build_sh,
        cluster,
        dereplicate_genomes,
        estimate_bf_size,
        filter_quality,
        get_file_info,
        get_level_boundaries,
        howdesbt,
        init_logger,
        load_boundaries,
        load_input_table,
        load_manifest,
        number,
        println,
        quality,
        run,
    )
except Exception:
    pass

# Define the module name
TOOL_ID = "update"

# Define the list of dependencies
DEPENDENCIES = [
    "checkm2",
    "checkv",
    "eukcc",
    "howdesbt",
    "ntcard"
]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Database folder path
    "--input-list",  # Input file path
    "--log",  # Path to the log file
    "--tmp-dir",  # Temporary folder path
]

# Define the software root directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def read_params(argv: List[str]):
    """Read and test the input arguments.

    Parameters
    ----------
    argv : list
        List with argument names and values.

    Returns
    -------
    argparse.ArgumentParser
        The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Update a specific database with a new set of reference genomes or metagenome-assembled genomes",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )

    # General arguments
    general_group = p.add_argument_group("General arguments")

    general_group.add_argument(
        "--boundaries",
        type=os.path.abspath,
        required = "--resume" not in argv,
        help="Path to the output table produced by the boundaries module",
    )
    general_group.add_argument(
        "--boundary-uncertainty",
        type=number(float, minv=0.0, maxv=100.0),
        default=0.0,
        dest="boundary_uncertainty",
        help="Define the percentage of kmers to enlarge and reduce boundaries",
    )
    general_group.add_argument(
        "--cluster-prefix",
        type=str,
        default="MSBT",
        dest="cluster_prefix",
        help="Prefix of clusters numerical identifiers",
    )
    general_group.add_argument(
        "--cleanup",
        action="store_true",
        default=False,
        help="Remove temporary data at the end of the pipeline",
    )
    general_group.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required = "--resume" not in argv,
        dest="db_dir",
        help="This is the database directory with the taxonomically organised sequence bloom trees",
    )
    general_group.add_argument(
        "--extension",
        type=str,
        required = "--resume" not in argv,
        choices=["fa", "fa.gz", "fasta", "fasta.gz", "fna", "fna.gz"],
        help=(
            "Specify the input genome files extension. "
            "All the input genomes must have the same file extension before running this module"
        ),
    )
    general_group.add_argument(
        "--input-list",
        type=os.path.abspath,
        required = "--resume" not in argv,
        dest="input_list",
        help=(
            "Path to the input table with a list of genome file paths and an optional column with their taxonomic labels. "
            "Please note that the input genome files must all have the same extension and can be Gzip compressed (e.g.: *.fna.gz)"
        ),
    )
    general_group.add_argument(
        "--log",
        type=os.path.abspath,
        help="Path to the log file. Used to keep track of messages and errors printed on the stdout and stderr"
    )
    general_group.add_argument(
        "--nproc",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        help="This argument refers to the number of processors used to run the pipeline in parallel",
    )
    general_group.add_argument(
        "--parallel",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        help=(
            "Used to rebuild multiple clusters in parallel. "
            "Warning: The total number of spawned processes is equals to --parallel * --nproc"
        ),
    )
    general_group.add_argument(
        "--parent-version",
        type=str,
        dest="parent_version",
        help="Version of MetaSBT that called this process"
    )
    general_group.add_argument(
        "--quality-control",
        type=str.lower,
        choices=["checkm2", "checkv", "eukcc"],
        dest="quality_control",
        help="Quality control method to be used in conjunction with --completeness and --contamination",
    )
    general_group.add_argument(
        "--resume",
        type=os.path.abspath,
        help=(
            "Path to the \"update.sh\" file with the configuration of a previous run. "
            "Used to resume the update process in case of an unexpected error"
        ),
    )
    general_group.add_argument(
        "--tmp-dir",
        type=os.path.abspath,
        required = "--resume" not in argv,
        dest="tmp_dir",
        help="Path to the folder for storing temporary data",
    )
    general_group.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print messages and errors on the stdout"
    )
    general_group.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )

    # Group of arguments for quality control
    qc_group = p.add_argument_group("Quality control")

    qc_group.add_argument(
        "--completeness",
        type=number(float, minv=0.0, maxv=100.0),
        default=0.0,
        help="Input genomes must have a minimum completeness percentage before being processed and added to the database",
    )
    qc_group.add_argument(
        "--contamination",
        type=number(float, minv=0.0, maxv=100.0),
        default=100.0,
        help="Input genomes must have a maximum contamination percentage before being processed and added to the database",
    )

    # Group of arguments for the dereplication
    dereplication_group = p.add_argument_group("Dereplication of genomes")

    dereplication_group.add_argument(
        "--dereplicate",
        action="store_true",
        default=False,
        help="Enable the dereplication of genomes",
    )
    dereplication_group.add_argument(
        "--similarity",
        type=number(float, minv=0.0, maxv=100.0),
        default=100.0,
        help=(
            "Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. "
            "This is used exclusively in conjunction with the --dereplicate argument"
        ),
    )

    return p.parse_args(argv)


def profile_and_assign(
    genome_path: os.path.abspath,
    input_type: str,
    tmp_dir: os.path.abspath,
    db_dir: os.path.abspath,
    tmp_genomes_dir: os.path.abspath,
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    boundary_uncertainty: float=0.0,
    taxonomies: Optional[dict]=None,
    dereplicate: bool=False,
    similarity: float=100.0,
    quality_dict: Optional[Dict[str, Dict[str, str]]]=None,
    use_representatives: bool=False,
    logger: Optional[Logger]=None,
    verbose: bool=False,
) -> Tuple[dict, list, list]:
    """Profile and assign one input genome.

    Parameters
    ----------
    genome_path : os.path.abspath
        Path to the input genome file.
    input_type : {"MAGs", "references"}
        The nature of the input genomes.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    db_dir : os.path.abspath
        Path to the database root folder.
    tmp_genomes_dir : os.path.abspath
        Path to the temporary folder for uncompressing input genomes.
    boundaries : dict
        Content of the boundaries table as produced by the boundaries module.
    boundary_uncertainty : float, default 0.0
        Percentage of kmers to enlarge and reduce boundaries.
    taxonomies : dict, optional
        Dictionary with mapping between input genome names and their taxonomic labels (for reference genomes only).
    dereplicate : bool, default False
        Enable dereplication of input genome against genomes in the closest cluster.
    similarity : float, default 100.0
        Exclude input genome if its similarity with the closest genome in the database exceed this threshold.
    quality_dict : dict, optional
        Dictionary with quality stats.
    use_representatives : bool, default False
        The database has been built by using only 3 representatives per species.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.

    Returns
    -------
    tuple
        A tuple with the list of taxonomies that now contain at least one reference genome,
        the list of taxonomies that must be rebuilt, and the list of unassigned genomes.
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    if not taxonomies:
        taxonomies = dict()

    if not quality_dict:
        quality_dict = dict()

    in_strains = "strains" if use_representatives else ""

    # The input genome is added to the "to_known_taxa" in case it is a reference genome and the closest cluster is unknown
    to_known_taxa: Dict[str, List[str]] = dict()

    # Add its closest taxonomic label in case the input genome is assigned to that taxonomy
    rebuild: List[str] = list()

    # Otherwise, the input genome is unassigned
    unassigned: List[str] = list()

    # Get genome file info
    _, genome_name, genome_ext, compression = get_file_info(genome_path)

    filepath = genome_path

    if compression:
        filepath = os.path.join(tmp_genomes_dir, "{}{}".format(genome_name, genome_ext))

        if not os.path.isfile(filepath):
            # Unzip the genome file into the tmp folder
            with open(filepath, "w+") as genome_file:
                run(["gunzip", "-c", genome_path], stdout=genome_file)

    printline("Profiling {}".format(genome_name))

    # Define the path to the profiler output file
    profile_path = os.path.join(tmp_dir, "profiling", "{}__profiles.tsv".format(genome_name))

    if not os.path.isfile(profile_path):
        # Run the profiler to establish the closest genome and the closest group
        # for each taxonomic level in the tree
        # The profile module is in the same folder of the update module
        profiler = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "profile.py"),
            "--input-file",
            filepath,
            "--input-id",
            genome_name,
            "--input-type",
            "genome",
            "--tree",
            os.path.join(db_dir, "index", "index.detbrief.sbt"),
            "--expand",
            "--output-dir",
            os.path.join(tmp_dir, "profiling"),
            "--output-prefix",
            genome_name,
        ]

        run(profiler, silence=True)

    printline(profile_path)

    # Check whether the profile exists
    if os.path.isfile(profile_path):
        # Load the profile
        profile_data: Dict[str, Dict[str, Any]] = dict()

        with open(profile_path) as profile:
            for line in profile:
                line = line.strip()
                if line:
                    if not line.startswith("#"):
                        line_split = line.split("\t")
                        profile_data[line_split[1]] = {
                            "taxonomy": line_split[2],
                            "common_kmers": int(line_split[3].split("/")[0]),
                            "score": float(line_split[4]),
                        }

        closest_genome: Dict[str, Any] = profile_data["genome"]

        # Get the closest species
        closest_taxa = profile_data["species"]["taxonomy"]
        closest_common_kmers = profile_data["species"]["common_kmers"]
        closest_score = profile_data["species"]["score"]

        printline("Closest lineage: {} (score {})".format(closest_taxa, closest_score))
        printline("Closest genome: {} (score {})".format(closest_genome["taxonomy"], closest_genome["score"]))

        # Define the folder path of the closest taxonomy under the database
        closest_taxadir = os.path.join(db_dir, closest_taxa.replace("|", os.sep))

        # Load the set of reference genomes that belongs to the closest species
        references_filepath = os.path.join(closest_taxadir, "references.txt")
        references = list()

        if os.path.isfile(references_filepath):
            references = [ref.strip() for ref in open(references_filepath).readlines() if ref.strip()]

        # Also load the set of MAGs that belongs to the closest species
        mags_filepath = os.path.join(closest_taxadir, "mags.txt")
        mags = list()

        if os.path.isfile(mags_filepath):
            mags = [mag.strip() for mag in open(mags_filepath).readlines() if mag.strip()]

        # Check whether the input genome must be discarded
        skip_genome = False

        # Check whether the quality table must be updated
        add_quality_info = False

        if dereplicate and closest_genome["score"] * 100.0 >= similarity:
            # Discard the genome if it is too similar with the closest genome according to the similarity score
            # Also, do not discard the input genome if it is a reference genome and the closest genome in the database is
            # a MAG with a high number of overlapped kmers according to the similarity score
            printline("Dereplicating genome")

            if input_type == "MAGs":
                # In case the input genome is a MAG
                if closest_genome["taxonomy"] in references:
                    # In case the closest genome is a reference genome
                    skip_genome = True
                    # Print the reason why the current genome is excluded
                    printline("Discarding genome\nInput genome is a MAG and the closest genome is a reference")

                elif closest_genome["taxonomy"] in mags:
                    # In case the closest genome is a MAG
                    skip_genome = True
                    # Print the reason why the current genome is excluded
                    printline("Discarding genome\nInput genome and the closest genome are both MAGs")

            elif input_type == "references":
                # In case the input genome is a reference genome
                if closest_genome["taxonomy"] in references:
                    # In case the closest genome is a reference genome
                    skip_genome = True
                    # Print the reason why the current genome is excluded
                    printline("Discarding genome\nInput genome and the closest genome are both reference genomes")

        # In case the input genome survived the dereplication with the closest genome in the database
        if not skip_genome:
            # Retrieve the minimum number of common kmers for the closest taxa
            min_bound, _ = get_level_boundaries(boundaries, closest_taxa)
            # Add uncertainty to the boundaries
            min_bound -= int(min_bound * boundary_uncertainty / 100.0)

            if input_type == "MAGs":
                # In case the input genome is a MAG
                if closest_common_kmers >= min_bound:
                    if not os.path.isfile(os.path.join(closest_taxadir, in_strains, "genomes", os.path.basename(genome_path))):
                        # Assign the current genome to the closest lineage
                        shutil.copy(genome_path, os.path.join(closest_taxadir, in_strains, "genomes"))

                    genome_filepath = os.path.join(closest_taxadir, in_strains, "genomes", "{}{}".format(genome_name, genome_ext))

                    if not compression:
                        if not os.path.isfile("{}{}".format(genome_filepath, ".gz")):
                            # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                            run(["gzip", genome_filepath], silence=True)

                        else:
                            # There is a genome with the same name or the update module has been run with the --resume option
                            # We can get rid of the uncompressed version of the genome file
                            os.unlink(genome_filepath)

                    mags_entries = list()

                    if os.path.isfile(mags_filepath):
                        mags_entries = [line.strip() for line in open(mags_filepath).readlines() if line.strip()]

                    if not genome_name in mags_entries:
                        # Add the current genome to the list of MAGs
                        with open(mags_filepath, "a+") as mags_file:
                            mags_file.write("{}\n".format(genome_name))

                    # Also add its quality statistics if available
                    add_quality_info = True

                    # Do not remove the index here because there could be other input genome
                    # with the same closest taxonomic label
                    rebuild.append(closest_taxa)

                    printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                else:
                    # Mark the input genome as unassigned
                    # The quality stats for this genome will be reported after the assignment
                    unassigned.append(genome_path)

                    printline("{} is still unassigned".format(genome_name))

            elif input_type == "references":
                # In case the input genome is a reference genome
                # Retrieve the taxonomic label from the input mapping file
                taxalabel = taxonomies[genome_name]

                if closest_common_kmers >= min_bound:
                    # Check whether the closest taxonomy contains any reference genome
                    how_many_references = 0
                    if os.path.isfile(references_filepath):
                        how_many_references = sum([1 for ref in open(references_filepath).readlines() if ref.strip()])

                    if how_many_references == 0:
                        # If the closest genome belongs to a new cluster with no reference genomes
                        # Assign the current reference genome to the new cluster and rename its lineage with the taxonomic label of the reference genome
                        # Do not change lineage here because there could be more than one reference genome assigned to the same unknown cluster
                        # Assign a new taxonomy according to a majority rule applied on the reference genomes taxa
                        if closest_taxa not in to_known_taxa:
                            to_known_taxa[closest_taxa] = list()
                        to_known_taxa[closest_taxa].append(genome_path)

                        printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                    elif how_many_references >= 1:
                        # If the closest genome belongs to a cluster with at least one reference genome
                        if taxalabel != closest_taxa:
                            # If the taxonomic labels of the current reference genome and that one of the closest genome do not match
                            # Report the inconsistency and continue
                            printline(
                                "Inconsistency found for genome {}:\nInput lineage: {}\nClosest lineage: {}".format(
                                    genome_name, taxalabel, closest_taxa
                                )
                            )

                        # Assign the current genome to the closest lineage
                        if not os.path.isfile(os.path.join(closest_taxadir, in_strains, "genomes", os.path.basename(genome_path))):
                            shutil.copy(genome_path, os.path.join(closest_taxadir, in_strains, "genomes"))

                        genome_filepath = os.path.join(closest_taxadir, in_strains, "genomes", "{}{}".format(genome_name, genome_ext))

                        if not compression:
                            if not os.path.isfile("{}{}".format(genome_filepath, ".gz")):
                                # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                                run(["gzip", genome_filepath], silence=True)

                            else:
                                # There is a genome with the same name or the update module has been run with the --resume option
                                # We can get rid of the uncompressed version of the genome file
                                os.unlink(genome_filepath)

                        references_entries = [line.strip() for line in open(references_filepath).readlines() if line.strip()]

                        if not genome_name in references_entries:
                            # Add the current genome to the list of reference genomes
                            with open(references_filepath, "a+") as references_file:
                                references_file.write("{}\n".format(genome_name))

                        # Do not remove the index here because there could be other input genome
                        # with the same closest taxonomic label
                        rebuild.append(closest_taxa)

                        printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                    # Also add its quality statistics if available
                    add_quality_info = True

                else:
                    # If nothing is closed enough to the current genome and its taxonomic label does not exist in the database
                    # Create a new branch with the new taxonomy for the current genome
                    taxdir = os.path.join(db_dir, taxalabel.replace("|", os.sep))
                    os.makedirs(os.path.join(taxdir, "genomes"), exist_ok=True)

                    # Also create the strains folder
                    os.makedirs(os.path.join(taxdir, in_strains, "genomes"), exist_ok=True)

                    # Assign the current genome to the closest lineage
                    if not os.path.isfile(os.path.join(taxdir, in_strains, "genomes", os.path.basename(genome_path))):
                        shutil.copy(genome_path, os.path.join(taxdir, in_strains, "genomes"))

                    genome_filepath = os.path.join(taxdir, in_strains, "genomes", "{}{}".format(genome_name, genome_ext))
                    if not compression and not os.path.isfile("{}{}".format(genome_filepath, compression)):
                        # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                        run(["gzip", genome_filepath], silence=True)

                    references_entries = [line.strip() for line in open(os.path.join(taxdir, "references.txt")).readlines() if line.strip()]

                    if not genome_name in references_entries:
                        # Add the current genome to the list of reference genomes
                        with open(os.path.join(taxdir, "references.txt"), "a+") as references_file:
                            references_file.write("{}\n".format(genome_name))

                    # Do not remove the index here because there could be other input genome
                    # with the same closest taxonomic label
                    rebuild.append(closest_taxa)

                    printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                    # Also add its quality statistics if available
                    add_quality_info = True

        if add_quality_info:
            # Also add its quality statistics if available
            if genome_name in quality_dict:
                quality_filepath = os.path.join(closest_taxadir, "quality.tsv")
                add_to_quality_table = True

                if os.path.isfile(quality_filepath):
                    genomes_in_quality_file = list()

                    with open(quality_filepath) as table:
                        header = table.readline().strip()
                        for line in table:
                            line = line.strip()
                            if line:
                                line_split = line.split("\t")
                                genomes_in_quality_file.append(line_split[header.index("Name")])
                    
                    if genome_name in genomes_in_quality_file:
                        add_to_quality_table = False

                if add_to_quality_table:
                    not_exists = not os.path.isfile(quality_filepath)

                    with open(quality_filepath, "a+") as quality_file:
                        header = sorted(list(quality_dict[genome_name].keys()))

                        if not_exists:                                
                            quality_file.write("{}\n".format("\t".join(header)))

                        quality_file.write("{}\n".format("\t".join([quality_dict[genome_name][h] for h in header])))

    # Remove the uncompressed version of the input genome in the temporary folder
    if os.path.isfile(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0])):
        os.unlink(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0]))

    return to_known_taxa, rebuild, unassigned


def build_cluster(
    taxonomy: str,
    db_dir: Optional[os.path.abspath]=None,
    tmp_dir: Optional[os.path.abspath]=None,
    kmer_len: int=4,
    filter_size: int=2,
    min_kmer_occurrences: int=2,
    nproc: int=1,
    use_representatives: bool=False,
    logger: Optional[Logger]=None,
    verbose: bool=False,
) -> Tuple[str, bool]:
    """Build a cluster under a specific taxonomic level.

    Parameters
    ----------
    taxonomy : str
        Taxonomic label of the cluster up to a specific taxonomic level.
    db_dir : os.path.abspath, optional
        Path to the root folder of the database.
    tmp_dir : os.path.abspath, optional
        Path to the temporary folder.
    kmer_len : int, default 4
        Size of the kmers.
    filter_size : int, default 2
        Global bloom filter size.
    min_kmer_occurrences : int, default 2
        Minimum number of kmers occurrences.
    nproc : int, default 1
        Make it parallel.
    use_representatives : bool, default False
        The database has been built by using only 3 representatives per species.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.

    Returns
    -------
    tuple
        A tuple with the taxonomic label and a boolean flag that is True in case the clusters
        at the upper taxonomic levels should not be rebuilt.
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    build = True
    remove_from_rebuild = False

    if verbose:
        printline("\t{}".format(taxonomy))
    
    tax_dir = os.path.join(db_dir, taxonomy.replace("|", os.sep))

    if taxonomy.split("|")[-1].startswith("s__") and use_representatives:
        # In case the taxonomy is complete at the species level
        # Load manifest
        strains_manifest_filepath = os.path.join(tax_dir, "strains", "manifest.txt")

        if not os.path.isfile(strains_manifest_filepath):
            # Estimate a proper bloom filter size
            strains_filter_size = estimate_bf_size(
                genomes=list(Path(os.path.join(tax_dir, "strains", "genomes")).glob("*.gz")),
                kmer_len=kmer_len,
                min_occurrences=min_kmer_occurrences,
                prefix="genomes",
                tmp_dir=os.path.join(tax_dir, "strains", "tmp"),
                nproc=nproc,
            )

            with open(strains_manifest_filepath, "w+") as strains_manifest_file:
                strains_manifest_file.write("--min-kmer-occurrences {}\n".format(min_kmer_occurrences))
                strains_manifest_file.write("--kmer-len {}\n".format(kmer_len))
                strains_manifest_file.write("--filter-size {}\n".format(strains_filter_size))

        strains_manifest = load_manifest(strains_manifest_filepath)

        # Rebuild the strains SBT first
        # Use the species-specific bloom filter size in manifest
        howdesbt(
            os.path.join(tax_dir, "strains"),
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
            filter_size=strains_manifest["filter_size"],
            nproc=nproc,
            flat_structure=True,
        )

        # Select the representative genomes
        selected_genomes = list()

        # Get the bloom filters file paths
        bf_filepaths = [str(path) for path in Path(os.path.join(tax_dir, "strains", "filters")).glob("*.bf")]

        if len(bf_filepaths) <= 3:
            # 3 is the maximum number of selected species
            # as it is also the minimum number of genomes for computing boundaries
            selected_genomes = [
                get_file_info(bf_path, check_supported=False, check_exists=False)[1] for bf_path in bf_filepaths
            ]

        else:
            # Compute the theta distance between genomes
            bfdistance_theta = bfaction(
                bf_filepaths,
                os.path.join(tmp_dir, "howdesbt_strains"),
                kmer_len,
                min_occurrences=min_kmer_occurrences,
                filter_size=strains_manifest["filter_size"],
                nproc=nproc,
                action="bfdistance",
                mode="theta"
            )

            # Sum the distances to get the final score
            bfdistance_sums = {genome: sum(bfdistance_theta[genome].values()) for genome in bfdistance_theta if genome in bfdistance_theta}

            # Sort genomes according to the sum of their distances with all the other genomes
            sorted_genomes = sorted(bfdistance_sums, key=lambda genome: bfdistance_sums[genome])

            # First and last genomes are those that minimize and maximize the distance with all the other genomes
            selected_genomes.append(sorted_genomes[0])
            selected_genomes.append(sorted_genomes[-1])

            # Also select a genome in the middle of min and max distances
            selected_genomes.append(sorted_genomes[math.ceil(int(len(sorted_genomes) / 2))])

        # Get the list of previously selected representative genomes
        previously_selected_genomes = [
            get_file_info(str(path))[1] for path in Path(os.path.join(tax_dir, "filters")).glob("*.bf")
        ]

        add_selected_genomes = set(selected_genomes).difference(set(previously_selected_genomes))
        remove_selected_genomes = set(previously_selected_genomes).difference(set(selected_genomes))                        

        if len(remove_selected_genomes) > 0:
            # Remove old representatives
            for genome_path in Path(os.path.join(tax_dir, "genomes")).glob("*.gz"):
                _, genome_name, _, _ = get_file_info(str(genome_path))

                if genome_name in remove_selected_genomes:
                    os.unlink(str(genome_path))

            for genome_name in remove_selected_genomes:
                os.unlink(os.path.join(tax_dir, "filters", "{}.bf".format(genome_name)))
        
        if len(add_selected_genomes) > 0:
            # Create the genomes folder under the species level in case of a new species
            os.makedirs(os.path.join(tax_dir, "genomes"), exist_ok=True)

            # Add new representatives
            for genome_path in Path(os.path.join(tax_dir, "strains", "genomes")).glob("*.gz"):
                _, genome_name, _, _ = get_file_info(str(genome_path))

                if genome_name in add_selected_genomes:
                    os.symlink(str(genome_path), os.path.join(tax_dir, "genomes", os.path.basename(str(genome_path))))

        if not add_selected_genomes and not remove_selected_genomes:
            # In case the representative genomes of the species did not change
            # Remove the taxonomy from the rebuild list
            remove_from_rebuild = True
            build = False

    if build:
        # Remove the old index if it exists
        shutil.rmtree(os.path.join(tax_dir, "index"), ignore_errors=True)

        if not taxonomy.split("|")[-1].startswith("s__"):
            # Remove filters
            shutil.rmtree(os.path.join(tax_dir, "filters"), ignore_errors=True)

        # Also remove the copy of the root node if it exists
        if os.path.isfile(os.path.join(tax_dir, "{}.bf".format(os.path.basename(tax_dir)))):
            os.unlink(os.path.join(tax_dir, "{}.bf".format(os.path.basename(tax_dir))))

        # Rebuild the index with HowDeSBT
        howdesbt(
            tax_dir,
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
            filter_size=filter_size,
            nproc=nproc,
            flat_structure=False,
        )

    return taxonomy, remove_from_rebuild


def update(
    input_list: os.path.abspath,
    extension: str,
    db_dir: os.path.abspath,
    tmp_dir: os.path.abspath,
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    boundary_uncertainty: float=0.0,
    qc_method: str="CheckM2",
    completeness: float=0.0,
    contamination: float=100.0,
    dereplicate: bool=False,
    similarity: float=100.0,
    cluster_prefix: str="MSBT",
    logger: Optional[Logger]=None,
    verbose: bool=False,
    nproc: int=1,
    parallel: int=1,
) -> None:
    """Update a database with a new set of metagenome-assembled genomes and reference genomes.
    Also create new clusters in case the input genomes result too far from everything in the database.

    Parameters
    ----------
    input_list : os.path.abspath
        Path to the input file with the list of input genome paths.
    extension : str
        File extension of the input genome files.
    db_dir : os.path.abspath
        Path to the database root folder.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    boundaries : dict
        Dictionary with boundaries produced by the boundaries module.
    boundary_uncertainty : float, default 0.0
        Percentage of kmers to enlarge and reduce boundaries.
    qc_method : {"CheckM2", "CheckV", "EukCC"}, default "CheckM2"
        Quality control method.
    completeness : float, default 0.0
        Threshold on the completeness.
    contamination : float, default 100.0
        Threshold on the contamination.
    dereplicate : bool, default False
        Enable the dereplication step to get rid of replicated genomes.
    similarity : float, default 100.0
        Get rid of genomes according to this threshold in case the dereplication step is enabled.
    cluster_prefix : str, default "MSBT"
        Prefix of clusters numerical identifiers.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.
    nproc : int, default 1
        Make the process parallel when possible.
    parallel : int, default 1
        Used to rebuild multiple clusters in parallel.

    Raises
    ------
    FileNotFoundError
        If the "manifest.txt" file is not available under the database root folder.
    ValueError
        - If the kmer length is 0, the bloom filter size is 0, or the minimum kmer occurrences is <1;
        - If the clusters counter is 0.
    Exception
        - If it is unable to retrieve data from the manifest file;
        - If there are no input genomes available.
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Check whether the manifest file exists in the database
    manifest_filepath = os.path.join(db_dir, "manifest.txt")
    if not os.path.isfile(manifest_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    # Retrieve both the kmer length and the filter size from the manifest file
    manifest = load_manifest(manifest_filepath)

    try:
        kmer_len = manifest["kmer_len"]
        filter_size = manifest["filter_size"]
        clusters_counter = manifest["clusters_counter"]
        min_kmer_occurrences = manifest["min_kmer_occurrences"]
        use_representatives = manifest["use_representatives"]

        # Check whether the kmer length and the filter size have been successfully retrieved
        if kmer_len == 0 or filter_size == 0 or min_kmer_occurrences < 1:
            raise ValueError("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath))

        if clusters_counter == 0:
            raise ValueError("There is nothing to update here!")

    except Exception as ex:
        raise Exception("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath)).with_traceback(
            ex.__traceback__
        )

    in_strains = "strains" if use_representatives else ""

    # Load the list of input genomes and eventually their taxonomic labels
    taxonomy2genomes = load_input_table(input_list, input_extension=extension)

    # Get the list of genome paths
    input_genomes_paths = set().union(*taxonomy2genomes.values())

    # Load the file with the mapping between genome names and taxonomic labels
    # Only in case of input reference genomes
    taxonomies = {
        get_file_info(genome)[1]: taxonomy for taxonomy, genomes in taxonomy2genomes.items() for genome in genomes
    } if "NA" not in taxonomy2genomes else dict()

    input_type = "references" if "NA" not in taxonomy2genomes else "MAGs"

    # Create a temporary folder in case input genomes are gzip compressed
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes")
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # Symlink input genomes
    genomes_paths = list()

    for genome_path in input_genomes_paths:
        genome_link = os.path.join(tmp_genomes_dir, os.path.basename(genome_path))

        if not os.path.islink(genome_link):
            os.symlink(genome_path, genome_link)

        genomes_paths.append(genome_link)

    printline("Processing {} input genomes ({})".format(len(genomes_paths), input_type))

    # Take track of the quality stats
    quality_dict = dict()

    # Check completeness and contamination requirements
    if genomes_paths and (completeness > 0.0 or contamination < 100.0) and qc_method:
        printline(
            "Quality controlling {} genomes (completeness >= {}; contamination <= {})".format(
                len(genomes_paths),
                completeness,
                contamination
            )
        )

        quality_tmp_dir = os.path.join(tmp_dir, "quality")
        os.makedirs(quality_tmp_dir, exist_ok=True)

        # Search for quality output files and try to resume the process
        quality_tables = [str(filepath) for filepath in Path(quality_tmp_dir).glob("run_*.tsv")]

        quality_dict = dict()

        genome_ids = [get_file_info(genome_path)[1] for genome_path in genome_paths]

        if quality_tables:
            # Load the quality tables
            for table_path in quality_tables:
                with open(table_path) as table:
                    header = table.readline().strip()

                    for line in table:
                        line = line.strip()
                        
                        if line:
                            line_split = line.split("\t")

                            genome = line_split[header.index("Name")]

                            if "." in extension:
                                if genome.lower().endswith(".{}".format(extension.lower().split(".")[0])):
                                    # Fix the genome name
                                    genome = ".".join(genome.split(".")[:-1])

                            quality_dict[genome] = dict()

                            for idx, value in enumerate(line_split):
                                quality_dict[genome][header[idx]] = value
        
        if len(quality_dict) < len(genomes_paths):
            # Run the quality control
            # A proper QC method cannot be inferred here because input can be MAGs with no taxonomic label
            quality_dict = quality(
                method=qc_method,
                args={
                    "genomes_paths": genomes_paths,
                    "tmp_dir": quality_tmp_dir,
                    "file_extension": extension,
                    "nproc": nproc
                }
            )

            # Filter genomes according to the input --completeness and --contamination thresholds
            genome_ids = filter_quality(quality_dict, completeness=completeness, contamination=contamination)

        # Build the new list of genome paths that passed the filter
        new_genomes_paths = list()
        for genome_path in genomes_paths:
            _, genome_id, _, _ = get_file_info(genome_path)

            if genome_id in genome_ids:
                new_genomes_paths.append(genome_path)

        printline("{} genomes have been excluded".format(len(genomes_paths) - len(new_genomes_paths)))

        # Define the new list of genomes
        genomes_paths = new_genomes_paths

    # Check whether the input genomes must be dereplicated
    if len(genomes_paths) > 1 and dereplicate:
        printline("Dereplicating {} genomes".format(len(genomes_paths)))
        before_dereplication = len(genomes_paths)

        # Search for the filtered file with the list of excluded genomes
        filtered_filepath = os.path.join(tmp_dir, "howdesbt", "filtered.txt")

        if os.path.isfile(filtered_filepath):
            filtered_genomes_paths = [line.strip() for line in open(filtered_filepath).readlines() if line.strip()]

            # Redefine the list of genomes by removing the filtered ones
            genomes_paths = list(set(genomes_paths).difference(set(filtered_genomes_paths)))

        else:
            genomes_paths = dereplicate_genomes(
                genomes_paths,
                "",
                tmp_dir,
                kmer_len,
                filter_size=filter_size,
                nproc=nproc,
                similarity=similarity,
            )

        printline("{} genomes have been excluded".format(before_dereplication - len(genomes_paths)))

    # Check whether at least one genome survived both the quality control and dereplication steps
    if not genomes_paths:
        raise Exception("No input genomes available!")

    # In case one or more reference genomes are assigned to an unknown cluster
    # Keep track of the unknown taxonomy and change the lineage with
    # the most occurring taxa among the assigned reference genomes
    to_known_taxa: Dict[str, List[str]] = dict()

    # Keep track of the lineages that must be rebuilt
    rebuild: List[str] = list()

    # Keep track of the unassigned genomes
    # They will be assigned to new groups
    unassigned: List[str] = list()

    # Define a partial function around profile_and_assign
    profile_and_assign_partial = partial(
        profile_and_assign,
        input_type=input_type,
        tmp_dir=tmp_dir,
        db_dir=db_dir,
        tmp_genomes_dir=tmp_genomes_dir,
        boundaries=boundaries,
        boundary_uncertainty=boundary_uncertainty,
        taxonomies=taxonomies,
        dereplicate=dereplicate,
        similarity=similarity,
        quality_dict=quality_dict,
        logger=logger if nproc == 1 else None,
        verbose=verbose if nproc == 1 else False,
    )

    with mp.Pool(processes=nproc) as pool, tqdm.tqdm(
        total=len(genomes_paths), disable=(not verbose or nproc == 1)
    ) as pbar:
        # Wrapper around the update function of tqdm
        def progress(*args):
            pbar.update()

        # Process the input genome files
        jobs = [
            pool.apply_async(profile_and_assign_partial, args=(genome_path,), callback=progress)
            for genome_path in genomes_paths
        ]

        # Get results from jobs
        for job in jobs:
            partial_to_known_taxa, partial_rebuild, partial_unassigned = job.get()

            # Merge partial results
            for tx in partial_to_known_taxa:
                if tx not in to_known_taxa:
                    to_known_taxa[tx] = list()
                to_known_taxa[tx].extend(partial_to_known_taxa[tx])

            rebuild.extend(partial_rebuild)
            unassigned.extend(partial_unassigned)

    # In case a reference has been assigned to an unknown cluster
    # Update the cluster taxonomy by applying a majority voting on the reference genomes taxa
    if to_known_taxa:
        printline("Characterising {} unknown taxa".format(len(to_known_taxa)))

        # Iterate over the unknown clusters in dictionary
        for unknown_taxonomy in to_known_taxa:
            # Rebuild the folder path to the unknown taxonomy in the database
            unknown_taxonomy_path = os.path.join(db_dir, unknown_taxonomy.replace("|", os.sep))

            # Take track of the reference genomes taxonomic labels
            known_taxa = list()

            # Retrieve the reference genomes taxa from the input mapping file
            # and finally move the genome files
            for genome_path in to_known_taxa[unknown_taxonomy]:
                _, genome_name, genome_ext, compression = get_file_info(genome_path)
                known_taxa.append(taxonomies[genome_name])

                if not os.path.isfile(os.path.join(unknown_taxonomy_path, in_strains, "genomes", os.path.basename(genome_path))):
                    # Copy the input genome file to the genomes folder
                    shutil.copy(genome_path, os.path.join(unknown_taxonomy_path, in_strains, "genomes"))

                # In case the input genome is not Gzip compressed
                genome_filepath = os.path.join(unknown_taxonomy_path, in_strains, "genomes", "{}{}".format(genome_name, genome_ext))
                if not compression and not os.path.isfile("{}{}".format(genome_filepath, compression)):
                    run(["gzip", genome_filepath], silence=True)

                references_entries = [line.strip() for line in open(os.path.join(unknown_taxonomy_path, "references.txt")).readlines() if line.strip()]

                if not genome_name in references_entries:
                    # Add the current genome to the list of reference genomes
                    with open(os.path.join(unknown_taxonomy_path, "references.txt"), "a+") as references_file:
                        references_file.write("{}\n".format(genome_name))

            # Get the most occurring taxonomic label
            known_taxa_counter: List[Tuple[str, int]] = Counter(known_taxa).most_common()  # type: ignore

            # In case of multiple most occurring taxa, pick the first one in lexicographic order
            most_common_taxa = list(
                takewhile(lambda x: x[1] >= known_taxa_counter[0][1], known_taxa_counter)  # noqa: B023
            )

            assigned_taxonomy = sorted(most_common_taxa, key=lambda x: x[0])[0][0]

            # Split taxonomies into levels
            unknown_levels: List[str] = unknown_taxonomy.split("|")
            assigned_levels: List[str] = assigned_taxonomy.split("|")

            # Characterise unknown clusters
            # Iterate backwards over levels
            for i in range(6, -1, -1):
                if unknown_levels[i] == assigned_levels[i]:
                    break

                else:
                    # Get the partial levels for the new taxonomy
                    new_levels = assigned_levels[: i + 1]
                    new_levels_subpath = os.path.join(db_dir, os.sep.join(new_levels))
                    os.makedirs(new_levels_subpath, exist_ok=True)

                    # Get the partial levels for the old taxonomy
                    old_levels = unknown_levels[: i + 1]
                    old_levels_subpath = os.path.join(db_dir, os.sep.join(old_levels))

                    # Start moving folders
                    for folder_path in Path(old_levels_subpath).glob("**/*"):
                        if os.path.isdir(str(folder_path)):
                            shutil.move(str(folder_path), new_levels_subpath)

                    # Also move remaining files
                    for file_path in Path(old_levels_subpath).glob("**/*"):
                        if os.path.isfile(str(file_path)):
                            shutil.move(str(file_path), new_levels_subpath)

                    # Remove the old bloom filter root node
                    bloom_filter_node = os.path.join(new_levels_subpath, "{}.bf".format(unknown_levels[i]))

                    if os.path.isfile(bloom_filter_node):
                        os.unlink(bloom_filter_node)

            # Add the new renamed cluster to the rebuild list of taxa
            rebuild.append(assigned_taxonomy)

    # Cluster unassigned genomes before rebuilding the updated lineages
    if unassigned:
        printline("Defining new clusters for {} unassigned genomes".format(len(unassigned)))

        assignments_filepath = os.path.join(db_dir, "assignments.txt")
        unassigned_filepath = os.path.join(db_dir, "unassigned.txt")

        if os.path.isfile(assignments_filepath) and os.path.isfile(unassigned_filepath):
            # Load the assignments and resume the process
            assignments = dict()

            not_assigned = list()

            with open(assignments_filepath) as assignments_file:
                for line in assignments_file:
                    line = line.strip()
                    if line:
                        if not line.startswith("#"):
                            line_split = line.split("\t")
                            assignments[line_split[0]] = {
                                "taxonomy": line_split[1],
                                "cluster": line_split[2]
                            }
            
            # Fix assignments keys with genome file paths
            for genome_filepath in unassigned:
                _, genome_name, _, _ = get_file_info(genome_filepath)

                if genome_name in assignments:
                    assignments[genome_filepath] = assignments.pop(genome_name)

                else:
                    not_assigned.append(genome_filepath)
        
        else:
            if os.path.isfile(assignments_filepath):
                os.unlink(assignments_filepath)

            if os.path.isfile(unassigned_filepath):
                os.unlink(unassigned_filepath)

            # Define a temporary folder for building bloom filters
            # This already exists in case the dereplication step has been executed
            howdesbt_tmp_dir = os.path.join(tmp_dir, "howdesbt_unassigned")

            # Cluster genomes according to the boundaries defined by the boundaries module
            # Define a cluster for each taxomomic level
            # Look at the genomes profiles and update or build new clusters
            assignments, not_assigned = cluster(
                unassigned,
                boundaries,
                manifest_filepath,
                os.path.join(tmp_dir, "profiling"),
                howdesbt_tmp_dir,
                os.path.join(db_dir, "assignments.txt"),
                cluster_prefix=cluster_prefix,
                min_occurrences=min_kmer_occurrences,
                nproc=nproc,
            )

            if not_assigned:
                printline("Unable to process {} genomes".format(len(not_assigned)))

                with open(os.path.join(db_dir, "unassigned.txt"), "w+") as unassigned_file:
                    for genome_path in not_assigned:
                        unassigned_file.write("{}\n".format(get_file_info(genome_path)[1]))

        # Iterate over the input genomes
        for genome_path in assignments:
            _, genome_name, genome_ext, compressed = get_file_info(genome_path)

            # Create the new cluster folder in the database
            tax_dir = os.path.join(db_dir, assignments[genome_path]["taxonomy"].replace("|", os.sep))
            os.makedirs(os.path.join(tax_dir, "genomes"), exist_ok=True)

            # Also create the strains folder
            os.makedirs(os.path.join(tax_dir, in_strains, "genomes"), exist_ok=True)

            # Copy the input genome into the genomes folder of the new cluster
            if not os.path.isfile(os.path.join(tax_dir, in_strains, "genomes", os.path.basename(genome_path))):
                shutil.copy(genome_path, os.path.join(tax_dir, in_strains, "genomes"))

            # In case the input genome is not gzip compressed
            if not compressed:
                genome_filepath = os.path.join(tax_dir, in_strains, "genomes", "{}{}".format(genome_name, genome_ext))
                if not os.path.isfile("{}{}".format(genome_filepath, compressed)):
                    run(["gzip", genome_filepath], silence=True)

            mags_entries = list()

            if os.path.isfile(os.path.join(tax_dir, "mags.txt")):
                mags_entries = [line.strip() for line in open(os.path.join(tax_dir, "mags.txt")).readlines() if line.strip()]

            if not genome_name in mags_entries:
                # Also update the mags.txt file
                with open(os.path.join(tax_dir, "mags.txt"), "a+") as mags_file:
                    mags_file.write("{}\n".format(genome_name))

            # Also report the quality stats of the genomes in the new clusters
            if genome_name in quality_dict:
                quality_filepath = os.path.join(tax_dir, "quality.tsv")
                add_to_quality_table = True

                if os.path.isfile(quality_filepath):
                    genomes_in_quality_file = list()

                    with open(quality_filepath) as table:
                        header = table.readline().strip()
                        for line in table:
                            line = line.strip()
                            if line:
                                line_split = line.split("\t")
                                genomes_in_quality_file.append(line_split[header.index("Name")])

                    if genome_name in genomes_in_quality_file:
                        add_to_quality_table = False

                if add_to_quality_table:
                    not_exists = not os.path.isfile(quality_filepath)

                    with open(quality_filepath, "a+") as quality_file:
                        header = sorted(list(quality_dict[genome_name].keys()))

                        if not_exists:                                
                            quality_file.write("{}\n".format("\t".join(header)))

                        quality_file.write("{}\n".format("\t".join([quality_dict[genome_name][h] for h in header])))

            # Also initialize the metadata table with the cluster ID
            metadata_filepath = os.path.join(tax_dir, "metadata.tsv")

            if not os.path.isfile(metadata_filepath):
                with open(metadata_filepath, "w+") as metadata_file:
                    metadata_file.write("# Cluster ID: {}".format(assignments[genome_path]["cluster"]))

            # Add the full taxonomy to the list of taxonomic labels that must be rebuilt
            rebuild.append(assignments[genome_path]["taxonomy"])

    # Check whether there is at least one lineage that must be rebuilt
    rebuild = list(set(rebuild))

    if not rebuild:
        printline("No lineages have been updated")

    else:
        printline("Updating the database")

        level_names = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        level_ids = [name[0] for name in level_names]

        # Extract the taxonomic levels from the list of taxa that must be rebuilt
        # Process all the species first, then all the genera, and so on up to the kingdom level
        for i in range(6, -1, -1):
            taxalist = set()
            for label in rebuild:
                levels = label.split("|")
                # Temporarily skip partial taxonomic labels if the number of levels is lower than i
                # Waiting for the right i
                if len(levels) < i + 1:
                    continue

                # Build the partial taxonomic label
                taxalist.add("|".join(levels[: i + 1]))

            # Get the level ID of the last characterized taxonomic level
            level_id = list(taxalist)[0].split("|")[-1][0]

            printline("Rebuilding {}".format(level_names[level_ids.index(level_id)]))

            build_cluster_partial = partial(
                build_cluster,
                db_dir=db_dir,
                tmp_dir=tmp_dir,
                kmer_len=kmer_len,
                filter_size=filter_size,
                min_kmer_occurrences=min_kmer_occurrences,
                nproc=nproc,
                use_representatives=use_representatives,
                logger=logger if parallel == 1 else None,
                verbose=verbose if parallel == 1 else False,
            )

            with mp.Pool(processes=parallel) as pool, tqdm.tqdm(
                total=len(taxalist), disable=(not verbose or parallel == 1)
            ) as pbar:
                # Wrapper around the update function of tqdm
                def progress(*args):
                    pbar.update()

                # Process the input genome files
                jobs = [
                    pool.apply_async(build_cluster_partial, args=(taxonomy,), callback=progress)
                    for taxonomy in taxalist
                ]

                # Get results from jobs
                for job in jobs:
                    taxonomy, remove_from_rebuild = job.get()

                    if remove_from_rebuild:
                        rebuild.remove(taxonomy)


def main() -> None:
    # Load command line parameters
    args = read_params(sys.argv[1:])

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether this run is a "resume"
    resume = args.resume

    # Take track of the framework version to check for code compatibility
    parent_version = args.parent_version

    if args.resume:
        # Try to resume the update process
        # Read the last command line from the sh script
        argv = [line.strip() for line in open(args.resume).readlines() if line.strip() and not line.startswith("#")][-1].split()

        # Reload command line parameters
        if argv[0] == TOOL_ID:
            # It was a direct call to the update module
            args = read_params(argv[1:])

        else:
            # The update module has been run through the controller
            argv = read_params(argv[2:])

        if args.parent_version != parent_version:
            println(
                "Warning: your MetaSBT version is v{} while the process you are trying to resume has been initially run with v{}".format(
                    parent_version,
                    args.parent_version
                ),
                logger=logger,
                verbose=args.verbose,
            )

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Check whether the database folder exists
    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    # Check whether the file with the list of input genome paths exists
    if not os.path.isfile(args.input_list):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.input_list)

    # Check whether the boundaries table exists
    if not os.path.isfile(args.boundaries):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.boundaries)
    
    # Load the boundaries table
    boundaries_table = load_boundaries(args.boundaries)

    if not resume:
        # Build a sh script with the command line used to launch the update module
        build_sh(sys.argv, TOOL_ID, args.db_dir)

    t0 = time.time()

    update(
        args.input_list,
        args.extension,
        args.db_dir,
        args.tmp_dir,
        boundaries_table,
        boundary_uncertainty=args.boundary_uncertainty,
        qc_method=args.quality_control,
        completeness=args.completeness,
        contamination=args.contamination,
        dereplicate=args.dereplicate,
        similarity=args.similarity,
        cluster_prefix=args.cluster_prefix,
        logger=logger,
        verbose=args.verbose,
        nproc=args.nproc,
        parallel=args.parallel,
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
