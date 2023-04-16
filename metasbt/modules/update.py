#!/usr/bin/env python3
"""
Update a specific database with a new set of reference genomes or metagenome-assembled genomes
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 15, 2023"

import argparse as ap
import errno
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
        checkm,
        cluster,
        dereplicate_genomes,
        filter_checkm_tables,
        get_file_info,
        get_level_boundaries,
        howdesbt,
        init_logger,
        load_boundaries,
        load_manifest,
        number,
        println,
        run,
    )
except Exception:
    pass

# Define the module name
TOOL_ID = "update"

# Define the list of dependencies
DEPENDENCIES = ["checkm", "howdesbt"]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Database folder path
    "--input-list",  # Input file path
    "--log",  # Path to the log file
    "--tmp-dir",  # Temporary folder path
]

# Define the software root directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
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
        required=True,
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
        required=True,
        dest="db_dir",
        help="This is the database directory with the taxonomically organised sequence bloom trees",
    )
    general_group.add_argument(
        "--extension",
        type=str,
        required=True,
        choices=["fa", "fa.gz", "fasta", "fasta.gz", "fna", "fna.gz"],
        help=(
            "Specify the input genome files extension. "
            "All the input genomes must have the same file extension before running this module"
        ),
    )
    general_group.add_argument(
        "--fast-profile",
        action="store_true",
        default=False,
        dest="fast_profile",
        help=(
            "Query the tree with the representative genomes at the species level "
            "instead of the strains tree (fast but low precision in establishing the closest genome)"
        ),
    )
    general_group.add_argument(
        "--input-list",
        type=os.path.abspath,
        required=True,
        dest="input_list",
        help="This file contains the list of paths to the new genomes that will be added to the database",
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
        help="This argument refers to the number of processors used for parallelizing the pipeline when possible",
    )
    general_group.add_argument(
        "--parallel",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        help="Maximum number of processors to process input genomes in parallel",
    )
    general_group.add_argument(
        "--taxa",
        type=os.path.abspath,
        help=(
            "Input file with the mapping between input genome IDs and their taxonomic label. "
            'This is used in case of reference genomes only "--type references"'
        ),
    )
    general_group.add_argument(
        "--tmp-dir",
        type=os.path.abspath,
        required=True,
        dest="tmp_dir",
        help="Path to the folder for storing temporary data",
    )
    general_group.add_argument(
        "--type",
        type=str,
        required=True,
        choices=["MAGs", "references"],
        help="Define the nature of the input genomes",
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

    # Group of arguments for CheckM
    qc_group = p.add_argument_group("CheckM: Quality control")

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
    qc_group.add_argument(
        "--pplacer-threads",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        dest="pplacer_threads",
        help="Maximum number of threads for pplacer. This is required to maximise the CheckM performances",
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

    return p.parse_args()


def load_taxa(taxa_filepath: Optional[str] = None, genomes_path: Optional[List[str]] = None) -> Dict[str, str]:
    """
    Load the taxa file with the mapping between the input genome names and their taxonomic labels

    :param taxa_filepath:   Path to the input mapping file with the taxonomic labels
    :param genomes_path:    List with paths to the input genome files
    :return:                Dictionary with mapping between genome names and their taxonomic labels
    """

    if taxa_filepath is None or genomes_path is None:
        return dict()

    if not os.path.isfile(taxa_filepath) or len(genomes_path) == 0:
        return dict()

    # Define a mapping between input genome names and taxonomic labels
    taxonomies = dict()

    # Load taxonomic labels
    with open(taxa_filepath) as taxa_file:
        for line in taxa_file:
            line = line.strip()
            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")
                    taxonomies[line_split[0]] = line_split[1]

    # Retrieve the input genome names
    genome_names = list()

    if genomes_path:
        for genome_path in genomes_path:
            # Get the genome name
            _, genome_name, _, _ = get_file_info(genome_path)

            genome_names.append(genome_name)

    # Check whether a taxonomic label is defined for each of the input genome
    if len(set(genome_names).intersection(set(taxonomies.keys()))) < len(set(genome_names)):
        raise Exception("Unable to find all the taxonomic labels for the input reference genomes")

    return taxonomies


def profile_and_assign(
    genome_path: str,
    input_type: str,
    tmp_dir: str,
    db_dir: str,
    tmp_genomes_dir: str,
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    boundary_uncertainty: float = 0.0,
    fast_profile: bool = False,
    taxonomies: Optional[dict] = None,
    dereplicate: bool = False,
    similarity: float = 100.0,
    checkm_header: Optional[str] = None,
    checkm_data: Optional[dict] = None,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> Tuple[dict, list, list]:
    """
    Profile and assign one input genome

    :param genome_path:             Path to the input genome file
    :param input_type:              Nature of the input genomes (MAGs or references)
    :param tmp_dir:                 Path to the temporary folder
    :param db_dir:                  Path to the database root folder
    :param tmp_genomes_dir:         Path to the temporary folder for uncompressing input genomes
    :param boundaries:              Boundaries table produced by the boundaries module
    :param boundary_uncertainty:    Percentage of kmers to enlarge and reduce boundaries
    :param fast_profile:            Query the tree with the representative genomes under the species level instead
                                    of the strains tree (fast but low precision in establishing the closest genome)
    :param taxonomies:              Dictionary with mapping between input genome names and their taxonomic labels (for reference genomes only)
    :param dereplicate:             Enable dereplication of input genome against genomes in the closest cluster
    :param similarity:              Exclude input genome if its similarity with the closest genome in the database exceed this threshold
    :param checkm_data:             Dictionary with CheckM statistics
    :param logger:                  Logger object
    :param verbose:                 Print messages on screen
    :return:                        The assignment
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    if taxonomies is None:
        taxonomies = dict()

    if checkm_data is None:
        checkm_data = dict()

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

        # Unzip the genome file into the tmp folder
        with open(filepath, "w+") as file:
            run(["gunzip", "-c", genome_path], stdout=file)

    printline("Profiling {}".format(genome_name))

    # Run the profiler to establish the closest genome and the closest group
    # for each taxonomic level in the tree
    # The profile module is in the same folder of the update module
    profiler = [
        sys.executable,
        os.path.join(SCRIPT_DIR, "profile.py"),
        "--input-file",
        filepath,
        "--input-id",
        filepath,
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

    if fast_profile:
        profiler.append("--fast")

    run(profiler, silence=True)

    # Define the path to the profiler output file
    profile_path = os.path.join(tmp_dir, "profiling", "{}__profiles.tsv".format(genome_name))
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
        closest_taxa = "|".join(profile_data["genome"]["taxonomy"].split("|")[:-1])
        closest_common_kmers = profile_data["genome"]["common_kmers"]
        closest_score = profile_data["genome"]["score"]

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
                    # Assign the current genome to the closest lineage
                    shutil.copy(genome_path, os.path.join(closest_taxadir, "strains", "genomes"))

                    if not compression:
                        # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                        run(["gzip", os.path.join(closest_taxadir, "strains", "genomes", "{}{}".format(genome_name, genome_ext))], silence=True)

                    # Add the current genome to the list of MAGs
                    with open(mags_filepath, "a+") as file:
                        file.write("{}\n".format(genome_name))

                    # Also add its CheckM statistics if available
                    if genome_name in checkm_data:
                        checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")
                        with open(checkm_filepath, "a+") as checkm_file:
                            if not os.path.isfile(checkm_filepath):
                                checkm_file.write("{}\n".format(checkm_header))
                            checkm_file.write("{}\n".format(checkm_data[genome_name]))

                    # Do not remove the index here because there could be other input genome
                    # with the same closest taxonomic label
                    rebuild.append(closest_taxa)

                    printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                else:
                    # Mark the input genome as unassigned
                    # The CheckM statistics for this genome will be reported after the assignment
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
                        shutil.copy(genome_path, os.path.join(closest_taxadir, "strains", "genomes"))

                        if not compression:
                            # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                            run(["gzip", os.path.join(closest_taxadir, "strains", "genomes", "{}{}".format(genome_name, genome_ext))], silence=True)

                        # Add the current genome to the list of reference genomes
                        with open(references_filepath, "a+") as file:
                            file.write("{}\n".format(genome_name))

                        # Do not remove the index here because there could be other input genome
                        # with the same closest taxonomic label
                        rebuild.append(closest_taxa)

                        printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                    # Also add its CheckM statistics if available
                    if genome_name in checkm_data:
                        checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")

                        with open(checkm_filepath, "a+") as checkm_file:
                            if not os.path.isfile(checkm_filepath):
                                checkm_file.write("{}\n".format(checkm_header))

                            checkm_file.write("{}\n".format(checkm_data[genome_name]))

                else:
                    # If nothing is closed enough to the current genome and its taxonomic label does not exist in the database
                    # Create a new branch with the new taxonomy for the current genome
                    taxdir = os.path.join(db_dir, taxalabel.replace("|", os.sep))
                    os.makedirs(os.path.join(taxdir, "genomes"), exist_ok=True)

                    # Also create the strains folder
                    os.makedirs(os.path.join(taxdir, "strains", "genomes"), exist_ok=True)

                    # Assign the current genome to the closest lineage
                    shutil.copy(genome_path, os.path.join(taxdir, "strains", "genomes"))

                    if not compression:
                        # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                        run(["gzip", os.path.join(taxdir, "strains", "genomes", "{}{}".format(genome_name, genome_ext))], silence=True)

                    # Add the current genome to the list of reference genomes
                    with open(os.path.join(taxdir, "references.txt"), "a+") as file:
                        file.write("{}\n".format(genome_name))

                    # Do not remove the index here because there could be other input genome
                    # with the same closest taxonomic label
                    rebuild.append(closest_taxa)

                    printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                    # Also add its CheckM statistics if available
                    if genome_name in checkm_data:
                        checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")

                        with open(checkm_filepath, "a+") as checkm_file:
                            if not os.path.isfile(checkm_filepath):
                                checkm_file.write("{}\n".format(checkm_header))

                            checkm_file.write("{}\n".format(checkm_data[genome_name]))

    # Remove the uncompressed version of the input genome in the temporary folder
    if os.path.isfile(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0])):
        os.unlink(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0]))

    return to_known_taxa, rebuild, unassigned


def update(
    input_list: str,
    input_type: str,
    extension: str,
    db_dir: str,
    tmp_dir: str,
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    boundary_uncertainty: float = 0.0,
    taxa_map: Optional[str] = None,
    fast_profile: bool = False,
    completeness: float = 0.0,
    contamination: float = 100.0,
    dereplicate: bool = False,
    similarity: float = 100.0,
    cluster_prefix: str = "MSBT",
    logger: Optional[Logger] = None,
    verbose: bool = False,
    nproc: int = 1,
    pplacer_threads: int = 1,
    parallel: int = 1,
) -> None:
    """
    Update a database with a new set of metagenome-assembled genomes and reference genomes.
    Also create new clusters in case the input genomes result too far from everything in the database

    :param input_list:              Path to the input file with the list of input genome paths
    :param input_type:              Nature of the input genomes (MAGs or references)
    :param extension:               File extension of the input genome files
    :param db_dir:                  Path to the database root folder
    :param tmp_dir:                 Path to the temporary folder
    :param boundaries:              Boundaries table produced by the boundaries module
    :param boundary_uncertainty:    Percentage of kmers to enlarge and reduce boundaries
    :param taxa_map:                Path to the file with the mapping between the input genome names and their taxonomic labels
                                    Used only in case of input reference genomes
    :param fast_profile:            Query the tree with the representative genomes under the species level instead
                                    of the strains tree (fast but low precision in establishing the closest genome)
    :param completeness:            Threshold on the CheckM completeness
    :param contamination:           Threshold on the CheckM contamination
    :param dereplicate:             Enable the dereplication step to get rid of replicated genomes
    :param similarity:              Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param cluster_prefix:          Prefix of clusters numerical identifiers
    :param logger:                  Logger object
    :param verbose:                 Print messages on screen
    :param nproc:                   Make the process parallel when possible
    :param pplacer_threads:         Maximum number of threads to make pplacer parallel with CheckM
    :param parallel:                Maximum number of processors to process input genomes in parallel
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

        # Check whether the kmer length and the filter size have been successfully retrieved
        if kmer_len == 0 or filter_size == 0 or min_kmer_occurrences < 2:
            raise ValueError("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath))

        if clusters_counter == 0:
            raise ValueError("There is nothing to update here!")

    except Exception as ex:
        raise Exception("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath)).with_traceback(
            ex.__traceback__
        )

    # Load the list of genome paths
    input_genomes_paths = [path.strip() for path in open(input_list).readlines() if path.strip()]

    # Create a temporary folder in case input genomes are gzip compressed
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes")
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # Symlink input genomes
    genomes_paths = list()
    for genome_path in input_genomes_paths:
        os.symlink(genome_path, os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))
        genomes_paths.append(os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))

    # Load the file with the mapping between genome names and taxonomic labels
    # Only in case of input reference genomes
    taxonomies: Dict[str, str] = load_taxa(taxa_map, genomes_paths) if input_type == "references" else dict()

    printline("Processing {} input genomes ({})".format(len(genomes_paths), input_type))

    # Take track of the CheckM data
    # Mapping between genome names and lines in the CheckM output tables
    checkm_header = ""
    checkm_data = dict()

    # Check completeness and contamination requirements
    if genomes_paths and (completeness > 0.0 or contamination < 100.0):
        printline("Quality controlling {} genomes".format(len(genomes_paths)))
        printline("Minimum completeness threshold: {}".format(completeness))
        printline("Maximum contamination threshold: {}".format(contamination))

        checkm_tmp_dir = os.path.join(tmp_dir, "checkm")
        os.makedirs(checkm_tmp_dir, exist_ok=True)

        # Run CheckM to quality control genomes
        checkm_tables = checkm(
            genomes_paths,
            checkm_tmp_dir,
            file_extension=extension,
            nproc=nproc,
            pplacer_threads=pplacer_threads,
        )

        # Load CheckM tables
        for checkm_table in checkm_tables:
            with open(checkm_table) as table:
                first_line = True
                for line in table:
                    line = line.strip()
                    if line:
                        if first_line:
                            checkm_header = line
                            first_line = False
                        else:
                            line_split = line.split("\t")
                            checkm_data[line_split[0]] = line

        # Filter genomes according to the input --completeness and --contamination thresholds
        genome_ids = filter_checkm_tables(checkm_tables, completeness=0.0, contamination=100.0)

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
        fast_profile=fast_profile,
        taxonomies=taxonomies,
        dereplicate=dereplicate,
        similarity=similarity,
        checkm_header=checkm_header,
        checkm_data=checkm_data,
        logger=logger if parallel == 1 else None,
        verbose=verbose if parallel == 1 else False,
    )

    with mp.Pool(processes=parallel) as pool, tqdm.tqdm(
        total=len(genomes_paths), disable=(not verbose or parallel == 1)
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
                _, genome_name, _, compression = get_file_info(genome_path)
                known_taxa.append(taxonomies[genome_name])

                # Copy the input genome file to the genomes folder
                shutil.copy(genome_path, os.path.join(unknown_taxonomy_path, "strains", "genomes"))

                # In case the input genome is not Gzip compressed
                if not compression:
                    gzip_genome_filepath = os.path.join(unknown_taxonomy_path, "strains", "genomes", os.path.basename(genome_path))
                    run(["gzip", gzip_genome_filepath], silence=True)

                # Also update the list of reference genomes
                with open(os.path.join(unknown_taxonomy_path, "references.txt"), "a+") as file:
                    file.write("{}\n".format(genome_name))

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
            _, genome_name, _, compressed = get_file_info(genome_path)

            # Create the new cluster folder in the database
            tax_dir = os.path.join(db_dir, assignments[genome_path]["taxonomy"].replace("|", os.sep))
            os.makedirs(os.path.join(tax_dir, "genomes"), exist_ok=True)

            # Also create the strains folder
            os.makedirs(os.path.join(tax_dir, "strains", "genomes"), exist_ok=True)

            # Copy the input genome into the genomes folder of the new cluster
            shutil.copy(genome_path, os.path.join(tax_dir, "strains", "genomes"))

            # In case the input genome is not gzip compressed
            if not compressed:
                gzip_genome_filepath = os.path.join(tax_dir, "strains", "genomes", os.path.basename(genome_path))
                run(["gzip", gzip_genome_filepath], silence=True)

            # Also update the mags.txt file
            with open(os.path.join(tax_dir, "mags.txt"), "a+") as file:
                file.write("{}\n".format(genome_name))

            # Also report the CheckM statistics of the genomes in the new clusters
            if genome_name in checkm_data:
                checkm_filepath = os.path.join(tax_dir, "checkm.tsv")

                with open(checkm_filepath, "a+") as checkm_file:
                    # Add the header line in case the CheckM table does not exist
                    if not os.path.isfile(checkm_filepath):
                        checkm_file.write("{}\n".format(checkm_header))

                    checkm_file.write("{}\n".format(checkm_data[genome_name]))
            
            # Also initialize the metadata table with the cluster ID
            metadata_filepath = os.path.join(tax_dir, "metadata.tsv")
            if not os.path.isfile(metadata_filepath):
                with open(metadata_filepath, "w+") as metadata_file:
                    metadata_file.write("# Cluster ID: {}".format(assignments[genome_name]["cluster"]))

            # Add the full taxonomy to the list of taxonomic labels that must be rebuilt
            rebuild.append(assignments[genome_name]["taxonomy"])

    # Check whether there is at least one lineage that must be rebuilt
    rebuild = list(set(rebuild))

    if not rebuild:
        printline("No lineages have been updated")

    else:
        printline("Updating the database")
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

            # Rebuild the sequence bloom trees
            for taxonomy in taxalist:
                build_cluster = True

                printline("\t{}".format(taxonomy))
                tax_dir = os.path.join(db_dir, taxonomy.replace("|", os.sep))

                if taxonomy.split("|")[-1].startswith("s__"):
                    # In case the taxonomy is complete at the species level
                    # Load manifest
                    strains_manifest = load_manifest(os.path.join(tax_dir, "strains", "manifest.txt"))
                    
                    # Remove the old index
                    if os.path.isdir(os.path.join(tax_dir, "strains", "index")):
                        shutil.rmtree(os.path.join(tax_dir, "strains", "index"), ignore_errors=True)

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
                    bf_filepaths = [str(path) for path in Path(os.path.join(tax_dir, "strains", "filters")).glob("*.bf.gz")]

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
                    bfdistance_sums = {genome: sum(bfdistance_theta[genome].values()) for genome in bfdistance_theta}

                    # Sort genomes according to the sum of their distances with all the other genomes
                    sorted_genomes = sorted(bfdistance_sums, key=lambda genome: bfdistance_sums[genome])

                    # First and last genomes are those that minimize and maximize the distance with all the other genomes
                    selected_genomes.append(sorted_genomes[0])
                    selected_genomes.append(sorted_genomes[-1])

                    # Also select a genome in the middles of min and max distances
                    selected_genomes.append(sorted_genomes[int(len(sorted_genomes) / 2)])

                    # Get the list of previously selected representative genomes
                    previously_selected_genomes = [
                        get_file_info(str(path))[1] for path in Path(os.path.join(tax_dir, "filters").glob("*.bf.gz"))
                    ]

                    add_selected_genomes = set(selected_genomes).difference(set(previously_selected_genomes))
                    remove_selected_genomes = set(previously_selected_genomes).difference(set(selected_genomes))

                    if len(remove_selected_genomes) > 0:
                        # In case the representative genomes of the species did not change
                        # Remove the taxonomy from the rebuild list
                        rebuild.remove(taxonomy)
                        build_cluster = False

                        # Remove old representatives
                        for genome_path in Path(os.path.join(tax_dir, "genomes")).glob("*.gz"):
                            _, genome_name, _, _ = get_file_info(str(genome_path))

                            if genome_name in remove_selected_genomes:
                                os.unlink(str(genome_path))

                                os.unlink(os.path.join(tax_dir, "filters", "{}.bf.gz".format(genome_name)))
                        
                        # Add new representatives
                        for genome_path in Path(os.path.join(tax_dir, "strains", "genomes")).glob("*.gz"):
                            _, genome_name, _, _ = get_file_info(str(genome_path))

                            if genome_name in add_selected_genomes:
                                os.symlink(str(genome_path), os.path.join(tax_dir, "genomes", os.path.basename(str(genome_path))))

                if build_cluster:
                    # Remove the old index if it exists
                    if os.path.isdir(os.path.join(tax_dir, "index")):
                        shutil.rmtree(os.path.join(tax_dir, "index"), ignore_errors=True)

                    # Remove filters
                    if os.path.isdir(os.path.join(tax_dir, "filters")):
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


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether the database folder exists
    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    # Check whether the file with the list of input genome paths exists
    if not os.path.isfile(args.input_list):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.input_list)

    # Check whether the input file with the mapping between genome name and taxonomic labels exists
    # Only in case of input reference genomes
    if args.type == "references" and not os.path.isfile(args.taxa):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.taxa)

    # Check whether the boundaries table exists
    if not os.path.isfile(args.boundaries):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.boundaries)
    
    # Load the boundaries table
    boundaries_table = load_boundaries(args.boundaries)

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Build a sh script with the command line used to launch the update module
    build_sh(sys.argv, TOOL_ID, args.db_dir)

    t0 = time.time()

    update(
        args.input_list,
        args.type,
        args.extension,
        args.db_dir,
        args.tmp_dir,
        boundaries_table,
        boundary_uncertainty=args.boundary_uncertainty,
        taxa_map=args.taxa,
        fast_profile=args.fast_profile,
        completeness=args.completeness,
        contamination=args.contamination,
        dereplicate=args.dereplicate,
        similarity=args.similarity,
        cluster_prefix=args.cluster_prefix,
        logger=logger,
        verbose=args.verbose,
        nproc=args.nproc,
        pplacer_threads=args.pplacer_threads,
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
