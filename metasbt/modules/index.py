#!/usr/bin/env python3
"""
Build a database with a set of genomes indexed with HowDeSBT.
Genomes are provided as inputs or automatically downloaded from NCBI GenBank
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 1, 2023"

import argparse as ap
import copy
import gzip
import hashlib
import math
import multiprocessing as mp
import os
import re
import shutil
import sys
import tarfile
import time
from functools import partial
from logging import Logger
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy  # type: ignore
import tqdm  # type: ignore

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import (  # type: ignore  # isort: skip
        build_sh,
        checkm,
        dereplicate_genomes,
        download,
        estimate_bf_size,
        filter_checkm_tables,
        get_file_info,
        howdesbt,
        init_logger,
        integrity_check,
        number,
        optimal_k,
        println,
        run,
        validate_url,
    )
except Exception:
    pass

# Define the module name
TOOL_ID = "index"

# Define the list of dependencies
DEPENDENCIES = [
    "checkm",
    "howdesbt",
    "kitsune",
    "ncbitax2lin",
    "ntcard",
    "wget",
]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Database folder path
    "--input-list",  # File with a list of paths to the input genomes
    "--log",  # Path to the log file
    "--tmp-dir",  # Temporary folder path
]

# Define the url to the NCBI taxdump
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Define the url to the NCBI GenBank Assembly Summary
# https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
# https://www.ncbi.nlm.nih.gov/assembly/help/
ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

# Consider a genome as a reference if it contains one of the following tags
# under the excluded_from_refseq column in the NCBI GenBank Assembly Summary table
# https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
REFERENCE_TAGS = [
    "derived from single cell",
    "derived from surveillance project",
    "assembly from type material",
    "assembly from synonym type material",
    "assembly designated as neotype",
    "assembly designated as reftype",
    "assembly from pathotype material",
    "assembly from proxytype material",
    "missing strain identifier",
    "genus undefined",
    "from large multi-isolate project"
]


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "Build a database with a set of genomes indexed with HowDeSBT. "
            "Genomes are provided as inputs or automatically downloaded from NCBI GenBank"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )

    # General arguments
    general_group = p.add_argument_group("General arguments")

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
        "--flat-structure",
        action="store_true",
        default=False,
        dest="flat_structure",
        help=(
            "Organize genomes without any taxonomic organization. "
            "This will lead to the creation of a single sequence bloom tree"
        ),
    )
    general_group.add_argument(
        "--input-list",
        type=os.path.abspath,
        dest="input_list",
        help=(
            "Path to the input table with a list of genome file paths and an optional column with their taxonomic labels. "
            "Please note that the input genome files must be gz compressed with fna extension (i.e.: *.fna.gz)"
        ),
    )
    general_group.add_argument(
        "--kingdom",
        type=str,
        help=(
            "Consider genomes whose lineage belongs to a specific kingdom. "
            "It is optional and must be provided in conjunction with --superkingdom"
        ),
    )
    general_group.add_argument(
        "--kmer-len",
        type=number(int, minv=4),
        dest="kmer_len",
        help="This is the length of the kmers used for building bloom filters",
    )
    general_group.add_argument(
        "--limit-estimation-number",
        type=number(int, minv=1),
        dest="limit_estimation_number",
        help=(
            "Limit the number of genomes per group to be considered as input for kitsune and ntCard. "
            "Must be used in conjunction with --estimate-kmer-size and/or --estimate-filter-size. "
            "It overrides --limit-estimation-percentage in case of a number > 0"
        ),
    )
    general_group.add_argument(
        "--limit-estimation-percentage",
        type=number(float, minv=sys.float_info.min, maxv=100.0),
        default=100.0,
        dest="limit_estimation_percentage",
        help=(
            "Percentage on the total number of genomes per group to be considered as input for kitsune and ntCard. "
            "Must be used in conjunction with --estimate-kmer-size and/or --estimate-filter-size"
        ),
    )
    general_group.add_argument(
        "--limit-genomes",
        type=number(int, minv=1),
        default=numpy.Inf,
        dest="limit_genomes",
        help=(
            "Limit the number of genomes per species. "
            "This will remove the exceeding number of genomes randomly to cut the overall number of genomes per species to this number. "
            "The number of genomes per species is not limited by default"
        ),
    )
    general_group.add_argument(
        "--log",
        type=os.path.abspath,
        help="Path to the log file. Used to keep track of messages and errors printed on the stdout and stderr"
    )
    general_group.add_argument(
        "--max-genomes",
        type=number(int, minv=1),
        default=numpy.Inf,
        dest="max_genomes",
        help=(
            "Consider species with this number of genomes at most. "
            "There is not a maximum number of genomes per species by default"
        ),
    )
    general_group.add_argument(
        "--min-genomes",
        type=number(int, minv=1),
        default=1,
        dest="min_genomes",
        help=(
            "Consider species with a minimum number of genomes. "
            "There is not a minimum number of genomes per species by default"
        ),
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
        help="Maximum number of processors to process each NCBI tax ID in parallel",
    )
    general_group.add_argument(
        "--superkingdom",
        type=str,
        choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
        help="Consider genomes whose lineage belongs to a specific superkingdom if --input-list is not provided",
    )
    general_group.add_argument(
        "--tmp-dir",
        type=os.path.abspath,
        required=True,
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
        type=number(float, minv=0.0, maxv=1.0),
        default=1.0,
        help=(
            "Dereplicate genomes if they have a theta distance greather than this threshold. "
            "This is used exclusively in conjunction with the --dereplicate argument"
        ),
    )

    # Group of arguments for estimating the bloom filter size
    filter_size_group = p.add_argument_group("Bloom filter size")

    filter_size_group.add_argument(
        "--estimate-filter-size",
        action="store_true",
        default=False,
        dest="estimate_filter_size",
        help="Automatically estimate the best bloom filter size with ntCard",
    )
    filter_size_group.add_argument(
        "--filter-size",
        type=number(int, minv=10000),
        dest="filter_size",
        help="This is the size of the bloom filters",
    )
    filter_size_group.add_argument(
        "--increase-filter-size",
        type=number(float, minv=0.0, maxv=100.0),
        default=0.0,
        dest="increase_filter_size",
        help=(
            "Increase the estimated filter size by the specified percentage. "
            "This is used in conjunction with the --estimate-filter-size argument only. "
            "It is highly recommended to increase the filter size by a good percentage in case you are planning to update the index with new genomes"
        ),
    )

    # Group of arguments for estimating the optimal kmer size
    kitsune_group = p.add_argument_group("Kitsune: Estimation of the optimal kmer size")

    kitsune_group.add_argument(
        "--closely-related",
        action="store_true",
        default=False,
        dest="closely_related",
        help="For closesly related genomes use this flag",
    )
    kitsune_group.add_argument(
        "--estimate-kmer-size",
        action="store_true",
        default=False,
        dest="estimate_kmer_size",
        help="Automatically estimate the optimal kmer size with kitsune",
    )
    kitsune_group.add_argument(
        "--jellyfish-threads",
        type=number(int, minv=1, maxv=os.cpu_count()),
        default=1,
        dest="jellyfish_threads",
        help="Maximum number of threads for Jellyfish. This is required to maximise the kitsune performances",
    )
    kitsune_group.add_argument(
        "--limit-kmer-size",
        type=number(int, minv=4),
        default=32,
        dest="limit_kmer_size",
        help="Limit the estimation of the optimal kmer size with kitsune to this value at most",
    )

    return p.parse_args()


def download_assembly_summary(assembly_summary_url: str, folder_path: str) -> Dict[str, List[str]]:
    """
    Download and read the NCBI GenBank Assembly Summary

    :param assembly_summary_url:    URL to the NCBI GenBank Assembly Summary table
    :param folder_path:             Path to the folder in which the assembly table will be downloaded
    :return:                        Dictionary with the list of URLs to the genome files indexed by the species taxid
    """

    assembly_summary_filepath = os.path.join(folder_path, os.path.basename(assembly_summary_url))
    
    # Check whether the assembly table exists and it is GZ compressed
    if os.path.isfile("{}.gz".format(assembly_summary_filepath)):
        with open(assembly_summary_filepath, "w+") as file:
            run(["gzip", "-dc", "{}.gz".format(assembly_summary_filepath)], stdout=file, stderr=file)

    # Check whether the assembly summary table already exists
    # Otherwise, download it
    if not os.path.isfile(assembly_summary_filepath):
        assembly_summary_filepath = download(url=assembly_summary_url, folder=folder_path)

    # Raise an exception in case the assembly_summary_genbank.txt file does not exist
    if not os.path.isfile(assembly_summary_filepath):
        raise Exception("Unable to retrieve data from remote location\n{}".format(assembly_summary))
    
    assembly_summary = dict()

    with open(assembly_summary_filepath) as asf:
        # Skip the first line, it is just a comment
        next(asf)

        # Load the header line
        header = next(asf)[1:].strip().split("\t")

        for line in asf:
            if line.strip():
                line_split = line.split("\t")

                species_taxid = line_split[header.index("species_taxid")]
                ftp_path = line_split[header.index("ftp_path")]
                genome_url = os.path.join(ftp_path, "{}_genomic.fna.gz".format(os.path.basename(ftp_path)))

                # Check whether the genome URL is valid
                # Otherwise, skip the genome
                if species_taxid and validate_url(genome_url):
                    if species_taxid not in assembly_summary:
                        assembly_summary[species_taxid] = list()

                    species_info = dict()
                    
                    for h in header:
                        species_info[h] = line_split[header.index(h)].strip()

                    species_info["ftp_filepath"] = genome_url

                    _, local_filename, _, _ = get_file_info(genome_url, check_supported=False, check_exists=False)
                    species_info["local_filename"] = os.path.basename(local_filename)

                    assembly_summary[species_taxid].append(species_info)

    if not os.path.isfile("{}.gz".format(assembly_summary_filepath)):
        with open("{}.gz".format(assembly_summary_filepath), "w+") as file:
            run(["gzip", "-c", assembly_summary_filepath], stdout=file, stderr=file)
        
    if os.path.isfile(assembly_summary_filepath):
        os.unlink(assembly_summary_filepath)

    return assembly_summary


def download_taxdump(taxdump_url: str, folder_path: str) -> Tuple[str, str]:
    """
    Download and extract the NCBI taxdump tarball

    :param taxdump_url:     URL to the NCBI taxdump
    :param folder_path:     Path to the folder in which the taxdump tarball will be unpacked
    :return:                The nodes.dmp and names.dmp file paths
    """

    taxdump = download(url=taxdump_url, folder=folder_path)

    # Raise an exception in case the taxdump.tar.gz file does not exist
    if not os.path.isfile(taxdump):
        raise Exception("Unable to retrieve data from remote location\n{}".format(taxdump_url))

    # Create the taxdump folder in the temporary directory
    taxdump_dir = os.path.join(folder_path, "taxdump")
    os.makedirs(taxdump_dir, exist_ok=True)

    # Decompress the archive
    with tarfile.open(taxdump, "r:gz") as tar:
        tar.extractall(taxdump_dir)

    return os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp")


def level_name(current_level: str, prev_level: str) -> str:
    """
    Define a taxonomic level name

    :param current_level:   Current level name
    :param prev_level:      Previous level name in case of unclassified
    :return:                The new level name
    """

    # Remove special characters from current and previous level names
    current_level = re.sub(r"_+", "_", re.sub(r"\W+", "_", current_level)).strip("_")
    prev_level = re.sub(r"_+", "_", re.sub(r"\W+", "_", prev_level)).strip("_")

    # Build the new level name
    level_prefix = current_level.strip()
    level_suffix = ""
    if not level_prefix:
        level_prefix = prev_level
        # Fill empty taxa levels with unclassified
        level_suffix = "_unclassified"

    return "{}{}".format(level_prefix, level_suffix)


def load_taxa(
    ncbitax2lin_table: str,
    superkingdom: Optional[str] = None,
    kingdom: Optional[str] = None,
    dump: Optional[str] = None
) -> Tuple[list, list]:
    """
    Load the ncbitax2lin output table

    :param ncbitax2lin_table:   Path to the output table produced by ncbitax2lin
    :param superkingdom:        Consider a specific superkingdom only
    :param kingdom:             Kingdom
    :param dump:                Path to the output table
    :return:                    The lists of NCBI tax IDs and full taxonomic labels
    """

    # Take track of NCBI tax IDs and full taxonomic labels
    taxa_map: Dict[str, List[int]] = dict()

    with gzip.open(ncbitax2lin_table, "rt") as ncbi_table:
        # Load the first line as header and search for "superkingdom" and "kingdom" columns
        header = ncbi_table.readline().split(",")
        superkingdom_pos = header.index("superkingdom")  # Archaea, Bacteria, Eukaryota, Viruses
        kingdom_pos = header.index("kingdom")

        for line in ncbi_table:
            line = line.strip()
            if line:
                line_split = line.split(",")

                # Check whether the current taxonomy must be processed
                skip = True
                
                if not superkingdom and line_split[superkingdom_pos].strip():
                    skip = False
                
                elif line_split[superkingdom_pos] == superkingdom:
                    if not kingdom:
                        skip = False
                    
                    elif line_split[kingdom_pos] == kingdom:
                        skip = False

                if not skip:
                    # Build the current full taxonomic label
                    label = "k__{}|p__{}|c__{}|o__{}|f__{}|g__{}|s__{}".format(
                        superkingdom,  # Superkingdom
                        level_name(line_split[2], superkingdom if superkingdom else line_split[1]),  # Phylum
                        level_name(line_split[3], line_split[2]),  # Class
                        level_name(line_split[4], line_split[3]),  # Order
                        level_name(line_split[5], line_split[4]),  # Family
                        level_name(line_split[6], line_split[5]),  # Genus
                        level_name(line_split[7], line_split[6]),  # Species
                    )

                    # Exclude unclassified taxonomic labels
                    if "unclassified" not in label:
                        if label not in taxa_map:
                            taxa_map[label] = list()
                        taxa_map[label].append(int(line_split[0]))

    # Use a single tax ID for each label
    tax_labels = list()
    tax_ids = list()
    for tax_label in taxa_map:
        tax_labels.append(tax_label)
        tax_ids.append(str(min(taxa_map[tax_label])))

    if dump:
        # Dump the tax map to file
        with open(dump, "w+") as tax_table:
            # Add the header line
            tax_table.write("# {}\t{}\n".format("tax_label", "tax_id"))
            for pos, label in enumerate(tax_labels):
                tax_table.write("{}\t{}\n".format(label, tax_ids[pos]))

    return tax_ids, tax_labels


def ncbitax2lin(nodes: str, names: str, out_dir: str) -> str:
    """
    Extract lineages from the NCBI taxdump

    :param nodes:       Path to the nodes.dmp
    :param names:       Path to the names.dmp
    :param out_dir:     Path to the output folder
    :return:            Output table file path
    """

    ncbitax2lin_table = os.path.join(out_dir, "ncbi_lineages.csv.gz")
    run(
        [
            "ncbitax2lin",
            "--nodes-file",
            nodes,
            "--names-file",
            names,
            "--output",
            ncbitax2lin_table,
        ],
        silence=True,
    )

    # Raise an exception in case the ncbi_lineages.csv.gz file does not exist
    if not os.path.isfile(ncbitax2lin_table):
        raise Exception("An error has occurred while running ncbitax2lin")

    return ncbitax2lin_table


def quality_control(
    genomes: list,
    tax_id: str,
    tmp_dir: str,
    completeness: float = 0.0,
    contamination: float = 100.0,
    nproc: int = 1,
    pplacer_threads: int = 1,
) -> Tuple[List[str], List[str]]:
    """
    Quality-control genomes with CheckM

    :param genomes:             List of genome file paths
    :param tax_id:              NCBI tax ID
    :param tmp_dir:             Path to the temporary folder
    :param completeness:        Completeness threshold
    :param contamination:       Contamination threshold
    :param nproc:               Make CheckM parallel
    :param pplacer_threads:     Max number of threads for pplacer
    :return:                    List of genome file paths for genomes that passed the quality-control
                                in addition to a list with the CheckM output table paths
    """

    # Define the CheckM temporary folder
    checkm_tmp_dir = os.path.join(tmp_dir, "checkm", tax_id)
    os.makedirs(checkm_tmp_dir, exist_ok=True)

    # Run CheckM on the current set of genomes
    # Genomes are always downloaded as fna.gz
    checkm_tables = checkm(
        genomes,
        checkm_tmp_dir,
        file_extension="fna.gz",
        nproc=nproc,
        pplacer_threads=pplacer_threads,
    )

    # Filter genomes according to the input --completeness and --contamination thresholds
    genome_ids = filter_checkm_tables(checkm_tables, completeness=completeness, contamination=contamination)

    # Rebuild the genome file paths
    genome_paths = [os.path.join(tmp_dir, "genomes", tax_id, "{}.fna.gz".format(genome_id)) for genome_id in genome_ids]

    return genome_paths, checkm_tables


def organize_data(
    genomes: list,
    db_dir: str,
    tax_label: str,
    tax_id: str,
    cluster_id: int = 1,
    cluster_prefix: str = "MSBT",
    metadata: Optional[List[Dict[str, str]]] = None,
    checkm_tables: Optional[list] = None,
    flat_structure: bool = False,
) -> List[str]:
    """
    Organize genome files

    :param genomes:         List with path to the genome files
    :param db_dir:          Path to the database folder
    :param tax_label:       Full taxonomic label
    :param tax_id:          NCBI tax ID
    :param cluster_id:      Numberical cluster ID
    :param cluster_prefix:  Cluster prefix
    :param metadata:        List of dictionaries with genomes information
    :param checkm_tables:   List with paths to the CheckM output tables
    :param flat_structure:  Organize genomes in the same folder without any taxonomic organization
    :return:                List with genome file paths
    """

    genomes_paths = list()

    # In case at least one genome survived both the quality control and dereplication steps
    # Define the taxonomy folder in database
    tax_dir = os.path.join(db_dir, tax_label.replace("|", os.sep)) if not flat_structure else db_dir
    genomes_dir = os.path.join(tax_dir, "genomes")
    os.makedirs(genomes_dir, exist_ok=True)

    references_path = "references.txt" if not flat_structure else "genomes_{}.txt".format(tax_id)
    with open(os.path.join(tax_dir, references_path), "w+") as refsfile:
        # Move the processed genomes to the taxonomy folder
        for genome_path in genomes:
            # Get the genome name from file path
            _, genome_name, _, _ = get_file_info(genome_path)

            # Move the genome file into the species folder
            shutil.move(genome_path, genomes_dir)

            # Also take track of the genome names in the references.txt file
            refsfile.write("{}\n".format(genome_name))

            # Add the genome to the full list of genomes in database
            genomes_paths.append(os.path.join(genomes_dir, os.path.basename(genome_path)))

    # Create the metadata table in the taxonomy folder
    with open(os.path.join(tax_dir, "metadata.tsv"), "w+") as metafile:
        if not flat_structure:
            metafile.write("# Cluster ID: {}{}\n".format(cluster_prefix, cluster_id))

        if metadata:
            header_list = list(metadata[0].keys())
            metafile.write("# {}\n".format("\t".join(header_list)))

            for genome_info in metadata:
                if os.path.join(genomes_dir, "{}.fna.gz".format(genome_info["local_filename"])) in genomes_paths:
                    line = list()

                    for h in header_list:
                        line.append(genome_info[h])

                    metafile.write("{}\n".format("\t".join(line)))

    if checkm_tables:
        # Also merge the CheckM output tables and move the result to the taxonomy folder
        checkm_path = "checkm.tsv" if not flat_structure else "checkm_{}.tsv".format(tax_id)
        with open(os.path.join(tax_dir, checkm_path), "w+") as table:
            header = True
            for table_path in checkm_tables:
                with open(table_path) as partial_table:
                    line_count = 0
                    for line in partial_table:
                        line = line.strip()
                        if line:
                            if line_count == 0:
                                if header:
                                    table.write("{}\n".format(line))
                                    header = False
                            else:
                                table.write("{}\n".format(line))
                            line_count += 1

    return genomes_paths


def process_input_genomes(
    genomes_list: list,
    taxonomic_label: str,
    cluster_id: int,
    db_dir: str,
    tmp_dir: str,
    kmer_len: int,
    cluster_prefix: str = "MSBT",
    nproc: int = 1,
    pplacer_threads: int = 1,
    completeness: float = 0.0,
    contamination: float = 100.0,
    dereplicate: bool = False,
    similarity: float = 1.0,
    flat_structure: bool = False,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> List[str]:
    """
    Process a provided list of input genomes
    Organize, quality-control, and dereplicate genomes

    :param genomes_list:        List of input genome file paths
    :param taxonomic_label:     Taxonomic label of the input genomes
    :param cluster_id:          Numberical cluster ID
    :param db_dir:              Path to the database root folder
    :param tmp_dir:             Path to the temporary folder
    :param kmer_len:            Length of the kmers
    :param cluster_prefix:      Cluster prefix
    :param nproc:               Make the process parallel when possible
    :param pplacer_threads:     Maximum number of threads to make pplacer parallel with CheckM
    :param completeness:        Threshold on the CheckM completeness
    :param contamination:       Threshold on the CheckM contamination
    :param dereplicate:         Enable the dereplication step to get rid of replicated genomes
    :param similarity:          Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param flat_structure:      Do not taxonomically organize genomes
    :param logger:              Logger object
    :param verbose:             Print messages on screen
    :return:                    The list of paths to the genome files
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    if verbose:
        printline("Processing {}".format(taxonomic_label))

    # Create a temporary folder to store genomes
    tax_id = hashlib.md5(taxonomic_label.encode("utf-8")).hexdigest()
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes", tax_id)
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # Iterate over the list of input genome file paths
    genomes = list()
    for genome_path in genomes_list:
        # Symlink input genomes into the temporary folder
        os.symlink(genome_path, os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))
        genomes.append(os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))

    # Take track of the paths to the CheckM output tables
    checkm_tables: List[str] = list()

    # Quality control
    if completeness > 0.0 or contamination < 100.0:
        before_qc = len(genomes)

        genomes, checkm_tables = quality_control(
            genomes,
            tax_id,
            tmp_dir,
            completeness=completeness,
            contamination=contamination,
            nproc=nproc,
            pplacer_threads=pplacer_threads,
        )

        if before_qc > len(genomes) and verbose:
            printline("Quality control: excluding {}/{} genomes".format(before_qc - len(genomes), before_qc))

    # Dereplication
    if len(genomes) > 1 and dereplicate:
        before_dereplication = len(genomes)

        genomes = dereplicate_genomes(
            genomes,
            tax_id,
            tmp_dir,
            kmer_len,
            filter_size=None,
            nproc=nproc,
            similarity=similarity,
        )

        if before_dereplication > len(genomes) and verbose:
            printline(
                "Dereplication: excluding {}/{} genomes".format(
                    before_dereplication - len(genomes), before_dereplication
                )
            )

    # Check whether no genomes survived the quality control and the dereplication steps
    if not genomes:
        if verbose:
            printline("No more genomes available for the NCBI tax ID {}".format(tax_id))
        return genomes

    # Organize genome files
    genomes_paths = organize_data(
        genomes,
        db_dir,
        taxonomic_label,
        tax_id,
        cluster_id=cluster_id,
        cluster_prefix=cluster_prefix,
        metadata=None,
        checkm_tables=checkm_tables,
        flat_structure=flat_structure,
    )

    return genomes_paths


def process_tax_id(
    tax_id: str,
    tax_label: str,
    genomes_info: List[Dict[str, str]],
    cluster_id: int,
    db_dir: str,
    tmp_dir: str,
    kmer_len: int,
    cluster_prefix: str = "MSBT",
    limit_genomes: float = numpy.Inf,
    max_genomes: float = numpy.Inf,
    min_genomes: float = 1.0,
    nproc: int = 1,
    pplacer_threads: int = 1,
    completeness: float = 0.0,
    contamination: float = 100.0,
    dereplicate: bool = False,
    similarity: float = 1.0,
    flat_structure: bool = False,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> List[str]:
    """
    Process a specific NCBI tax ID
    Download reference genomes only, perform quality-control and dereplication

    :param tax_id:              NCBI tax ID of a species
    :param tax_label:           Full taxonomic label
    :param genomes_info:        List of dictionaries with genomes information
    :param cluster_id:          Numberical cluster ID
    :param db_dir:              Path to the database root folder
    :param tmp_dir:             Path to the temporary folder
    :param kmer_len:            Length of the kmers
    :param cluster_prefix:      Cluster prefix
    :param limit_genomes:       Limit the number of genomes per species
    :param max_genomes:         Consider species with this number of genomes at most
    :param min_genomes:         Consider species with a minimum number of genomes
    :param nproc:               Make the process parallel when possible
    :param pplacer_threads:     Maximum number of threads to make pplacer parallel with CheckM
    :param completeness:        Threshold on the CheckM completeness
    :param contamination:       Threshold on the CheckM contamination
    :param dereplicate:         Enable the dereplication step to get rid of replicated genomes
    :param similarity:          Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param flat_structure:      Do not taxonomically organize genomes
    :param logger:              Logger object
    :param verbose:             Print messages on screen
    :return:                    The list of paths to the genome files
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Create a temporary folder to store the downloaded genomes
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes", tax_id)
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # Take track of the paths to the genome files
    genomes_paths = list()

    # Take track of genomes URLs
    genomes_urls = list()

    for genome_info in genomes_info:
        if not genome_info["excluded_from_refseq"].strip() or all([ex.strip() in REFERENCE_TAGS for ex in genome_info["excluded_from_refseq"].split(";")]):
            genomes_urls.append(genome_info["ftp_filepath"])

    if len(genomes_urls) > 0:
        if verbose:
            printline("Processing NCBI tax ID {} with {} genomes".format(tax_id, len(genomes_urls)))

        # Check whether the number of genomes falls in the --min-genomes --max-genomes interval
        if len(genomes_urls) < min_genomes or len(genomes_urls) > max_genomes:
            if verbose:
                printline(
                    "WARNING: the number of genomes does not respect the minimum and maximum number of genomes allowed"
                )
            return list()

        # Check whether the number of genomes is greater than the number provided with --limit-genomes
        if len(genomes_urls) > limit_genomes:
            if verbose:
                printline("WARNING: the number of genomes per species is limited to {} at most".format(limit_genomes))

            # Always use the same seed for reproducibility
            rng = numpy.random.default_rng(0)

            # Subsampling genomes as input for kitsune
            rng.shuffle(genomes_urls)
            genomes_urls = genomes_urls[:limit_genomes]

        # Keep track of the paths to the downloaded genome files
        genomes = list()

        for genome_url in genomes_urls:
            genome_filepath = download(url=genome_url, folder=tmp_genomes_dir, raise_exception=False)

            if genome_filepath and integrity_check(genome_filepath):
                genomes.append(genome_filepath)

        if verbose:
            printline("Retrieved {}/{} genomes".format(len(genomes), len(genomes_urls)))

        # Take track of the paths to the CheckM output tables
        checkm_tables: List[str] = list()

        # Quality control
        if genomes and (completeness > 0.0 or contamination < 100.0):
            before_qc = len(genomes)

            genomes, checkm_tables = quality_control(
                genomes,
                tax_id,
                tmp_dir,
                completeness=completeness,
                contamination=contamination,
                nproc=nproc,
                pplacer_threads=pplacer_threads,
            )

            if before_qc > len(genomes) and verbose:
                printline("Quality control: excluding {}/{} genomes".format(before_qc - len(genomes), before_qc))

        # Dereplication
        if len(genomes) > 1 and dereplicate:
            before_dereplication = len(genomes)

            genomes = dereplicate_genomes(
                genomes,
                tax_id,
                tmp_dir,
                kmer_len,
                filter_size=None,
                nproc=nproc,
                similarity=similarity,
            )

            if before_dereplication > len(genomes) and verbose:
                printline(
                    "Dereplication: excluding {}/{} genomes".format(
                        before_dereplication - len(genomes), before_dereplication
                    )
                )

        # Check whether no genomes survived the quality control and the dereplication steps
        if not genomes:
            if verbose:
                printline("No more genomes available for the NCBI tax ID {}".format(tax_id))

            return list()

        # Organize genome files
        genomes_paths = organize_data(
            genomes,
            db_dir,
            tax_label,
            tax_id,
            cluster_id=cluster_id,
            cluster_prefix=cluster_prefix,
            metadata=genomes_info,
            checkm_tables=checkm_tables,
            flat_structure=flat_structure,
        )

    # Return the list of paths to the genome files
    return genomes_paths


def get_sublist(genomes, limit_number=None, limit_percentage=100.0, flat_structure=False) -> List[str]:
    """
    Given a list of genomes, define a sublist with a limited number of genomes taken randomly

    :param genomes:             Input list of paths to the genomes files
    :param limit_number:        Limit the number of elements in the resulting list to this number
                                If defined, it overrides the percentage
    :param limit_percentage:    Limit the size of the resulting list to this percentage computed 
                                on the total number of elements in the input list
    :param flat_structure:      Genomes organization
    :return:                    The sublist of genomes
    """

    if not genomes or limit_number == 0 or limit_percentage == 0.0:
        return list()

    # Always use the same seed for reproducibility
    rng = numpy.random.default_rng(0)

    # Number of genomes to be considered for estimating the optimal kmer size with kitsune
    select_at_most = 0
    is_percentage = False

    if isinstance(limit_number, int):
        # Select a specific number of genomes at most
        select_at_most = limit_number

    else:
        # Remember that the number of genomes will be selected based on a percentage
        is_percentage = True

    if is_percentage and limit_percentage == 100.0:
        # There is no need for subsampling here
        return genomes

    genomes_sub = list()

    if flat_structure:
        if select_at_most == 0:
            # Select genomes as a percentage of the total number of genomes
            select_at_most = int(math.ceil(len(genomes) * limit_percentage / 100.0))

        # Subsampling genomes as input for kitsune
        rng.shuffle(genomes)
        genomes_sub = genomes[:select_at_most]

    else:
        # Genomes must be taxonomically reorganized
        species2genomes = dict()

        for genome_path in genomes:
            # Get the genomes folder
            dirpath, _, _, _ = get_file_info(genome_path)

            # Retrieve the species id
            # Genomes are located under the "genomes" folder, under the species directory
            species = dirpath.split(os.sep)[-2]

            # Group genomes according to their species
            if species not in species2genomes:
                species2genomes[species] = list()
            species2genomes[species].append(genome_path)

        for species in species2genomes:
            species_genomes = species2genomes[species]

            # Cluster specific maximum number of genomes to be considered for kitsune
            select_at_most_in_species = select_at_most

            if select_at_most_in_species == 0:
                # Select genomes as a percentage of the total number of genomes
                select_at_most_in_species = int(math.ceil(len(species_genomes) * limit_percentage / 100.0))

            # Subsampling genomes as input for kitsune
            rng.shuffle(species_genomes)
            genomes_sub.extend(species_genomes[:select_at_most_in_species])

    return genomes_sub


def index(
    db_dir: str,
    input_list: str,
    superkingdom: str,
    tmp_dir: str,
    kingdom: Optional[str] = None,
    cluster_prefix: str = "MSBT",
    kmer_len: Optional[int] = None,
    filter_size: Optional[int] = None,
    flat_structure: bool = False,
    limit_genomes: float = numpy.Inf,
    max_genomes: float = numpy.Inf,
    min_genomes: float = 1.0,
    estimate_filter_size: bool = True,
    increase_filter_size: float = 0.0,
    estimate_kmer_size: bool = True,
    limit_estimation_number: Optional[int] = None,
    limit_estimation_percentage: float = 100.0,
    limit_kmer_size: int = 32,
    completeness: float = 0.0,
    contamination: float = 100.0,
    dereplicate: bool = False,
    similarity: float = 1.0,
    closely_related: bool = False,
    logger: Optional[Logger] = None,
    verbose: bool = False,
    nproc: int = 1,
    pplacer_threads: int = 1,
    jellyfish_threads: int = 1,
    parallel: int = 1,
) -> None:
    """
    Build the database baseline

    :param db_dir:                      Path to the database root folder
    :param input_list:                  Path to the file with a list of input genome paths
    :param superkingdom:                Retrieve genomes that belong to a specific superkingdom
    :param tmp_dir:                     Path to the temporary folder
    :param kingdom:                     Kingdom
    :param cluster_prefix:              Prefix of clusters numerical identifiers
    :param kmer_len:                    Length of the kmers
    :param filter_size:                 Size of the bloom filters
    :param flat_structure:              Do not taxonomically organize genomes
    :param limit_genomes:               Limit the number of genomes per species
    :param max_genomes:                 Consider species with this number of genomes at most
    :param min_genomes:                 Consider species with a minimum number of genomes
    :param estimate_filter_size:        Run ntCard to estimate the most appropriate bloom filter size
    :param increase_filter_size:        Increase the estimated bloom filter size by the specified percentage
    :param estimate_kmer_size:          Run kitsune to estimate the best kmer size
    :param limit_estimation_number:     Number of genomes per group as input for kitsune and ntCard
    :param limit_estimation_percentage: Percentage on the total number of genomes per group as input for kitsune and ntCard
    :param limit_kmer_size:             Maximum kmer size for kitsune kopt
    :param completeness:                Threshold on the CheckM completeness
    :param contamination:               Threshold on the CheckM contamination
    :param dereplicate:                 Enable the dereplication step to get rid of replicated genomes
    :param closely_related:             For closesly related genomes use this flag
    :param similarity:                  Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param logger:                      Logger object
    :param verbose:                     Print messages on screen
    :param nproc:                       Make the process parallel when possible
    :param pplacer_threads:             Maximum number of threads to make pplacer parallel with CheckM
    :param jellyfish_threads:           Maximum number of threads to make Jellyfish parallel with kitsune
    :param parallel:                    Maximum number of processors to process each NCBI tax ID in parallel
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    if input_list:
        # Load the set of input genomes
        taxonomy2genomes: Dict[str, List[str]] = dict()
        with open(input_list) as file:
            for line in file:
                line = line.strip()
                if line:
                    if not line.startswith("#"):
                        line_split = line.split("\t")

                        taxonomy = "NA"
                        if len(line_split) == 2:
                            taxonomy = line_split[1]

                            if len(taxonomy.split("|")) != 7:
                                # Taxonomic labels must have 7 levels
                                raise Exception(
                                    "Invalid taxonomic label! Please note that taxonomies must have 7 levels:\n{}".format(
                                        line_split[1]
                                    )
                                )

                        # This automatically check whether extension and compression are supported
                        dirpath, genome_name, extension, compression = get_file_info(line_split[0])

                        if taxonomy not in taxonomy2genomes:
                            taxonomy2genomes[taxonomy] = list()

                        taxonomy2genomes[taxonomy].append(
                            os.path.join(dirpath, "{}{}{}".format(
                                genome_name, extension, compression if compression else ""
                            ))
                        )

        # Force flat structure in case of genomes with no taxonomic label
        if "NA" in taxonomy2genomes:
            flat_structure = True

        # Build a partial function around process_input_genomes
        process_partial = partial(
            process_input_genomes,
            db_dir=db_dir,
            tmp_dir=tmp_dir,
            kmer_len=kmer_len,
            cluster_prefix=cluster_prefix,
            nproc=nproc,
            pplacer_threads=pplacer_threads,
            completeness=completeness,
            contamination=contamination,
            dereplicate=dereplicate,
            similarity=similarity,
            flat_structure=flat_structure,
            logger=logger,
            verbose=False,
        )

        printline("Processing clusters")

        # Define the length of the progress bar
        pbar_len = len(taxonomy2genomes)

    else:
        # Download the NCBI taxdump
        printline("Downloading taxonomy dump from NCBI")

        # Recover already processed nodes.dmp and names.dmp
        nodes_dmp = os.path.join(tmp_dir, "taxdump", "nodes.dmp")
        names_dmp = os.path.join(tmp_dir, "taxdump", "names.dmp")
        if not os.path.isfile(nodes_dmp) or not os.path.isfile(names_dmp):
            nodes_dmp, names_dmp = download_taxdump(TAXDUMP_URL, tmp_dir)

        # Run ncbitax2lin to extract lineages
        printline("Exporting lineages")
        ncbitax2lin_table = os.path.join(tmp_dir, "ncbi_lineages.csv.gz")
        if not os.path.isfile(ncbitax2lin_table):
            ncbitax2lin_table = ncbitax2lin(nodes_dmp, names_dmp, tmp_dir)

        # Build a mapping between tax IDs and full taxonomic labels
        # Take track of the taxonomic labels to remove duplicates
        # Taxonomy IDs are already sorted in ascending order in the ncbitax2lin output table
        printline("Building species tax ID to full taxonomy mapping")
        tax_ids, tax_labels = load_taxa(
            ncbitax2lin_table,
            superkingdom=superkingdom,
            kingdom=kingdom,
            dump=os.path.join(tmp_dir, "taxa.tsv") if not os.path.isfile(os.path.join(tmp_dir, "taxa.tsv")) else None,
        )

        # Download the NCBI GenBank Assembly Summary table
        # and load the list of genomes grouped by species taxid
        printline("Downloading NCBI GenBank Assembly Summary table")
        assembly_summary = download_assembly_summary(ASSEMBLY_SUMMARY_URL, tmp_dir)

        # Build a partial function around process_tax_id
        process_partial = partial(
            process_tax_id,
            db_dir=db_dir,
            tmp_dir=tmp_dir,
            kmer_len=kmer_len,
            cluster_prefix=cluster_prefix,
            limit_genomes=limit_genomes,
            max_genomes=max_genomes,
            min_genomes=min_genomes,
            nproc=nproc,
            pplacer_threads=pplacer_threads,
            completeness=completeness,
            contamination=contamination,
            dereplicate=dereplicate,
            similarity=similarity,
            flat_structure=flat_structure,
            logger=logger,
            verbose=False,
        )

        printline("Processing clusters")

        # Reshape tax_ids and tax_labels according to the assembly summary
        for txid in copy.deepcopy(tax_ids):
            if txid not in assembly_summary:
                tax_pos = tax_ids.index(txid)

                # Delete a tax ID
                del tax_ids[tax_pos]
                del tax_labels[tax_pos]
        
        # Define the length of the progress bar
        pbar_len = len(tax_ids)

    # Define a cluster counter
    clusters_counter = pbar_len

    # Take track of all the genomes paths
    genomes_paths = list()

    with mp.Pool(processes=parallel) as pool, tqdm.tqdm(total=pbar_len, disable=(not verbose)) as pbar:
        # Wrapper around the update function of tqdm
        def progress(*args):
            pbar.update()

        if input_list:
            # Process input genomes
            jobs = [
                pool.apply_async(
                    process_partial,
                    args=(taxonomy2genomes[taxonomy], taxonomy, pos + 1),
                    callback=progress,
                )
                for pos, taxonomy in enumerate(taxonomy2genomes)
            ]

        else:
            # Process the NCBI tax IDs
            jobs = [
                pool.apply_async(
                    process_partial,
                    args=(tax_id, tax_labels[pos], assembly_summary[tax_id], pos + 1),
                    callback=progress,
                )
                for pos, tax_id in enumerate(tax_ids) if tax_id in assembly_summary
            ]

        # Get results from jobs
        for job in jobs:
            genomes_paths.extend(job.get())

    if not genomes_paths:
        raise Exception("No input genomes found")

    printline("Processing {} genomes".format(len(genomes_paths)))

    # Define the manifest file path
    # This is used in case the --filter-size and/or --kmer-len must be estimated
    manifest_filepath = os.path.join(db_dir, "manifest.txt")

    # Add cluster counter
    with open(manifest_filepath, "a+") as manifest:
        manifest.write("--clusters-counter {}\n".format(clusters_counter))

    # Limited set of genomes in case of --estimate-kmer-size and/or --estimate-filter-size
    # The set is limited to the number of genomes specified with --limit-estimation-number or limit_estimation_percentage
    genomes_paths_sub = list()

    # Check whether the kmer size must be estimated
    if estimate_kmer_size and not kmer_len:
        printline("Estimating the best k-mer size. Please be patient, this may take a while")

        # Define a subset of genomes
        genomes_paths_sub = get_sublist(
            genomes_paths,
            limit_number=limit_estimation_number,
            limit_percentage=limit_estimation_percentage,
            flat_structure=flat_structure
        )

        if len(genomes_paths_sub) < 2:
            raise Exception("Not enough genomes for estimating the optimal kmer size with kitsune!")

        # Estimate a kmer size
        kmer_len = optimal_k(
            genomes_paths_sub,
            limit_kmer_size,
            os.path.join(tmp_dir, "kitsune"),
            closely_related=closely_related,
            nproc=nproc,
            threads=jellyfish_threads
        )

        # Update the manifest file with the --filter-size
        with open(manifest_filepath, "a+") as manifest:
            manifest.write("--kmer-len {}\n".format(kmer_len))

    # Check whether the bloom filter size must be estimated
    if estimate_filter_size and not filter_size:
        printline("Estimating the bloom filter size")

        # Use the precomputed subset of genomes in case it has already been
        # defined because of the --estimate-kmer-size
        if not genomes_paths_sub:
            # Define a subset of genomes
            genomes_paths_sub = get_sublist(
                genomes_paths,
                limit_number=limit_estimation_number,
                limit_percentage=limit_estimation_percentage,
                flat_structure=flat_structure
            )

        # Estimate the bloom filter size
        filter_size = estimate_bf_size(
            genomes_paths_sub,
            kmer_len,
            os.path.join(tmp_dir, "genomes"),
            tmp_dir,
            nproc=nproc,
        )

        # Increment the estimated bloom filter size
        increment = int(math.ceil(filter_size * increase_filter_size / 100.0))
        filter_size += increment

        # Update the manifest file with the --filter-size
        with open(manifest_filepath, "a+") as manifest:
            manifest.write("--filter-size {}\n".format(filter_size))

    if filter_size and kmer_len:
        # Retrieve the current working directory
        current_folder = os.getcwd()

        printline("Indexing genomes")

        # Define a partial function around the howdesbt wrapper
        howdesbt_partial = partial(
            howdesbt,
            kmer_len=kmer_len,
            filter_size=filter_size,
            nproc=nproc,
            flat_structure=flat_structure,
        )

        if not flat_structure:
            # Iterate over all the taxonomic levels from the species up to the superkingdom
            for level in [
                "species",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "kingdom",
            ]:
                folders = [str(path) for path in Path(db_dir).glob("**/{}__*".format(level[0]))]
                printline("Running HowDeSBT at the {} level ({} clusters)".format(level, len(folders)))

                with mp.Pool(processes=parallel) as pool, tqdm.tqdm(total=len(folders), disable=(not verbose)) as pbar:
                    # Wrapper around the update function of tqdm
                    def progress(*args):
                        pbar.update()

                    # Process the NCBI tax IDs
                    jobs = [
                        pool.apply_async(howdesbt_partial, args=(level_dir,), callback=progress)
                        for level_dir in folders
                    ]

                    for job in jobs:
                        job.get()

        # Also run HowDeSBT on the database folder to build
        # the bloom filter representation of the superkingdom
        if not flat_structure:
            printline("Building the database root bloom filter with HowDeSBT")

        howdesbt_partial(db_dir)

        if flat_structure:
            # Merge all the genomes_*.txt into a single file
            gen = Path(db_dir).glob("genomes_*.txt")
            with open(os.path.join(db_dir, "genomes.txt"), "w+") as genomes_file:
                for filepath in gen:
                    with open(str(filepath)) as file:
                        for line in file:
                            line = line.strip()
                            genomes_file.write("{}\n".format(line))

                    # Get rid of the partial genomes file
                    os.unlink(str(filepath))

        # The howdesbt function automatically set the current working directory to
        # the index folder of the taxonomic labels
        # Come back to the original folder
        os.chdir(current_folder)


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    if not args.input_list and not args.superkingdom:
        raise Exception("Please specify --input-list or a valid --superkingdom")

    if args.superkingdom:
        # Check whether the database folder exists
        if os.path.isdir(os.path.join(args.db_dir, "k__{}".format(args.superkingdom))):
            raise Exception(
                (
                    "An indexed version of the {} superkingdom already exists in the database!\n"
                    "Please use the update module to add new genomes"
                ).format(args.superkingdom)
            )

    # Create the database folder
    os.makedirs(args.db_dir, exist_ok=True)

    # Skip the bloom filter size estimation if the filter size is passed as input
    if not args.filter_size and not args.estimate_filter_size:
        raise Exception(
            (
                "Please specify a bloom filter size with the --filter-size option or "
                "use the --estimate-filter-size flag to automatically estimate the best bloom filter size with ntCard"
            )
        )

    # Skip the kmer size estimation if the kmer length is passed as input
    if not args.kmer_len and not args.estimate_kmer_size:
        raise Exception(
            (
                "Please specify a kmer size with the --kmer-len option or "
                "use the --estimate-kmer-size flag to automatically estimate the optimal kmer size with kitsune"
            )
        )

    superkingdoms = set()
    if not args.flat_structure and not args.superkingdom:
        # In this case --input-list must be specified
        # Build as many manifest file as the number of superkingdoms in the input list of taxonomic labels
        try:
            with open(args.input_list) as inlist:
                for line in inlist:
                    line = line.strip()
                    if line:
                        superkingdoms.add(line.split("\t")[1].split("|")[0].split("__")[-1])

        except Exception as ex:
            raise Exception("Input file is not correctly formatted: {}".format(args.input_list)).with_traceback(
                ex.__traceback__
            )

    # We can process a single superkingdom per run
    if len(superkingdoms) > 1:
        raise Exception("The {} module cannot process more than 1 superkingdom per run!".format(TOOL_ID))

    if not args.superkingdom and len(superkingdoms):
        args.superkingdom = superkingdoms[0]

    if args.superkingdom:
        # Create the superkingdom directory
        os.makedirs(os.path.join(args.db_dir, "k__{}".format(args.superkingdom)), exist_ok=True)

    # Define the database manifest file
    manifest_filepath = os.path.join(args.db_dir, "manifest.txt")

    if os.path.isfile(manifest_filepath):
        # Load and compare --kmer-len and --filter-size
        # This is require to resume the process
        with open(manifest_filepath) as manifest:
            for line in manifest:
                line = line.strip()
                if line:
                    line_split = line.split(" ")

                    if line_split[0] == "--kmer-len":
                        if not args.kmer_len:
                            args.kmer_len = int(line_split[1])

                        elif args.kmer_len != int(line_split[1]):
                            raise Exception("The kmer length is not compatible with the specified database")

                    elif line_split[0] == "--filter-size":
                        if not args.filter_size:
                            args.filter_size = int(line_split[1])

                        elif args.filter_size != int(line_split[1]):
                            raise Exception("The bloom filter size is not compatible with the specified database")

    else:
        # Initialize manifest file
        with open(manifest_filepath, "w+") as manifest:
            if args.kmer_len:
                manifest.write("--kmer-len {}\n".format(args.kmer_len))
            
            if args.filter_size:
                manifest.write("--filter-size {}\n".format(args.filter_size))

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Build a sh script with the command line used to launch the index module
    build_sh(sys.argv, TOOL_ID, args.db_dir)

    t0 = time.time()

    index(
        args.db_dir,
        args.input_list,
        args.superkingdom,
        args.tmp_dir,
        kingdom=args.kingdom,
        cluster_prefix=args.cluster_prefix,
        kmer_len=args.kmer_len,
        filter_size=args.filter_size,
        flat_structure=args.flat_structure,
        limit_genomes=args.limit_genomes,
        max_genomes=args.max_genomes,
        min_genomes=args.min_genomes,
        estimate_filter_size=args.estimate_filter_size,
        increase_filter_size=args.increase_filter_size,
        estimate_kmer_size=args.estimate_kmer_size,
        limit_estimation_number=args.limit_estimation_number,
        limit_estimation_percentage=args.limit_estimation_percentage,
        limit_kmer_size=args.limit_kmer_size,
        completeness=args.completeness,
        contamination=args.contamination,
        dereplicate=args.dereplicate,
        similarity=args.similarity,
        closely_related=args.closely_related,
        logger=logger,
        verbose=args.verbose,
        nproc=args.nproc,
        pplacer_threads=args.pplacer_threads,
        jellyfish_threads=args.jellyfish_threads,
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
