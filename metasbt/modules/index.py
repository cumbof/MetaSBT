#!/usr/bin/env python3
"""Build a database with a set of genomes indexed with HowDeSBT.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.3"
__date__ = "Sep 24, 2024"

import argparse as ap
import errno
import itertools
import math
import multiprocessing as mp
import os
import shutil
import sys
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
        bfaction,
        build_sh,
        dereplicate_genomes,
        estimate_bf_size,
        filter_quality,
        get_file_info,
        howdesbt,
        init_logger,
        load_input_table,
        load_manifest,
        number,
        optimal_k,
        println,
        strand,
        quality,
    )
except Exception:
    pass

# Define the module name
TOOL_ID = "index"

# Define the list of dependencies
DEPENDENCIES = [
    "checkm2",
    "checkv",
    "eukcc",
    "howdesbt",
    "kitsune",
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
        description="Build a database with a set of genomes indexed with HowDeSBT",
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
        required = "--resume" not in argv,
        dest="input_list",
        help=(
            "Path to the input table with a list of genome file paths and an optional column with their taxonomic labels. "
            "Please note that the input genome files must all have the same extension and can be Gzip compressed (e.g.: *.fna.gz)"
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
        help="Maximum number of processors to process each cluster in parallel",
    )
    general_group.add_argument(
        "--parent-version",
        type=str,
        dest="parent_version",
        help="Version of MetaSBT that called this process"
    )
    general_group.add_argument(
        "--resume",
        type=os.path.abspath,
        help=(
            "Path to the \"index.sh\" file with the configuration of a previous run. "
            "Used to resume the index process in case of an unexpected error. "
            "Not available with --flat-structure"
        ),
    )
    general_group.add_argument(
        "--uniform-strand",
        action="store_true",
        default=False,
        dest="uniform_strand",
        help=(
            "Preprocess the input fasta sequences to make them all on the same strand. "
            "Must be used in case of viral sequences"
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
        "--use-representatives",
        action="store_true",
        default=False,
        dest="use_representatives",
        help=(
            "Use only 3 representative genomes per species choosen as the ones that minimize and maximize the "
            "distance with all the other genomes in the same cluster, plus the \"centroid\". "
            "These may change in case of an update. "
            "Not available with --flat-structure"
        )
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
    filter_size_group.add_argument(
        "--min-kmer-occurrences",
        type=number(int, minv=1),
        default=2,
        dest="min_kmer_occurrences",
        help=(
            "Minimum number of occurrences of kmers to be considered for estimating the bloom filter size "
            "and for building the bloom filter files"
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

    return p.parse_args(argv)


def quality_control(
    genomes: List[os.path.abspath],
    tax_id: str,
    tmp_dir: os.path.abspath,
    method: str,
    input_extension: str="fna.gz",
    completeness: float=0.0,
    contamination: float=100.0,
    nproc: int=1,
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """Quality-control genomes.

    Parameters
    ----------
    genomes : list
        List of genome file paths.
    tax_id : str
        NCBI tax ID.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    method : {"checkm2", "checkv", "eukcc"}
        Quality-control method.
    input_extension : str, default "fna.gz"
        File extension of the input genome files.
    completeness : float, default 0.0
        Completeness threshold.
    contamination : float, default 100.0
        Contamination threshold.
    nproc : int, default 1
        Make it paralell when possible.

    Returns
    -------
    tuple
        A tuple with the list of genome file paths for the genomes that passes the quality-control,
        in addition to a dictionary with quality information indexed by the genome names.
    """

    # Define the temporary folder
    quality_tmp_dir = os.path.join(tmp_dir, "quality", tax_id)
    os.makedirs(quality_tmp_dir, exist_ok=True)

    # Run quality on the current set of genomes
    quality_dict = quality(
        method=method,
        args={
            "genomes_paths": genomes,
            "tmp_dir": quality_tmp_dir,
            "file_extension": input_extension,
            "nproc": nproc
        }
    )

    # Filter genomes according to the input --completeness and --contamination thresholds
    genome_ids = filter_quality(quality_dict, completeness=completeness, contamination=contamination)

    # Rebuild the genome file paths
    genome_paths = [os.path.join(tmp_dir, "genomes", tax_id, "{}.{}".format(genome_id, input_extension)) for genome_id in genome_ids]

    return genome_paths, quality_dict


def organize_data(
    genomes: list,
    db_dir: os.path.abspath,
    tax_label: str,
    cluster_id: int=1,
    cluster_prefix: str="MSBT",
    metadata: Optional[List[Dict[str, str]]]=None,
    quality_dict: Optional[Dict[str, Dict[str, str]]]=None,
    flat_structure: bool=False,
) -> List[os.path.abspath]:
    """Organize genome files.

    Parameters
    ----------
    genomes : list
        List with the path to the genome files.
    db_dir : os.path.abspath
        Path to the database root folder.
    tax_label : str
        Full taxonomic label.
    cluster_id: int, default 1
        Numerical cluster ID.
    cluster_prefix : str, default "MSBT"
        Cluster prefix.
    metadata : list, optional
        List of dictionaries with genomes information.
    quality_dict : dict, optional
        Dictionary with the quality stats.
    flat_structure : bool, default False
        Organize genomes in the same folder without any taxonomic organization if True.

    Returns
    -------
    list
        A list with the genome file paths.
    """

    genomes_paths = list()

    # In case at least one genome survived both the quality control and dereplication steps
    # Define the taxonomy folder in database
    tax_dir = os.path.join(db_dir, tax_label.replace("|", os.sep)) if not flat_structure else db_dir
    genomes_dir = os.path.join(tax_dir, "genomes")
    os.makedirs(genomes_dir, exist_ok=True)

    if not flat_structure:
        # Create the strains folder
        genomes_dir = os.path.join(tax_dir, "strains", "genomes")
        os.makedirs(genomes_dir, exist_ok=True)

    references_path = "references.txt" if not flat_structure else "genomes.txt"
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

    if quality_dict:
        # Also merge the quality tables and move the result to the taxonomy folder
        with open(os.path.join(tax_dir, "quality.tsv"), "w+") as table:
            header = None

            for genome_id in quality_dict:
                if not header:
                    header = sorted(list(quality_dict[genome_id].keys()))
                    table.write("{}\n".format("\t".join(header)))

                table.write("{}\n".format("\t".join([quality_dict[genome_id][h] for h in header])))

    return genomes_paths


def process_input_genomes(
    genomes_list: list,
    taxonomic_label: str,
    cluster_id: int,
    db_dir: os.path.abspath,
    tmp_dir: os.path.abspath,
    kmer_len: int,
    min_kmer_occurrences: int=2,
    input_extension: str="fna.gz",
    cluster_prefix: str="MSBT",
    nproc: int=1,
    completeness: float=0.0,
    contamination: float=100.0,
    dereplicate: bool=False,
    similarity: float=1.0,
    flat_structure: bool=False,
    logger: Optional[Logger]=None,
    verbose: bool=False,
) -> List[os.path.abspath]:
    """Process a list of input genomes. Organize, quality-control, and dereplicate genomes.

    Parameters
    ----------
    genomes_list : list
        List of input genome file paths.
    taxonomic_label : str
        Taxonomic label of the input genomes.
    cluster_id : int
        Numerical cluster ID.
    db_dir : os.path.abspath
        Path to the database root folder.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    kmer_len : int
        Length of the kmers.
    min_kmer_occurrences : int, default 2
        Minimum number of kmer occurrences.
    input_extension : str, default "fna.gz"
        File extension of the input files in `genomes_list`.
    cluster_prefix : str, default "MSBT"
        Cluster prefix.
    nproc : int, default 1
        Make the process parallel when possible.
    completeness : float, default 0.0
        Threshold on the completeness.
    contamination : float, default 100.0
        Threshold on the contamination.
    dereplicate : bool, default False
        Enable the dereplication step to get rid of replicated genomes if True.
    similarity : float, default 1.0
        Get rid of genomes according to this similarity threshold in case `dereplicate` is True.
    float_structure : bool, default False
        Taxonomically organize genomes if True.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.

    Returns
    -------
    list
        A list of paths to the genome files.
    """

    # Check whether this set of genomes must be processed
    process_genomes = True

    # Define the cluster folder path
    tax_dir = os.path.join(db_dir, taxonomic_label.replace("|", os.sep)) if not flat_structure else db_dir

    # Search for the metadata file
    metadata_table = os.path.join(tax_dir, "metadata.tsv")

    # Also search for the quality-control table
    quality_table = os.path.join(tax_dir, "quality.tsv")

    if os.path.isfile(metadata_table):
        if not flat_structure:
            os.unlink(metadata_table)

            # Recreate the metadata table
            with open(metadata_table, "w+") as metafile:
                metafile.write("# Cluster ID: {}{}\n".format(cluster_prefix, cluster_id))

        # Genomes have been already processed
        process_genomes = False

    if os.path.isfile(quality_table):
        os.unlink(quality_table)

        # Locate partial quality tables in the tmp folder
        quality_tmp_dir = os.path.join(tmp_dir, "quality", str(cluster_id))
        quality_tables = [str(filepath) for filepath in Path(quality_tmp_dir).glob("run_*.tsv")]

        # Load the quality tables
        quality_dict = dict()

        for table_path in quality_tables:
            with open(table_path) as table:
                header = table.readline().strip()

                for line in table:
                    line = line.strip()
                    
                    if line:
                        line_split = line.split("\t")

                        genome = line_split[header.index("Name")]

                        if "." in input_extension:
                            if genome.lower().endswith(".{}".format(input_extension.lower().split(".")[0])):
                                # Fix the genome name
                                genome = ".".join(genome.split(".")[:-1])

                        quality_dict[genome] = dict()

                        for idx, value in enumerate(line_split):
                            quality_dict[genome][header[idx]] = value

        if quality_dict:
            # Recreate the quality table
            with open(quality_table, "w+") as table:
                header = None

                for genome_id in quality_dict:
                    if not header:
                        header = sorted(list(quality_dict[genome_id].keys()))
                        table.write("{}\n".format("\t".join(header)))

                    table.write("{}\n".format("\t".join([quality_dict[genome_id][h] for h in header])))

        # Genomes have been already processed
        process_genomes = False

    if not process_genomes:
        genomes_dir = os.path.join(tax_dir, "genomes")

        if not flat_structure:
            genomes_dir = os.path.join(tax_dir, "strains", "genomes")

        # Retrieve the genomes
        genomes_paths = [str(filepath) for filepath in Path(genomes_dir).glob("*") if os.path.isfile(filepath)]

    else:
        # Remove data from a previous run
        shutil.rmtree(os.path.join(tmp_dir, "genomes", str(cluster_id)), ignore_errors=True)
        shutil.rmtree(os.path.join(tmp_dir, "howdesbt", str(cluster_id)), ignore_errors=True)
        shutil.rmtree(os.path.join(tmp_dir, "quality", str(cluster_id)), ignore_errors=True)
        shutil.rmtree(tax_dir, ignore_errors=True)

        # Define a partial println function to avoid specifying logger and verbose
        # every time the println function is invoked
        printline = partial(println, logger=logger, verbose=verbose)

        if verbose:
            printline("Processing {}".format(taxonomic_label))

        # Create a temporary folder to store genomes
        tmp_genomes_dir = os.path.join(tmp_dir, "genomes", str(cluster_id))
        os.makedirs(tmp_genomes_dir, exist_ok=True)

        # Iterate over the list of input genome file paths
        genomes = list()

        for genome_path in genomes_list:
            # Symlink input genomes into the temporary folder
            os.symlink(genome_path, os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))
            genomes.append(os.path.join(tmp_genomes_dir, os.path.basename(genome_path)))

        # Take track of the quality dictionary
        quality_dict: Dict[str, Dict[str, str]] = dict()

        # Quality control
        if completeness > 0.0 or contamination < 100.0:
            supported_superkingdoms = ["bacteria", "archaea", "eukaryota", "viruses"]

            superkingdom = taxonomic_label.split("|")[0][3:]

            if superkingdom.lower() in supported_superkingdoms:
                # Try to automatically infer the QC method from the superkingdom
                qc_method = None

                if superkingdom.lower() in ["bacteria", "archaea"]:
                    qc_method = "CheckM2"

                elif superkingdom.lower() == "viruses":
                    qc_method = "CheckV"

                elif superkingdom.lower() == "eukaryota":
                    qc_method = "EukCC"

                if qc_method:
                    before_qc = len(genomes)

                    genomes, quality_dict = quality_control(
                        genomes,
                        cluster_id,
                        tmp_dir,
                        qc_method,
                        input_extension=input_extension,
                        completeness=completeness,
                        contamination=contamination,
                        nproc=nproc,
                    )

                    if before_qc > len(genomes) and verbose:
                        printline("Quality control: excluding {}/{} genomes".format(before_qc - len(genomes), before_qc))

        # Dereplication
        if len(genomes) > 1 and dereplicate:
            before_dereplication = len(genomes)

            genomes = dereplicate_genomes(
                genomes,
                cluster_id,
                tmp_dir,
                kmer_len,
                min_occurrences=min_kmer_occurrences,
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
                printline("No more genomes available for the NCBI tax ID {}".format(cluster_id))
            return genomes

        # Organize genome files
        genomes_paths = organize_data(
            genomes,
            db_dir,
            taxonomic_label,
            cluster_id=cluster_id,
            cluster_prefix=cluster_prefix,
            metadata=None,
            quality_dict=quality_dict,
            flat_structure=flat_structure,
        )

    return genomes_paths


def get_sublist(
    genomes: List[os.path.abspath],
    limit_number: Optional[int]=None,
    limit_percentage: float=100.0,
    flat_structure: bool=False
) -> List[str]:
    """Given a list of genomes, define a sublist with a limited number of genomes taken randomly.

    Parameters
    ----------
    genomes : list
        List of paths to the genome files.
    limit_number : int, optional
        Limit the number of elements in the resulting list to this number.
        It overrides the `limit_percentage` if not None.
    limit_percentage : float, default 100.0
        Limit the size of the resulting list to this percentage computed on the total number of elements in the input list.
    flat_structure : bool, default False
        Taxonomically organize genomes if True.

    Returns
    -------
    list
        The sublist of genomes.
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


def estimate_bf_size_and_howdesbt(
    strains_dir: os.path.abspath,
    tmp_dir: os.path.abspath,
    extension: str="fna.gz",
    kmer_len: int=21,
    min_kmer_occurrences: int=2,
    prefix: str="genomes",
    limit_number: Optional[int]=None,
    limit_percentage: float=100.0,
    increase_filter_size: float=0.0,
    nproc: int=1,
) -> Tuple[os.path.abspath, List[os.path.abspath]]:
    """Wrapper around ntCard and HowDeSBT for building strains SBTs with species-specific bloom filter sizes.
    The structure is always flat here.

    Parameters
    ----------
    strains_dir : os.path.abspath
        Path to the strains folder.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    extension : str, default "fna.gz"
        Input file extension.
    kmer_len : int, default 21
        The length of the kmers.
    min_kmer_occurrences : int, default 2
        Minimum number of occurrences of kmers for estimating the bloom filter size and for building the bloom filters.
    prefix : str, default "genomes"
        Prefix of the ntCard output histogram file.
    limit_number : int, optional
        Maximum number of genomes as input for ntCard.
    limit_percentage : float, default 100.0
        Maximum number of genomes as input for ntCard (in percentage).
    increase_filter_size : float, default 0.0
        Increase the estimated bloom filter size by the specified percentage.
    nproc : int, default 1
        Make ntCard and HowDeSBT parallel.

    Returns
    -------
    tuple
        A tuple with the strains folder path and the list of paths to the selected representative genomes.
    """

    genomes_dir = os.path.join(strains_dir, "genomes")

    if not os.path.isdir(genomes_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), genomes_dir)

    os.makedirs(tmp_dir, exist_ok=True)

    genomes_paths = [str(path) for path in Path(genomes_dir).glob("*.{}".format(extension))]

    # Path to the manifest with the species-specific bloom filter size
    manifest_filepath = os.path.join(strains_dir, "manifest.txt")

    filter_size = None

    if os.path.isfile(manifest_filepath):
        # Load the bloom filter size
        manifest = load_manifest(manifest_filepath)

        filter_size = manifest["filter_size"] if "filter_size" in manifest else None

    if not filter_size:
        # Define a subset of genomes
        genomes_paths_sub = get_sublist(
            genomes_paths,
            limit_number=limit_number,
            limit_percentage=limit_percentage,
            flat_structure=True
        )

        # Estimate the bloom filter size
        filter_size = estimate_bf_size(
            genomes_paths_sub,
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
            prefix=prefix,
            tmp_dir=tmp_dir,
            nproc=nproc,
        )

        # Increment the estimated bloom filter size
        increment = int(math.ceil(filter_size * increase_filter_size / 100.0))
        filter_size += increment

        # Run HowDeSBT
        howdesbt(
            strains_dir,
            extension=extension,
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
            filter_size=filter_size,
            nproc=nproc,
            flat_structure=True,
        )

        # Define the manifest file
        with open(manifest_filepath, "w+") as manifest:
            manifest.write("--min-kmer-occurrences {}\n".format(min_kmer_occurrences))
            manifest.write("--kmer-len {}\n".format(kmer_len))
            manifest.write("--filter-size {}\n".format(filter_size))

    # Select the representative genomes
    selected_genomes = list()

    if len(genomes_paths) <= 3:
        # 3 is the maximum number of selected species
        # as it is also the minimum number of genomes for computing boundaries
        selected_genomes = [
            get_file_info(genome_path, check_supported=False, check_exists=False)[1] for genome_path in genomes_paths
        ]

    else:
        # Check whether the representative genomes have already been selected in a previous run
        species_dir = os.path.dirname(strains_dir)

        if os.path.isdir(os.path.join(species_dir, "genomes")):
            selected_genomes = [os.path.realpath(str(filepath)) for filepath in Path(os.path.join(species_dir, "genomes")).glob("*")
                                if os.path.isfile(str(filepath))]

        if not selected_genomes:
            # Get the bloom filters file paths
            bf_filepaths = [str(path) for path in Path(os.path.join(strains_dir, "filters")).glob("*.bf")]

            # Compute the theta distance between genomes
            bfdistance_theta = bfaction(
                bf_filepaths,
                os.path.join(tmp_dir, "filters"),
                kmer_len,
                min_occurrences=min_kmer_occurrences,
                filter_size=filter_size,
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

            # Also select a genome in the middles of min and max distances
            selected_genomes.append(sorted_genomes[math.ceil(int(len(sorted_genomes) / 2))])

    return strains_dir, selected_genomes


def index(
    db_dir: os.path.abspath,
    input_list: os.path.abspath,
    tmp_dir: os.path.abspath,
    input_extension: str="fna.gz",
    uniform_strand: bool=False,
    cluster_prefix: str="MSBT",
    kmer_len: Optional[int]=None,
    filter_size: Optional[int]=None,
    flat_structure: bool=False,
    estimate_filter_size: bool=True,
    increase_filter_size: float=0.0,
    min_kmer_occurrences: int=2,
    estimate_kmer_size: bool=True,
    limit_estimation_number: Optional[int]=None,
    limit_estimation_percentage: float=100.0,
    limit_kmer_size: int=32,
    completeness: float=0.0,
    contamination: float=100.0,
    dereplicate: bool=False,
    similarity: float=1.0,
    closely_related: bool=False,
    use_representatives: bool=False,
    logger: Optional[Logger]=None,
    verbose: bool=False,
    nproc: int=1,
    jellyfish_threads: int=1,
    parallel: int=1,
) -> None:
    """Build the database baseline.

    Parameters
    ----------
    db_dir : os.path.abspath
        Path to the database root folder.
    input_list : os.path.abspath
        Path to the file with a list of input genome paths.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    input_extension : str, default "fna.gz"
        File extension of the input files whose paths are defined into the `input_list` file.
    uniform_strand : bool, default False
        Preprocess the input fasta sequences to make them all on the same strand.
    cluster_prefix : str, default "MSBT",
        Prefix of clusters.
    kmer_len : int, optional
        Length of the kmers.
    filter_size : int, optional
        Size of the bloom filters.
    flat_structure : bool, default False
        Taxonomically organize genomes if True.
    estimate_filter_size : bool, default True
        Run ntCard to estimate the most appropriate bloom filter size.
    increase_filter_size : float, default 0.0
        Increase the estimated bloom filter size by the specified percentage.
    min_kmer_occurrences : int, default 2
        Minimum number of occurrences of kmers for estimating the bloom filter size and for building bloom filters.
    estimate_kmer_size : bool, default True
        Run kitsune to estimate the best kmer size.
    limit_estimation_number : int, optional
        Number of genomes per group as input for kitsune and ntCard.
    limit_estimation_percentage : float, default 100.0
        Percentage on the total number of genomes per group as input for kitsune and ntCard.
    limit_kmer_size : int, default 32
        Maximum kmer size for kitsune kopt.
    completeness : float, default 0.0
        Threshold on the completeness.
    contamination : float, default 100.0
        Threshold on the contamination.
    dereplicate : bool, default False
        Enable the dereplication step to get rid of replicated genomes if True.
    similarity : float, default 1.0
        Get rid of genomes according to this threshold in case the dereplication step is enabled.
    closely_related : bool, default False
        Used to process closesly related genomes with kitsune kopt.
    use_representatives : bool, default False
        Use 3 representative genomes for each species.
    logger : logging.Logger, optional
        The Logger object.
    verbose : bool, default False
        Print messages on the stdout if True.
    nproc : int, default 1
        Make the process parallel when possible.
    jellyfish_threads : int, default 1
        Maximum number of threads to make Jellyfish parallel with kitsune.
    parallel : int, default 1
        Maximum number of processes to process each NCBI tax ID in parallel.

    Raises
    ------
    Exception
        - If there are no input genomes;
        - If there are not enough genomes for estimating the optimal kmer size with kitsune.
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Load the list of input genomes and eventually their taxonomic labels
    taxonomy2genomes = load_input_table(input_list, input_extension=input_extension)

    if uniform_strand:
        genome_paths = [genome_path for taxonomy in taxonomy2genomes for genome_path in taxonomy2genomes[taxonomy]]

        stranded_genome_paths = dict()

        with mp.Pool(processes=nproc) as pool:
            jobs = [
                pool.apply_async(strand, args=(genome_path, os.path.join(tmp_dir, "uniform_strand"),))
                    for genome_path in genome_paths
            ]

            for job in jobs:
                stranded_genome_path = job.get()

                # Get the genome name from file path
                _, stranded_genome_name, _, _ = get_file_info(stranded_genome_path)

                stranded_genome_paths[stranded_genome_name] = stranded_genome_path

        preprocessed_taxonomy2genomes = dict()

        for taxonomy in taxonomy2genomes:
            preprocessed_genomes = list()

            for genome in taxonomy2genomes[taxonomy]:
                # Get the genome name from file path
                _, genome_name, _, _ = get_file_info(genome)

                preprocessed_genomes.append(stranded_genome_paths[genome_name])

            preprocessed_taxonomy2genomes[taxonomy] = preprocessed_genomes

        taxonomy2genomes = preprocessed_taxonomy2genomes

    # Force flat structure in case of genomes with no taxonomic label
    if "NA" in taxonomy2genomes:
        flat_structure = True

    # Build a partial function around process_input_genomes
    process_partial = partial(
        process_input_genomes,
        db_dir=db_dir,
        tmp_dir=tmp_dir,
        kmer_len=kmer_len,
        min_kmer_occurrences=min_kmer_occurrences,
        input_extension=input_extension,
        cluster_prefix=cluster_prefix,
        nproc=nproc,
        completeness=completeness,
        contamination=contamination,
        dereplicate=dereplicate,
        similarity=similarity,
        flat_structure=flat_structure,
        logger=logger,
        verbose=False,
    )

    printline("Processing clusters")

    # Define a cluster counter
    # This is also the size of the progress bar
    clusters_counter = len(taxonomy2genomes)

    # Take track of all the genomes paths
    genomes_paths = list()

    with mp.Pool(processes=parallel) as pool, tqdm.tqdm(total=clusters_counter, disable=(not verbose)) as pbar:
        # Wrapper around the update function of tqdm
        def progress(*args):
            pbar.update()

        # Process input genomes
        jobs = [
            pool.apply_async(
                process_partial,
                args=(taxonomy2genomes[taxonomy], taxonomy, pos + 1),
                callback=progress,
            )
            for pos, taxonomy in enumerate(taxonomy2genomes)
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
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
            prefix="genomes",
            tmp_dir=tmp_dir,
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
            extension=input_extension,
            kmer_len=kmer_len,
            min_occurrences=min_kmer_occurrences,
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
                folders = [str(path) for path in Path(db_dir).glob("**/{}__*".format(level[0])) if os.path.isdir(str(path))]
                printline("Running HowDeSBT at the {} level ({} clusters)".format(level, len(folders)))

                if level == "species":
                    if use_representatives:
                        estimate_bf_size_and_howdesbt_partial = partial(
                            estimate_bf_size_and_howdesbt,
                            extension=input_extension,
                            kmer_len=kmer_len,
                            min_kmer_occurrences=min_kmer_occurrences,
                            prefix="genomes",
                            limit_number=limit_estimation_number,
                            limit_percentage=limit_estimation_percentage,
                            increase_filter_size=increase_filter_size,
                            nproc=nproc,
                        )

                        # Establish a species-specific bloom filter size
                        # Run howdesbt in flat mode
                        # Search for genomes that minimize and maximize the genetic distance versus all the other genomes in the same cluster
                        # Use only these genomes to build the species tree
                        with mp.Pool(processes=parallel) as strains_pool, tqdm.tqdm(total=len(folders), disable=(not verbose)) as pbar:
                            # Wrapper around the update function of tqdm
                            def progress(*args):
                                pbar.update()

                            # Process strains
                            jobs = [
                                strains_pool.apply_async(
                                    estimate_bf_size_and_howdesbt_partial,
                                    args=(
                                        os.path.join(species_dir, "strains"),
                                        os.path.join(species_dir, "strains", "tmp"),
                                    ),
                                    callback=progress
                                )
                                for species_dir in folders
                            ]

                            for job in jobs:
                                strains_dir, selected_genomes = job.get()
                                species_dir = os.sep.join(strains_dir.split(os.sep)[:-1])

                                # Populate the genomes folder at the species level with the selected genomes
                                for genome in selected_genomes:
                                    genome_link = os.path.join(species_dir, "genomes", "{}.{}".format(genome, input_extension))

                                    if not os.path.islink(genome_link):
                                        os.symlink(
                                            os.path.join(strains_dir, "genomes", "{}.{}".format(genome, input_extension)),
                                            genome_link
                                        )

                    else:
                        for species_dir in folders:
                            # Use all the genomes under this species
                            for genome_path in Path(os.path.join(species_dir, "strains", "genomes")).glob("*.{}".format(input_extension)):
                                genome_link = os.path.join(species_dir, "genomes", os.path.basename(genome_path))

                                if not os.path.islink(genome_link):
                                    os.symlink(
                                        genome_path,
                                        genome_link
                                    )

                with mp.Pool(processes=parallel) as pool, tqdm.tqdm(total=len(folders), disable=(not verbose)) as pbar:
                    # Wrapper around the update function of tqdm
                    def progress(*args):
                        pbar.update()

                    # Process clusters under a specific taxonomic level
                    # Run HowDeSBT only in case the bloom filter representation of the current cluster does not exist
                    jobs = [
                        pool.apply_async(howdesbt_partial, args=(level_dir,), callback=progress)
                            for level_dir in folders if not os.path.isfile(os.path.join(level_dir, "{}.bf".format(os.path.basename(level_dir))))
                    ]

                    for job in jobs:
                        job.get()

        # Also run HowDeSBT on the database root folder to build
        # the bloom filter representation of the superkingdom
        if not os.path.isfile(os.path.join(db_dir, "{}.bf".format(os.path.basename(db_dir)))):
            if not flat_structure:
                printline("Building the database root bloom filter with HowDeSBT")

            howdesbt_partial(db_dir)

        # The howdesbt function automatically set the current working directory to
        # the index folder of the taxonomic labels
        # Come back to the original folder
        os.chdir(current_folder)


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
        # Try to resume the index process
        # Read the last command line from the sh script
        argv = [line.strip() for line in open(args.resume).readlines() if line.strip() and not line.startswith("#")][-1].split()

        # Reload command line parameters
        if argv[0] == TOOL_ID:
            # It was a direct call to the index module
            args = read_params(argv[1:])

        else:
            # The index module has been run through the controller
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

        if args.flat_structure and resume:
            raise Exception("The --resume option cannot be used in conjunction with --flat-structure")

    # Create the database folder
    os.makedirs(args.db_dir, exist_ok=True)

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

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
    if not args.flat_structure:
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

        if superkingdoms:
            for superkingdom in superkingdoms:
                if os.path.isdir(os.path.join(args.db_dir, "k__{}".format(superkingdom))) and not resume:
                    raise Exception(
                        (
                            "An indexed version of the {} superkingdom already exists in the database!\n"
                            "Please use the update module to add new genomes"
                        ).format(superkingdom)
                    )

    # Define the database manifest file
    manifest_filepath = os.path.join(args.db_dir, "manifest.txt")

    if os.path.isfile(manifest_filepath):
        # Load and compare --kmer-len, --filter-size, --min-kmer-occurrences, and --use-representatives
        manifest = load_manifest(manifest_filepath)

        for arg in ["kmer_len", "filter_size", "min_kmer_occurrences", "use_representatives"]:
            if not hasattr(args, arg):
                setattr(args, arg, manifest[arg])

            elif getattr(args, arg) != manifest[arg]:
                raise ValueError(
                    "The specified \"--{}\" is not compatible with the selected database".format(
                        arg.replace("_", "-")
                    )
                )

    else:
        # Initialize manifest file
        with open(manifest_filepath, "w+") as manifest:
            if args.kmer_len:
                # This will be added later if it does not exist
                manifest.write("--kmer-len {}\n".format(args.kmer_len))

            if args.filter_size:
                # This will be added later if it does not exist
                manifest.write("--filter-size {}\n".format(args.filter_size))

            # Add --min-kmer-occurrences
            manifest.write("--min-kmer-occurrences {}\n".format(args.min_kmer_occurrences))

            # Add --use-representatives
            manifest.write("--use-representatives {}\n".format(args.use_representatives))

    if not resume:
        # Build a sh script with the command line used to launch the index module
        build_sh(sys.argv, TOOL_ID, args.db_dir)

    t0 = time.time()

    index(
        args.db_dir,
        args.input_list,
        args.tmp_dir,
        input_extension=args.extension,
        uniform_strand=args.uniform_strand,
        cluster_prefix=args.cluster_prefix,
        kmer_len=args.kmer_len,
        filter_size=args.filter_size,
        flat_structure=args.flat_structure,
        estimate_filter_size=args.estimate_filter_size,
        increase_filter_size=args.increase_filter_size,
        min_kmer_occurrences=args.min_kmer_occurrences,
        estimate_kmer_size=args.estimate_kmer_size,
        limit_estimation_number=args.limit_estimation_number,
        limit_estimation_percentage=args.limit_estimation_percentage,
        limit_kmer_size=args.limit_kmer_size,
        completeness=args.completeness,
        contamination=args.contamination,
        dereplicate=args.dereplicate,
        similarity=args.similarity,
        closely_related=args.closely_related,
        use_representatives=args.use_representatives,
        logger=logger,
        verbose=args.verbose,
        nproc=args.nproc,
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
