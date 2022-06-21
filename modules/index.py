#!/usr/bin/env python3

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 17, 2022"

import sys, os, time, tarfile, gzip, re, shutil, numpy, tqdm
import argparse as ap
import multiprocessing as mp
from pathlib import Path
from functools import partial
from logging import Logger
from typing import List, Tuple

try:
    # Load utility functions
    from utils import checkm, download, filter_checkm_tables, filter_genomes, howdesbt, init_logger, it_exists, kmtricks_matrix, number, println, run
except:
    pass

# Define the module name
TOOL_ID = "index"

# Define the list of dependencies
DEPENDENCIES = ["checkm", "howdesbt", "kmtricks", "ncbi-genome-download", "ncbitax2lin", "ntcard", "wget"]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",   # Database folder path
    "--log",      # Path to the log file
    "--tmp-dir"   # Temporary folder path
]

# Define the url to the NCBI taxdump
TAXDUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

def read_params():
    p = ap.ArgumentParser(prog=TOOL_ID,
                          description=("Retrieve complete reference genomes from the NCBI GenBank and build a "
                                       "Sequence Bloom Tree for each taxonomic level"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--completeness",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 0.0,
                    help = "Input genomes must have a minimum completeness percentage before being processed and added to the database" )
    p.add_argument( "--contamination",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 100.0,
                    help = "Input genomes must have a maximum contamination percentage before being processed and added to the database" )
    p.add_argument( "--cleanup",
                    action = "store_true",
                    default = False,
                    help = "Remove temporary data at the end of the pipeline" )
    p.add_argument( "--db-dir",
                    type = os.path.abspath,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--dereplicate",
                    action = "store_true",
                    default = False,
                    help = "Enable the dereplication of genomes" )
    p.add_argument( "--estimate-filter-size",
                    action = "store_true",
                    default = False,
                    dest = "estimate_filter_size",
                    help = "Automatically estimate the best bloom filter size" )
    p.add_argument( "--filter-size",
                    type = number(int, minv=10000),
                    dest = "filter_size",
                    help = "This is the size of the bloom filters" )
    p.add_argument( "--increase-filter-size",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 0.0,
                    dest = "increase_filter_size",
                    help = ("Increase the estimated filter size by the specified percentage. "
                            "This is used in conjunction with the --estimate_filter_size argument only. "
                            "It is highly recommended to increase the filter size by a good percentage in case you are planning to update the index with new genomes") )
    p.add_argument( "--kingdom",
                    type = str,
                    required = True,
                    choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
                    help = "Consider genomes whose lineage belongs to a specific kingdom" )
    p.add_argument( "--kmer-len",
                    type = number(int, minv=6),
                    required = True,
                    dest = "kmer_len",
                    help = "This is the length of the kmers used for building bloom filters" )
    p.add_argument( "--limit-genomes",
                    type = number(int, minv=1),
                    default = numpy.Inf,
                    dest = "limit_genomes",
                    help = ("Limit the number of genomes per species. "
                            "This will remove the exceeding number of genomes randomly to cut the overall number of genomes per species to this number. "
                            "The number of genomes per species is not limited by default") )
    p.add_argument( "--log",
                    type = os.path.abspath,
                    help = "Path to the log file" )
    p.add_argument( "--max-genomes",
                    type = number(int, minv=1),
                    default = numpy.Inf,
                    dest = "max_genomes",
                    help = ("Consider species with this number of genomes at most. "
                            "There is not a maximum number of genomes per species by default") )
    p.add_argument( "--min-genomes",
                    type = number(int, minv=1),
                    default = 1,
                    dest = "min_genomes",
                    help = ("Consider species with a minimum number of genomes. "
                            "There is not a minimum number of genomes per species by default") )
    p.add_argument( "--nproc",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "This argument refers to the number of processors used for parallelizing the pipeline when possible" )
    p.add_argument( "--parallel",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "Maximum number of processors to process each NCBI tax ID in parallel" )
    p.add_argument( "--pplacer-threads",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    dest = "pplacer_threads",
                    help = "Maximum number of threads for pplacer. This is required to maximise the CheckM performances" )
    p.add_argument( "--similarity",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 100.0,
                    help = ("Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. "
                            "This is used exclusively in conjunction with the --dereplicate argument") )
    p.add_argument( "--tmp-dir",
                    type = os.path.abspath,
                    required = True,
                    dest = "tmp_dir",
                    help = "Path to the folder for storing temporary data" )
    p.add_argument( "--verbose",
                    action = "store_true",
                    default = False,
                    help = "Print results on screen" )
    p.add_argument( "-v",
                    "--version",
                    action = "version",
                    version = "{} version {} ({})".format(TOOL_ID, __version__, __date__),
                    help = "Print the current {} version and exit".format(TOOL_ID) )
    return p.parse_args()

def count_genomes(tax_id: str, kingdom: str, out_dir: str) -> int:
    """
    Count the number of genomes classified as the provided NCBI tax ID

    :param tax_id:      NCBI tax ID
    :param kingdom:     Search for a specific kingdom
    :param out_dir:     Folder path for storing temporary data
    :return:            The genomes count
    """

    # Define the ncbi-genome-download command line
    ncbi_genome_download = ["ncbi-genome-download", "--assembly-levels", "complete",              # Search for complete genomes only
                                                    "--section", "genbank",                       # Search for genomes in GenBank
                                                    "--dry-run",                                  # Retrieve the list of genomes
                                                    "--species-taxids", tax_id, kingdom.lower()]  # Search for a specific NCBI tax ID

    # First check how many genomes belong to the specific tax ID
    # Define the log file path with the list of genomes retrieved by ncbi-genome-download
    genomes_list_filepath = os.path.join(out_dir, "genomes.txt")
    
    # Run ncbi-genome-download
    with open(genomes_list_filepath, "w+") as genomes_list:
        try:
            run(ncbi_genome_download, stdout=genomes_list, stderr=genomes_list)
        except:
            # This is required in case there are no genomes matching the search criteria
            pass
    
    # Count how many genomes belongs to the current tax ID
    genomes_counter = 0
    with open(genomes_list_filepath) as genomes_list:
        for l in genomes_list:
            l = l.strip()
            if l:
                if l.startswith("GCA_"):
                    genomes_counter += 1
    
    return genomes_counter

def download_taxdump(taxdump_url: str, folder_path: str) -> Tuple[str, str]:
    """
    Download and extract the NCBI taxdump tarball

    :param taxdump_url:     URL to the NCBI taxdump
    :param folder_path:     Path to the folder in which the taxdump tarball will be unpacked
    :return:                The nodes.dmp and names.dmp file paths
    """

    taxdump = download(taxdump_url, folder_path)
    
    # Raise an exception in case the taxdump.tar.gz file does not exist
    if not os.path.exists(taxdump):
        raise Exception("Unable to retrieve data from remote location\n{}".format(taxdump_url))

    # Create the taxdump folder in the temporary directory
    taxdump_dir = os.path.join(folder_path, "taxdump")
    os.makedirs(taxdump_dir, exist_ok=True)

    # Decompress the archive
    with tarfile.open(taxdump, "r:gz") as tar:
        tar.extractall(taxdump_dir)
    
    return os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp")

def estimate_bf_size(genomes_filepath: str, kmer_len: int, prefix: str, nproc: int=1) -> int:
    """
    Estimate the bloom filter size with ntCard

    :param genomes_filepath:    Path to the file with the list of paths to the genome files
    :param kmer_len:            Length of the kmers
    :param prefix:              Prefix of the output histogram file
    :param nproc:               Make it parallel
    :return:                    The estimated bloom filter size
    """

    # Estimate the bloom filter size with ntCard
    run(["ntcard", "--kmer={}".format(kmer_len),    # Kmer length
                   "--threads={}".format(nproc),    # Maximum number of threads
                   "--pref={}".format(prefix),      # Prefix of the output histogram file
                   "@{}".format(genomes_filepath)], # Input file with the list of genomes paths
        silence=True)

    F0 = 0
    f1 = 0

    # Read the ntCard output hist file
    with open(os.path.join(tmp_dir, "genomes_k{}.hist".format(kmer_len))) as histfile:
        F0_found = False
        f1_found = False
        for line in histfile:
            line = line.strip()
            if line:
                if line.startswith("F0"):
                    F0 = int(line.split()[-1])
                    F0_found = True
                elif line.startswith("1"):
                    f1 = int(line.split()[-1])
                    f1_found = True
                if F0_found and f1_found:
                    break
    
    # Compute the bloom filter size
    return F0-f1

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

def load_taxa(ncbitax2lin_table: str, kingdom: str=None, dump: str=None) -> Tuple[list, list]:
    """
    Load the ncbitax2lin output table

    :param ncbitax2lin_table:   Path to the output table produced by ncbitax2lin
    :param kingdom:             Consider a specific kingdom only
    :param dump:                Path to the output table
    :return:                    The lists of NCBI tax IDs and full taxonomic labels
    """
    
    # Take track of NCBI tax IDs and full taxonomic labels
    taxa_map = dict()

    with gzip.open(ncbitax2lin_table, "rt") as ncbi_table:
        next(ncbi_table)
        for line in ncbi_table:
            line = line.strip()
            if line:
                line_split = line.split(",")

                # Check whether the current taxonomy must be processed                
                skip = True
                if not kingdom:
                    if line_split[1].strip():
                        skip = False
                elif line_split[1] == kingdom:
                    skip = False

                if not skip:
                    # Build the current full taxonomic label
                    label = "k__{}|p__{}|c__{}|o__{}|f__{}|g__{}|s__{}".format(
                        kingdom,                                    # Kingdom
                        level_name(line_split[2], kingdom),         # Phylum
                        level_name(line_split[3], line_split[2]),   # Class
                        level_name(line_split[4], line_split[3]),   # Order
                        level_name(line_split[5], line_split[4]),   # Family
                        level_name(line_split[6], line_split[5]),   # Genus
                        level_name(line_split[7], line_split[6])    # Species
                    )

                    # Exclude unclassified taxonomic labels
                    if "unclassified" not in label:
                        if label not in taxa_map:
                            taxa_map[label] = list()
                        taxa_map[label].append(int(line_split[0]))
    
    # Use a single tax ID for each label
    tax_labels = [tax_label for tax_label in taxa_map]
    tax_ids = [str(min(taxa_map[tax_label])) for tax_label in taxa_map]
    
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
    run(["ncbitax2lin", "--nodes-file", nodes,
                        "--names-file", names,
                        "--output", ncbitax2lin_table],
        silence=True)

    # Raise an exception in case the ncbi_lineages.csv.gz file does not exist
    if not os.path.exists(ncbitax2lin_table):
        raise Exception("An error has occurred while running ncbitax2lin")
    
    return ncbitax2lin_table

def process_tax_id(tax_id: str, tax_label: str, kingdom: str, db_dir: str, tmp_dir: str, limit_genomes: float=numpy.Inf, 
                   max_genomes: float=numpy.Inf, min_genomes: float=1.0, nproc: int=1, pplacer_threads: int=1, completeness: float=0.0, 
                   contamination: float=100.0, dereplicate: bool=False, similarity: float=100.0, logger: Logger=None, verbose: bool=False) -> List[str]:
    """
    Process a specific NCBI tax ID
    Download, quality-control, and dereplicate genomes

    :param tax_id:              NCBI tax ID of a species
    :param tax_label:           Full taxonomic label
    :param kingdom:             Retrieve genomes that belong to a specific kingdom
    :param db_dir:              Path to the database root folder
    :param tmp_dir:             Path to the temporary folder
    :param limit_genomes:       Limit the number of genomes per species
    :param max_genomes:         Consider species with this number of genomes at most
    :param min_genomes:         Consider species with a minimum number of genomes
    :param nproc:               Make the process parallel when possible
    :param pplacer_threads:     Maximum number of threads to make pplacer parallel with CheckM
    :param completeness:        Threshold on the CheckM completeness
    :param contamination:       Threshold on the CheckM contamination
    :param dereplicate:         Enable the dereplication step to get rid of replicated genomes
    :param similarity:          Get rid of genomes according to this threshold in case the dereplication step is enabled
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

    # Count how many genomes belong to the same NCBI tax ID
    genomes_counter = count_genomes(tax_id, kingdom, tmp_genomes_dir)

    # Take track of the paths to the genome files
    genomes_paths = list()

    if genomes_counter > 0:
        if verbose:
            printline("Processing NCBI tax ID {} with {} genomes".format(tax_id, genomes_counter))

        # Check whether the number of genomes falls in the --min-genomes --max-genomes interval            
        if genomes_counter < min_genomes or genomes_counter > max_genomes:
            if verbose:
                printline("WARNING: the number of genomes does not respect the minimum and maximum number of genomes allowed")
            return list()

        # Check whether the number of genomes is greater than the number provided with --limit-genomes
        if genomes_counter > limit_genomes and verbose:
            printline("WARNING: the number of genomes per species is limited to {} at most".format(limit_genomes))

        # Retrieve the genome files from NCBI
        genomes, metadata = retrieve_genomes(tax_id, kingdom, tmp_genomes_dir, limit_genomes=limit_genomes, nproc=nproc, retries=1)
        if verbose:
            printline("Retrieved {}/{} genomes".format(len(genomes), genomes_counter))

        # Take track of the paths to the CheckM output tables
        checkm_tables = []

        # Quality control
        if genomes and (completeness > 0.0 or contamination < 100.0):
            # Define the CheckM temporary folder
            checkm_tmp_dir = os.path.join(tmp_dir, "checkm", tax_id)
            os.makedirs(checkm_tmp_dir, exist_ok=True)
            # Run CheckM on the current set of genomes
            checkm_tables = checkm(genomes, checkm_tmp_dir, file_extension="fna.gz", 
                                                            nproc=nproc, 
                                                            pplacer_threads=pplacer_threads)
            # Filter genomes according to the input --completeness and --contamination thresholds
            genome_ids = filter_checkm_tables(checkm_tables, completeness=0.0, contamination=100.0)
            
            before_qc = len(genomes)
            genomes = [os.path.join(tmp_genomes_dir, "{}.fna.gz".format(genome_id)) for genome_id in genome_ids]

            if before_qc > len(genomes) and verbose:
                printline("Quality control: excluding {}/{} genomes".format(before_qc-len(genomes), before_qc))

        # Dereplication
        if len(genomes) > 1 and dereplicate:
            # Define the kmtricks temporary folder
            kmtricks_tmp_dir = os.path.join(tmp_dir, "kmtricks", tax_id)
            os.makedirs(kmtricks_tmp_dir, exist_ok=True)

            # Create a fof file with the list of genomes
            genomes_fof_filepath = os.path.join(kmtricks_tmp_dir, "genomes.fof")
            ordered_genome_names = list()
            with open(genomes_fof_filepath, "w+") as genomes_fof:
                for genome_path in genomes:
                    genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                    if genome_path.endswith(".gz"):
                        genome_name = os.path.splitext(genome_name)[0]
                    genomes_fof.write("{} : {}\n".format(genome_name, genome_path))
                    ordered_genome_names.append(genome_name)

            # Run kmtricks to produce the kmer matrix
            output_table = os.path.join(kmtricks_tmp_dir, "matrix.txt")
            kmtricks_matrix(genomes_fof_filepath, kmtricks_tmp_dir, nproc, output_table)

            # Add the header line to the kmer matrix
            with open(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), "w+") as file1:
                file1.write("#kmer {}\n".format(" ".join(ordered_genome_names)))
                with open(output_table) as file2:
                    for line in file2:
                        line = line.strip()
                        if line:
                            file1.write("{}\n".format(line))
            
            # Get rid of the matrix withour header line
            os.unlink(output_table)
            shutil.move(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), output_table)

            # Filter genomes according to their percentage of common kmers defined with --similarity
            filtered_genomes_filepath = os.path.join(tmp_dir, "kmtricks", tax_id, "filtered.txt")
            filter_genomes(output_table, filtered_genomes_filepath, similarity=similarity)
            
            before_dereplication = len(genomes)

            # Iterate over the list of dereplicated genome and rebuild the list of genome file paths
            with open(filtered_genomes_filepath) as filtered_genomes:
                for l in filtered_genomes:
                    l = l.strip()
                    if l:
                        gpath = os.path.join(tmp_genomes_dir,"{}.fna.gz".format(l))
                        if gpath in genomes:
                            genomes.remove(gpath)
            
            if before_dereplication > len(genomes) and verbose:
                printline("Dereplication: excluding {}/{} genomes".format(before_dereplication-len(genomes), before_dereplication))
        
        # Check whether no genomes survived the quality control and the dereplication steps
        if not genomes:
            if verbose:
                printline("No more genomes available for the NCBI tax ID {}".format(tax_id))
            return genomes

        # In case at least one genome survived both the quality control and dereplication steps
        # Define the taxonomy folder in database
        tax_dir = os.path.join(db_dir, tax_label.replace("|", os.sep))
        genomes_dir = os.path.join(tax_dir, "genomes")
        os.makedirs(genomes_dir, exist_ok=True)

        with open(os.path.join(tax_dir, "references.txt"), "w+") as refsfile:
            # Move the processed genomes to the taxonomy folder
            for genome_path in genomes:
                shutil.move(genome_path, genomes_dir)
                
                # Also take track of the genome names in the references.txt file
                genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                if genome_path.endswith(".gz"):
                    genome_name = os.path.splitext(genome_name)[0]
                refsfile.write("{}\n".format(genome_name))

                # Add the genome to the full list of genomes in database
                genomes_paths.append(os.path.join(genomes_dir, os.path.basename(genome_path)))
        
        # Move the metadata table to the taxonomy folder
        if os.path.exists(metadata):
            shutil.move(metadata, tax_dir)
        
        if checkm_tables:
            # Also merge the CheckM output tables and move the result to the taxonomy folder
            with open(os.path.join(tax_dir, "checkm.tsv"), "w+") as table:
                header = True
                for table_path in checkm_tables:
                    with open(table_path) as partial_table:
                        line_count = 0
                        for l in partial_table:
                            l = l.strip()
                            if l:
                                if line_count == 0:
                                    if header:
                                        table.write("{}\n".format(l))
                                        header = False
                                else:
                                    table.write("{}\n".format(l))
                                line_count += 1
    
    # Return the list of paths to the genome files
    return genomes_paths

def retrieve_genomes(tax_id: str, kingdom: str, out_dir: str, limit_genomes: float=numpy.Inf, nproc: int=1, retries: int=1) -> Tuple[List[str], str]:
    """
    Retrieve genomes from NCBI under a specific tax ID

    :param tax_id:          NCBI tax ID
    :param kingdom:         Search for a specific kingdom
    :param out_dir:         Folder path for storing temporary data
    :param limit_genomes:   Limit the number of genomes
    :param nproc:           Make it parallel
    :param retries:         Retry downloading genomes
    :return:                List of paths to the genome files and path to the metadata table
    """

    ncbi_genome_download = ["ncbi-genome-download", "--assembly-levels", "complete",                            # Search for complete genomes only
                                                    "--section", "genbank",                                     # Search for genomes in GenBank
                                                    "--output-folder", out_dir, "--flat-output",                # Add the --output-folder
                                                    "--formats", "fasta",                                       # Retrieve the fna.gz files
                                                    "--metadata-table", os.path.join(out_dir, "metadata.tsv"),  # Also retrieve genomes metadata
                                                    "--parallel", str(nproc),                                   # Make it parallel
                                                    "--retries", str(retries),                                  # Retry downloading genomes
                                                    "--species-taxids", tax_id, kingdom.lower()]                # Search for a specific NCBI tax ID

    # Run ncbi-genome-download to retrieve all the complete genomes 
    # related to the current tax ID from NCBI GenBank
    run(ncbi_genome_download, silence=True)
    
    # Define the list of paths to the genome files
    genomes = list()
    for genome_path in Path(out_dir).glob("*.fna.gz"):
        # Retrieve the genome name
        genome_name = os.path.splitext(os.path.splitext(os.path.basename(str(genome_path)))[0])[0].split(".")[0]
        # Rename the genome file
        shutil.move(str(genome_path), os.path.join(out_dir, "{}.fna.gz".format(genome_name)))
        # Add the genome path to the list of genome paths
        genomes.append(os.path.join(out_dir, "{}.fna.gz".format(genome_name)))
    
    if len(genomes) > limit_genomes:
        # ncbi-genome-download does not allow to limit the number of genomes that must be retrieved
        # Selecte a random set of at most "--limit-genomes" genomes
        selected = list()
        random_choice = numpy.random.choice(len(genomes), size=limit_genomes, replace=False)
        for pos, value in enumerate(random_choice):
            if value == 0:
                os.unlink(genomes[pos])
            else:
                selected.append(genomes[pos])
        genomes = selected
    
    # Return the list of paths to the genome files
    return genomes, os.path.join(out_dir, "metadata.tsv")

def index(db_dir: str, kingdom: str, tmp_dir: str, kmer_len: int, filter_size: int, limit_genomes: float=numpy.Inf,
          max_genomes: float=numpy.Inf, min_genomes: float=1.0, estimate_filter_size: bool=False, increase_filter_size: float=0.0, 
          completeness: float=0.0, contamination: float=100.0, dereplicate: bool=False, similarity: float=100.0, logger: Logger=None, 
          verbose: bool=False, nproc: int=1, pplacer_threads: int=1, parallel: int=1) -> None:
    """
    Build the database baseline

    :param db_dir:                  Path to the database root folder
    :param kingdom:                 Retrieve genomes that belong to a specific kingdom
    :param tmp_dir:                 Path to the temporary folder
    :param kmer_len:                Length of the kmers
    :param filter_size:             Size of the bloom filters
    :param limit_genomes:           Limit the number of genomes per species
    :param max_genomes:             Consider species with this number of genomes at most
    :param min_genomes:             Consider species with a minimum number of genomes
    :param estimate_filter_size:    Run ntCard to estimate the most appropriate bloom filter size
    :param increase_filter_size:    Increase the estimated bloom filter size by the specified percentage
    :param completeness:            Threshold on the CheckM completeness
    :param contamination:           Threshold on the CheckM contamination
    :param dereplicate:             Enable the dereplication step to get rid of replicated genomes
    :param similarity:              Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param logger:                  Logger object
    :param verbose:                 Print messages on screen
    :param nproc:                   Make the process parallel when possible
    :param pplacer_threads:         Maximum number of threads to make pplacer parallel with CheckM
    :param parallel:                Maximum number of processors to process each NCBI tax ID in parallel
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Download the NCBI taxdump
    printline("Downloading taxonomy dump from NCBI")
    # Recover already processed nodes.dmp and names.dmp
    nodes_dmp = os.path.join(tmp_dir, "taxdump", "nodes.dmp")
    names_dmp = os.path.join(tmp_dir, "taxdump", "names.dmp")
    if not os.path.exists(nodes_dmp) or not os.path.exists(names_dmp):
        nodes_dmp, names_dmp = download_taxdump(TAXDUMP_URL, tmp_dir)

    # Run ncbitax2lin to extract lineages
    printline("Exporting lineages")
    ncbitax2lin_table = os.path.join(tmp_dir, "ncbi_lineages.csv.gz")
    if not os.path.exists(ncbitax2lin_table):
        ncbitax2lin_table = ncbitax2lin(nodes_dmp, names_dmp, tmp_dir)

    # Build a mapping between tax IDs and full taxonomic labels
    # Take track of the taxonomic labels to remove duplicates
    # Taxonomy IDs are already sorted in ascending order in the ncbitax2lin output table
    printline("Building species tax ID to full taxonomy mapping")
    dump_exists = os.path.exists(os.path.join(tmp_dir, "taxa.tsv"))
    tax_ids, tax_labels = load_taxa(ncbitax2lin_table, kingdom=kingdom, dump=os.path.join(tmp_dir, "taxa.tsv") if not dump_exists else None)

    # Take track of all the genomes paths
    genomes_paths = list()
    
    printline("Retrieving genomes from NCBI GenBank")
    if completeness > 0.0 or contamination < 100.0:
        printline("Enabling quality control with the following thresholds: completeness={}, contamination={}"
                  .format(completeness, contamination))
    if dereplicate:
        printline("Enabling the dereplication of the input genomes with the following threshold: similarity={}"
                  .format(similarity))

    # Build a partial function around process_tax_id
    process_tax_id_partial = partial(process_tax_id, kingdom=kingdom, db_dir=db_dir, tmp_dir=tmp_dir, limit_genomes=limit_genomes, 
                                                     max_genomes=max_genomes, min_genomes=min_genomes, nproc=nproc, pplacer_threads=pplacer_threads, 
                                                     completeness=completeness, contamination=contamination, dereplicate=dereplicate, similarity=similarity, 
                                                     logger=logger, verbose=False)

    # Initialise the progress bar
    pbar = tqdm.tqdm(total=len(tax_ids), disable=(not verbose))

    with mp.Pool(processes=parallel) as pool:
        def update_bar(*args):
            # Update the progress bar
            pbar.update()
            return
        
        # Process the NCBI tax IDs
        jobs = [pool.apply_async(process_tax_id_partial, args=(tax_id, tax_labels[pos], ), callback=update_bar) 
                    for pos, tax_id in enumerate(tax_ids)]

        # Get results from jobs
        for job in jobs:
            genomes_paths.extend(job.get())
    
    # Close the progress bar
    pbar.close()

    # Check whether the bloom filter size must be estimated
    if genomes_paths and estimate_filter_size and not filter_size:
        # Dump the global list of genome path to file
        genomes_counter = 0
        with open(os.path.join(tmp_dir, "genomes.txt"), "w+") as genomesfile:
            for path in genomes_paths:
                if os.path.exists(path):
                    genomesfile.write("{}\n".format(path))
                    genomes_counter += 1
        
        if genomes_counter > 0:
            printline("Estimating the bloom filter size on {} genomes".format(genomes_counter))

            # Estimate the bloom filter size
            filter_size = estimate_bf_size(os.path.join(tmp_dir, "genomes.txt"), kmer_len, 
                                           os.path.join(tmp_dir, "genomes"), nproc=nproc)

            # Increment the estimated bloom filter size
            increment = int(math.ceil(filter_size*increase_filter_size/100.0))
            filter_size += increment

            with open(os.path.join(db_dir, "manifest.txt"), "a+") as manifest:
                manifest.write("--filter-size {}\n".format(filter_size))
    
    if genomes_paths and filter_size and kmer_len:
        # Retrieve the current working directory
        current_folder = os. getcwd()

        printline("Indexing all the taxonomic levels")

        # Define a partial function around the howdesbt wrapper
        howdesbt_partial = partial(howdesbt, kmer_len=kmer_len, filter_size=filter_size, nproc=nproc)

        # Iterate over all the taxonomic levels from the species up to the kingdom
        for level in ["species", "genus", "family", "order", "class", "phylum", "kingdom"]:
            printline("Running HowDeSBT at the {} level".format(level))
            folders = [str(path) for path in Path(db_dir).glob("**/{}__*".format(level[0]))]
            
            # Initialise the progress bar
            pbar = tqdm.tqdm(total=len(folders), disable=(not verbose))

            with mp.Pool(processes=parallel) as pool:
                def update_bar(*args):
                    # Update the progress bar
                    pbar.update()
                    return
                
                # Process the NCBI tax IDs
                jobs = [pool.apply_async(howdesbt_partial, args=(level_dir, ), callback=update_bar) 
                            for level_dir in folders]
            
            # Close the progress bar
            pbar.close()

        # Also run HowDeSBT on the database folder to build 
        # the bloom filter representation of the kingdom
        printline("Building the database root bloom filter with HowDeSBT")
        howdesbt_partial(db_dir)

        # The howdesbt function automatically set the current working directory to 
        # the index folder of the taxonomic labels
        # Come back to the original folder
        os.chdir(current_folder)

def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether the database folder exists
    if it_exists(os.path.join(args.db_dir, "k__{}".format(args.kingdom)), path_type="folder"):
        raise Exception(("An indexed version of the {} kingdom already exists in the database!\n"
                         "Please use the update module to add new genomes").format(args.kingdom))
    
    # Create the database folder
    os.makedirs(args.db_dir, exist_ok=True)

    # Skip the bloom filter size estimation if the filter size is passed as input
    if args.filter_size:
        args.estimate_filter_size = False
    else:
        if not args.estimate_filter_size:
            raise Exception(("Please specify a bloom filter size with the --filter-size option or "
                             "use the --estimate-filter-size flag to automatically estimate the bloom filter size with ntCard"))

    # Define the database manifest file
    manifest_filepath = os.path.join(db_dir, "manifest.txt")
    if os.path.exists(manifest_filepath):
        # Load and compare kmer-len and filter-size
        with open(manifest_filepath) as manifest:
            for line in manifest:
                line = line.strip()
                if line:
                    line_split = line.split(" ")
                    if line_split[0] == "--kmer-len":
                        if args.kmer_len != int(line_split[1]):
                            raise Exception("The kmer length is not compatible with the specified database")
                    elif line_split[0] == "--filter-size":
                        if args.filter_size:
                            if args.filter_size != int(line_split[1]):
                                raise Exception("The bloom filter size is not compatible with the specified database")
    else:
        with open(manifest_filepath, "w+") as manifest:
            manifest.write("--kmer-len {}\n".format(kmer_len))
            if args.filter_size:
                manifest.write("--filter-size {}\n".format(filter_size))

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    t0 = time.time()

    index(args.db_dir, args.kingdom, args.tmp_dir, args.kmer_len, args.filter_size, limit_genomes=args.limit_genomes,
          max_genomes=args.max_genomes, min_genomes=args.min_genomes, estimate_filter_size=args.estimate_filter_size, 
          increase_filter_size=args.increase_filter_size, completeness=args.completeness, contamination=args.contamination, 
          dereplicate=args.dereplicate, similarity=args.similarity, logger=logger, verbose=args.verbose, nproc=args.nproc, 
          pplacer_threads=args.pplacer_threads, parallel=args.parallel)

    if args.cleanup:
        # Remove the temporary folder
        println("Cleaning up temporary space".format(level), logger=logger, verbose=args.verbose)
        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()
    println("Total elapsed time {}s".format(int(t1 - t0)), 
            logger=logger, verbose=args.verbose)

if __name__ == "__main__":
    main()