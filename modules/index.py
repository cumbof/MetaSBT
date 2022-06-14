#!/usr/bin/env python

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 10, 2022"

import sys, os, time, tarfile, gzip, re, shutil
import argparse as ap
from pathlib import Path
from itertools import partial
from utils import checkm, download, filter_genomes, howdesbt, init_logger, it_exists, kmtricks_matrix, number, println, run

# Define the module name
TOOL_ID = "index"

# Define the url to the NCBI taxdump
TAXDUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

def read_params():
    p = ap.ArgumentParser(description=("Retrieve complete reference genomes from the NCBI GenBank and build a "
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
                    type = str,
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
    p.add_argument( "--how-many",
                    type = number(int, minv=1),
                    dest = "how_many",
                    help = ("Limit the number of genomes per species. "
                            "The number of genomes per species is not limited by default") )
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
    p.add_argument( "--log",
                    type = str,
                    help = "Path to the log file" )
    p.add_argument( "--nproc",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "This argument refers to the number of processors used for parallelizing the pipeline when possible" )
    p.add_argument( "--pplacer_threads",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "Maximum number of threads for pplacer. This is required to maximise the CheckM performances" )
    p.add_argument( "--similarity",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 100.0,
                    help = ("Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. "
                            "This is used exclusively in conjunction with the --dereplicate argument") )
    p.add_argument( "--tmp-dir",
                    type = str,
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

def level_name(current_level, prev_level):
    """
    Define a taxonomic level name

    :param current_level:   Current level name
    :param prev_level:      Previous level name in case of unclassified
    """

    # Remove special characters from current and previous level names
    regex = "[^A-Za-z0-9\-]+"
    current_level = re.sub(regex, "", current_level)
    prev_level = re.sub(regex, "", prev_level)

    # Build the new level name
    level_prefix = current_level.strip()
    level_suffix = ""
    if not level_prefix:
        level_prefix = prev_level
        # Fill empty taxa levels with unclassified
        level_suffix = "unclassified"
    
    return "{}_{}".format(level_prefix, level_suffix)

def index(db_dir, kingdom, tmp_dir, kmer_len, filter_size, how_many=None,
          estimate_filter_size=False, increase_filter_size=0.0, completeness=0.0, contamination=100.0,
          dereplicate=False, similarity=100.0, logger=None, verbose=False, nproc=1, pplacer_threads=1):
    """
    Build the database baseline

    :param db_dir:                  Path to the database root folder
    :param kingdom:                 Retrieve genomes that belong to a specific kingdom
    :param tmp_dir:                 Path to the temporary folder
    :param kmer_len:                Length of the kmers
    :param filter_size:             Size of the bloom filters
    :param how_many:                Consider species with at most this number of genomes
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
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Download the NCBI taxdump
    printline("Downloading taxonomy dump from NCBI")
    taxdump = download(TAXDUMP_URL, tmp_dir)
    
    # Raise an exception in case the taxdump.tar.gz file does not exist
    if not os.path.exists(taxdump):
        raise Exception("Unable to retrieve data from remote location\n{}".format(TAXDUMP_URL))

    # Create the taxdump folder in the temporary directory
    taxdump_dir = os.path.join(tmp_dir, "taxdump")
    os.makedirs(taxdump_dir, exist_ok=True)

    # Decompress the archive
    with tarfile.open(source, "r:gz") as tar:
        tar.extractall(taxdump_dir)

    # Run ncbitax2lin to extract lineages
    printline("Exporting lineages")
    ncbitax2lin_table = os.path.join(tmp_dir, "ncbi_lineages.csv.gz")
    run(["ncbitax2lin", "--nodes-file", os.path.join(taxdump_dir, "nodes.dmp"),
                        "--names-file", os.path.join(taxdump_dir, "names.dmp"),
                        "--output", ncbitax2lin_table],
        silence=True)

    # Raise an exception in case the ncbi_lineages.csv.gz file does not exist
    if not os.path.exists(ncbitax2lin_table):
        raise Exception("An error has occurred while running ncbitax2lin")

    # Build a mapping between tax IDs and full taxonomic labels
    # Take track of the taxonomic labels to remove duplicates
    # Taxonomy IDs are already sorted in ascending order in the ncbitax2lin output table
    printline("Building tax ID to full taxonomy mapping")
    taxonomies = list()
    tax_table_filepath = os.path.join(tmp_dir, "taxa.tsv")
    with open(tax_table_filepath, "w+") as tax_table:
        # Add the header line
        tax_table.write("# {}\t{}\n".format("tax_id", "tax_label"))
        with gzip.open(ncbitax2lin_table, "rt") as ncbi_table:
            line_count = 0
            for line in ncbi_table:
                line = line.strip()
                if line:
                    # Exclude the header line
                    if line_count > 0:
                        line_split = line.split(",")
                        # Consider a specific kingdom only
                        if line_split[1] == kingdom:
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
                            # Exclude already processed taxonomic labels
                            if label not in taxonomies:
                                # Exclude unclassified genomes
                                if not "unclassified" in label:
                                    tax_label.write("{}\t{}\n".format(line_split[0], label))
                                # Update the list of processed taxonomic labels
                                taxonomies.append(label)
                    line_count += 1
    
    # Create a temporary folder to store the downloaded genomes
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes")
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # Take track of all the genomes paths for estimating the bloom filter size
    if estimate_filter_size:
        genome_paths = list()
    
    printline("Downloading genomes from NCBI GenBank")
    if completeness > 0.0 or contamination < 100.0:
        printline("Running CheckM to quality control input genomes")
        printline("Minimum completeness: {}".format(completeness))
        printline("Maximum contamination: {}".format(contamination))
    if dereplicate:
        printline("Dereplicating input genomes")
        printline("Similarity threshold: {}".format(similarity))
    
    # Iterate over the taxa.tsv file and download genomes from NCBI
    with open(tax_table_filepath) as tax_table:
        for line in tax_table:
            line = line.strip()
            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")
                    
                    # Define the ncbi-genome-download command line
                    cmd_line = ["ncbi-genome-download", "--assembly-levels", "complete",
                                                        "--section", "genbank",
                                                        "--species-taxids", line_split[0], kingdom.lower()]

                    # First check how many genomes belong to the specific tax ID
                    # Define the log file path with the list of genomes retrieved by ncbi-genome-download
                    genomes_list = open(os.path.join(tmp_genomes_dir, "genomes.txt"), "w+")
                    # Append the --dry-run argument to the command line
                    cmd_line.append("--dry-run")
                    # Run ncbi-genome-download
                    run(cmd_line, stdout=genomes_list, stderr=genomes_list)
                    # Close the handler
                    genomes_list.close()
                    # Remove the --dry-run argument from the command line
                    cmd_line.remove("--dry-run")
                    
                    # Count how many genomes belongs to the current tax ID
                    genomes_counter = 0
                    with open(os.path.join(tmp_genomes_dir, "genomes.txt")) as genomes_list:
                        for l in genomes_list:
                            l = l.strip()
                            if l:
                                if l.startswith("GCA_"):
                                    genomes_counter += 1
                    
                    # Check whether the number of genomes is greater than 
                    # the number provided with --how_many
                    retrieve_genomes = True
                    if how_many:
                        if genomes_counter < how_many:
                            retrieve_genomes = False

                    if retrieve_genomes:
                        # Add arguments to the command line
                        cmd_line.extend([
                            "--output-folder", tmp_genomes_dir, "--flat-output",                # Add the --output-folder
                            "--metadata-table", os.path.join(tmp_genomes_dir, "metadata.tsv"),  # Also retrieve genomes metadata
                            "--parallel", str(nproc),                                           # Make it parallel
                            "--retries", "1"                                                    # Retry downloading genomes
                        ])
                        # Run ncbi-genome-download to retrieve all the complete genomes 
                        # related to the current tax ID from NCBI GenBank
                        run(cmd_line, silence=True)
                        
                        # Define the list of paths to the genome files
                        genomes = list()
                        for genome_path in Path(tmp_genomes_dir).glob("*.fna.gz"):
                            # Retrieve the genome name
                            genome_name = os.path.splitext(os.path.splitext(os.path.basename(str(genome_path)))[0])[0].split(".")[0]
                            # Rename the genome file
                            shutil.move(str(genome_path), os.path.join(tmp_genomes_dir, "{}.fna.gz".format(genome_name)))
                            # Add the genome path to the list of genome paths
                            genomes.append(os.path.join(tmp_genomes_dir, "{}.fna.gz".format(genome_name)))

                        # Take track of the paths to the CheckM output tables
                        checkm_tables = []

                        # Quality control
                        if completeness > 0.0 or contamination < 100.0:
                            # Define the CheckM temporary folder
                            checkm_tmp_dir = os.path.join(tmp_dir, "checkm", line_split[0])
                            os.makedirs(checkm_tmp_dir, exist_ok=True)
                            # Run CheckM on the current set of genomes
                            checkm_tables = checkm(genomes, checkm_tmp_dir, file_extension="fna.gz", 
                                                                            nproc=nproc, 
                                                                            pplacer_threads=pplacer_threads)
                            # Filter genomes according to the input --completeness and --contamination thresholds
                            genome_ids = filter_checkm_tables(checkm_tables, completeness=0.0, contamination=100.0)
                            genomes = [os.path.join(tmp_genomes_dir, "{}.fna.gz".format(genome_id)) for genome_id in genome_ids]

                        # Dereplication
                        if genomes and dereplicate:
                            # Define the kmtricks temporary folder
                            kmtricks_tmp_dir = os.path.join(tmp_dir, "kmtricks", line_split[0])
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
                                file1.write("# kmer\t{}\n".format("\t".join(ordered_genome_names)))
                                with open(output_table) as file2:
                                    for line in file2:
                                        line = line.strip()
                                        if line:
                                            file1.write("{}\n".format(line))
                            
                            # Get rid of the matrix withour header line
                            os.unlink(output_table)
                            shutil.move(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), output_table)

                            # Filter genomes according to their percentage of common kmers defined with --similarity
                            dereplicated_genomes_filepath = os.path.join(tmp_dir, "kmtricks", line_split[0], "dereplicated.txt")
                            filter_genomes(output_table, dereplicated_genomes_filepath, similarity=similarity)
                            
                            # Iterate over the list of dereplicated genome and rebuild the list of genome file paths
                            genomes = list()
                            with open(dereplicated_genomes_filepath) as dereplicated_genomes:
                                for l in dereplicated_genomes:
                                    l = l.strip()
                                    if l:
                                        genomes.append(os.path.join(tmp_genomes_dir, "{}.fna.gz".format(l)))
                        
                        # In case at least one genome survived both the 
                        # quality control and dereplication steps
                        if genomes:
                            # Define the taxonomy folder in database
                            tax_dir = os.path.join(db_dir, line_split[1].replace("|", os.sep))
                            genomes_dir = os.path.join(tax_dir, "genomes")
                            os.makedirs(genomes_dir, exist_ok=True)

                            with open(os.path.join(tax_dir, "references.txt"), "w+") as file:
                                # Move the processed genomes to the taxonomy folder
                                for genome_path in genomes:
                                    shutil.move(genome_path, genomes_dir)
                                    # Also update the global list of genome paths for estimating the bloom filter size
                                    if estimate_filter_size:
                                        genome_paths.append(os.path.join(genomes_dir, os.path.basename(genome_path)))
                                    # Also take track of the genome names in the references.txt file
                                    genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                                    if genome_path.endswith(".gz"):
                                        genome_name = os.path.splitext(genome_name)[0]
                                    file.write("{}\n".format(genome_name))
                            
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

    # Check whether the bloom filter size must be estimated
    if genome_paths and estimate_filter_size and not filter_size:
        # Dump the global list of genome path to file
        count_genomes = 0
        with open(os.path.join(tmp_dir, "genomes.txt"), "w+") as file:
            for path in genome_paths:
                if os.path.exists(path):
                    file.write("{}\n".format(path))
                    count_genomes += 1
        
        if count_genomes > 0:
            printline("Running ntCard for estimating the bloom filter size")

            # Estimate the bloom filter size with ntCard
            run(["ntcard", "--kmer={}".format(kmer_len),
                           "--threads={}".format(nproc),
                           "--prefix={}".format(os.path.join(tmp_dir, "genomes")),
                           "@{}".format(os.path.join(tmp_dir, "genomes.txt"))],
                silence=True)

            # Read the ntcard output
            F0 = 0
            f1 = 0
            with open(os.path.join(tmp_dir, "genomes_k{}.hist".format(kmer_len))) as file:
                F0_found = False
                f1_found = False
                for line in file:
                    line = line.strip()
                    if line:
                        if line.startswith("F0\t"):
                            F0 = int(line.split()[-1])
                            F0_found = True
                        elif line.startswith("1\t"):
                            f1 = int(line.split()[-1])
                            f1_found = True
                        if F0_found and f1_found:
                            break
            
            # Compute the bloom filter size
            filter_size = F0-f1

            # Increment the estimated bloom filter size
            increment = int(math.ceil(filter_size*increase_filter_size/100.0))
            filter_size += increment
            
    # Define the database manifest file
    with open(os.path.join(db_dir, "manifest.txt"), "w+") as manifest:
        # Take track of the filter size, kingdom, and kmer length
        manifest.write("--filter-size {}\n".format(filter_size))
        manifest.write("--kingdom {}\n".format(kingdom))
        manifest.write("--kmer-len {}\n".format(kmer_len))
    
    # Retrieve the current working directory
    current_folder = os. getcwd()

    # Iterate over all the taxonomic levels from the species up to the kingdom
    for level in ["species", "genus", "family", "order", "class", "phylum", "kingdom"]:
        for level_dir in Path(db_dir).glob("**/{}__*".format(level[0])):
            if os.path.isdir(str(level_dir)):
                printline("Running HowDeSBT at the {} level".format(level))
                # Index the current taxonomic level with HowDeSBT
                howdesbt(level_dir, kmer_len, filter_size, nproc)

    # Also run HowDeSBT on the database folder to build 
    # the bloom filter representation of the kingdom
    printline("Building the database root bloom filter with HowDeSBT")
    howdesbt(db_dir, kmer_len, filter_size, nproc)

    # The howdesbt function automatically set the current working directory to 
    # the index folder of the taxonomic labels
    # Come back to the original folder
    os.chdir(current_folder)

def main():
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, verbose=args.verbose)

    # Check whether the database folder exists
    if it_exists(args.db_dir):
        raise Exception("The database folder already exists")
    
    # Create the database folder
    os.makedirs(args.db_dir)

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    t0 = time.time()

    index(args.db_dir, args.kingdom, args.tmp_dir, args.kmer_len, args.filter_size, how_many=args.how_many,
          estimate_filter_size=estimate_filter_size, increase_filter_size=increase_filter_size, 
          completeness=completeness, contamination=contamination, dereplicate=dereplicate, 
          similarity=similarity, logger=logger, verbose=args.verbose, nproc=args.nproc, 
          pplacer_threads=args.pplacer_threads)

    if args.cleanup:
        # Remove the temporary folder
        println("Cleaning up temporary space".format(level), logger=logger, verbose=verbose)
        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()
    println("Total elapsed time {}s".format(int(t1 - t0)), 
            logger=logger, verbose=args.verbose)

if __name__ == "__main__":
    main()