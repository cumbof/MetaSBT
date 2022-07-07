#!/usr/bin/env python3
"""
Define cluster-specific boundaries as the minimum and maximum number of common kmers among all the genomes under a specific cluster
"""

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jul 6, 2022"

import sys, os, time, errno, shutil
import argparse as ap
from pathlib import Path
from functools import partial
from logging import Logger

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import get_boundaries, init_logger, it_exists, load_manifest, kmtricks_matrix, number, println
except:
    pass

# Define the module name
TOOL_ID = "boundaries"

# Define the list of dependencies
DEPENDENCIES = ["kmtricks"]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",   # Database folder path
    "--log",      # Path to the log file
    "--output",   # Output table path
    "--tmp-dir"   # Temporary folder path
]

def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(prog=TOOL_ID,
                          description="Define taxonomy-specific boundaries based on kmers for the definition of new clusters",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--cleanup",
                    action = "store_true",
                    default = False,
                    help = "Remove temporary data at the end of the pipeline" )
    p.add_argument( "--db-dir",
                    type = os.path.abspath,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--flat-structure",
                    action = "store_true",
                    default = False,
                    dest = "flat_structure",
                    help = "Genomes in the database have been organized without a taxonomic structure" )
    p.add_argument( "--kingdom",
                    type = str,
                    help = "Consider genomes whose lineage belongs to a specific kingdom" )
    p.add_argument( "--log",
                    type = os.path.abspath,
                    help = "Path to the log file" )
    p.add_argument( "--min-genomes",
                    type = number(int, minv=3),
                    default = 3,
                    dest = "min_genomes",
                    help = "Consider clusters with at least this number of genomes" )
    p.add_argument( "--nproc",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "This argument refers to the number of processors used for parallelizing the pipeline when possible" )
    p.add_argument( "--output",
                    type = os.path.abspath,
                    required = True,
                    help = "Output file with kmer boundaries for each of the taxonomic labels in the database" )
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

def define_boundaries(level_dir: str, level_id: str, tmp_dir: str, output: str, kmer_len: int, filter_size: int, min_genomes: int=3, nproc: int=1) -> None:
    """
    Compute boundaries for the specified taxonomic level

    :param level_dir:       Path to the taxonomic level folder
    :param level_id:        ID of the taxonomic level
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param kmer_len:        Length of the kmers
    :param filter_size:     Size of the bloom filters
    :param min_genomes:     Consider clusters with at least this number of genomes
    :param nproc:           Make the process parallel when possible
    """

    # Search and merge all the reference genomes paths under all references.txt files in the current taxonomic level
    samples = dict()
    references_paths = list(Path(level_dir).glob("**/references.txt"))      # Genomes are usually listed in references.txt files
    references_paths.extend(list(Path(level_dir).glob("**/genomes.txt")))   # Databases without a taxonomic structure use genomes.txt
    
    for references_path in references_paths:
        path_split = str(references_path).split(os.sep)
        next_level = "NA" # In case the current level_id is species
        for path_pos, path_level in enumerate(path_split):
            if path_level.strip():
                if ("{}__".format(path_level[0]) == "{}__".format(level_id[0])) and level_id != "species":
                    next_level = path_split[path_pos+1]
                    break
        
        if level_id != "species":
            samples[next_level] = list()
        with open(str(references_path)) as references:
            for line in references:
                line = line.strip()
                if line:
                    if level_id == "species":
                        samples[line] = [os.path.join(os.path.dirname(str(references_path)), "genomes", "{}.fna.gz".format(line))]
                    else:
                        samples[next_level].append(os.path.join(os.path.dirname(str(references_path)), "genomes", "{}.fna.gz".format(line)))
    
    # In case the current taxonomic level is not the species level
    if level_id != "species":
        # Get rid of clusters with not enough genomes according to min_genomes
        for sample_id in list(samples.keys()):
            if len(samples[sample_id]) < min_genomes:
                del samples[sample_id]
    
    # In case the number of genomes in the current taxonomic level
    # is greater than or equals to the minimum number of genomes specified in input
    if len(samples) >= min_genomes:
        # Create a temporary folder for the specific taxonomic level
        tmp_level_dir = os.path.join(tmp_dir, "boundaries", level_id, os.path.basename(level_dir))
        os.makedirs(tmp_level_dir, exist_ok=True)

        with open(os.path.join(tmp_level_dir, "genomes.fof"), "w+") as genomes_fof:
            for sample_id in samples:
                genomes_fof.write("{} : {}\n".format(sample_id, " ; ".join(samples[sample_id])))

        # Run kmtricks to build the kmers matrix
        kmtricks_matrix(os.path.join(tmp_level_dir, "genomes.fof"), 
                        tmp_level_dir, 
                        kmer_len,
                        filter_size,
                        nproc,
                        os.path.join(tmp_level_dir, "kmers_matrix.txt"))
        
        # Extract boundaries from the kmtricks kmers matrix
        all_kmers, min_kmers, max_kmers = get_boundaries(os.path.join(tmp_level_dir, "kmers_matrix.txt"))

        # Get the full lineage from the level folder path
        lineage = list()
        kingdom_found = False
        for level in level_dir.split(os.sep):
            if level.startswith("k__"):
                kingdom_found = True
            if kingdom_found:
                lineage.append(level)
        lineage = "|".join(lineage)

        # In case of --flat-structure
        if len(lineage.strip()) == 0:
            lineage = level_dir

        # Dump results to the boundaries table
        with open(output, "a+") as table:
            table.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, len(samples), all_kmers, min_kmers, max_kmers, 
                                                              round(min_kmers/all_kmers, 3), round(max_kmers/all_kmers, 3)))

def boundaries(db_dir: str, tmp_dir: str, output: str, flat_structure: bool=False, min_genomes: int=3, kingdom: str=None, 
               logger: Logger=None, verbose: bool=False, nproc: int=1) -> None:
    """
    Define boundaries for each of the taxonomic levels in the database
    Boundaries are defined as the minimum and maximum number of common kmers among all the reference genomes under a specific taxonomic level

    :param db_dir:          Path to the database root folder
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param flat_structure:  Genomes in the database have been organized without a taxonomic structure
    :param min_genomes:     Consider clusters with at least this number of genomes
    :param kingdom:         Retrieve genomes that belong to a specific kingdom
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
        file.write("# --db-dir {}\n".format(db_dir))
        if kingdom:
            file.write("# --kingdom {}\n".format(kingdom))
        file.write("# --min-genomes {}\n".format(min_genomes))
        file.write("# {}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Lineage",       # Taxonomic label
                                                           "References",    # Number of reference genomes or clustrs under a specific taxonomic level
                                                           "Kmers",         # Total number of kmers
                                                           "Min kmers",     # Minimum number of common kmers among the reference genomes/clusters
                                                           "Max kmers",     # Maximum number of common kmers among the reference genomes/clusters
                                                           "Min score",     # Percentage of min kmers on the total number of kmers/clusters
                                                           "Max score"))    # Percentage of max kmers on the total number of kmers/clusters
    
    # Check whether the manifest file exists
    manifest_filepath = os.path.join(db_dir, "manifest.txt")
    if not it_exists(manifest_filepath, path_type="file"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)
    
    # Load the manifest file
    manifest = load_manifest(manifest_filepath)
    if "kmer_len" not in manifest or "filter_size" not in manifest:
        raise Exception("Manifest file does not contain --kmer-len and --filter-size information: {}".format(manifest_filepath))

    # Check whether the genomes folder exists under the database root directory
    if flat_structure:
        # This means that the database has been build with the --flat-structure option
        printline("Defining boundaries")

        # Treat the database as the species level
        define_boundaries(db_dir, "species", tmp_dir, output, manifest["kmer_len"], manifest["filter_size"],
                          min_genomes=min_genomes, nproc=nproc)

    else:
        # Genomes have been taxonomically organized
        target_dir = db_dir if not kingdom else os.path.join(db_dir, "k__{}".format(kingdom))
        levels = ["species", "genus", "family", "order", "class", "phylum"]
        if not kingdom:
            levels.append("kingdom")

        # Iterate over the taxonomic levels
        for level in levels:
            printline("Defining {} boundaries".format(level))
            for level_dir in Path(target_dir).glob("**/{}__*".format(level[0])):
                if os.path.isdir(str(level_dir)):
                    # Define boundaries for the current taxonomic level
                    define_boundaries(str(level_dir), level, tmp_dir, output, manifest["kmer_len"], manifest["filter_size"],
                                      min_genomes=min_genomes, nproc=nproc)
        
        if kingdom:
            # Also define boundaries for the specified kingdom
            define_boundaries(os.path.join(db_dir, kingdom), "kingdom", tmp_dir, output, manifest["kmer_len"], manifest["filter_size"],
                              min_genomes=min_genomes, nproc=nproc)

    # Report the path to the output boundaries table file
    printline("Output table: {}".format(output))

def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether the database folder exists
    target_dir = args.db_dir if not args.kingdom else os.path.join(args.db_dir, "k__{}".format(args.kingdom))
    if not it_exists(target_dir, path_type="folder"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), target_dir)
    
    # Check whether the output boundaries table alrady exists
    if it_exists(args.output, path_type="file"):
        raise Exception("The output boundaries table already exists")

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    t0 = time.time()

    boundaries(args.db_dir, args.tmp_dir, args.output, flat_structure=args.flat_structure, min_genomes=args.min_genomes, 
               kingdom=args.kingdom, logger=logger, verbose=args.verbose, nproc=args.nproc)

    if args.cleanup:
        # Remove the temporary folder
        println("Cleaning up temporary space".format(level), logger=logger, verbose=args.verbose)
        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()
    println("Total elapsed time {}s".format(int(t1 - t0)), 
            logger=logger, verbose=args.verbose)

if __name__ == "__main__":
    main()