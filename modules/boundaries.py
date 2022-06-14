#!/usr/bin/env python

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 12, 2022"

import sys, os, time, errno, shutil
import argparse as ap
from pathlib import Path
from itertools import partial
from utils import init_logger, it_exists, println

# Define the module name
TOOL_ID = "boundaries"

def read_params():
    p = ap.ArgumentParser(description="Define taxonomy-specific boundaries based on kmers for the definition of new clusters",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--cleanup",
                    action = "store_true",
                    default = False,
                    help = "Remove temporary data at the end of the pipeline" )
    p.add_argument( "--db-dir",
                    type = str,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--kingdom",
                    type = str,
                    required = True,
                    choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
                    help = "Consider genomes whose lineage belongs to a specific kingdom" )
    p.add_argument( "--log",
                    type = str,
                    help = "Path to the log file" )
    p.add_argument( "--min-genomes",
                    type = number(int, minv=3),
                    default = 3,
                    dest = "min_genomes"
                    help = "Consider clusters with at least this number of genomes" )
    p.add_argument( "--nproc",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "This argument refers to the number of processors used for parallelizing the pipeline when possible" )
    p.add_argument( "--output",
                    type = str,
                    required = True,
                    help = "Output file with kmer boundaries for each of the taxonomic labels in the database" )
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

def define_boundaries(level_dir, level_id, tmp_dir, output, min_genomes=3, nproc=1):
    """
    Compute boundaries for the specified taxonomic level

    :param level_dir:       Path to the taxonomic level folder
    :param level_id:        ID of the taxonomic level
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param min_genomes:     Consider clusters with at least this number of genomes
    :param nproc:           Make the process parallel when possible
    """

    # Create a temporary folder for the specific taxonomic level
    tmp_level_dir = os.path.join(tmp_dir, "boundaries", level_id, os.path.basename(level_dir))
    makedirs(tmp_level_dir, exist_ok=True)

    # Search and merge all the reference genomes paths under all references.txt files in the current taxonomic level
    how_many = 0
    for references_path in Path(level_dir).glob("**/references.txt"):
        with open(os.path.join(tmp_level_dir, "genomes.fof"), "a+") as genomes_fof:
            with open(str(references_path)) as references:
                for line in references:
                    line = line.strip()
                    if line:
                        genomes_fof.write("{} : {}\n".format(line, os.path.join(os.path.dirname(str(references_path)), "genomes", "{}.fna.gz".format(line))))
                        how_many += 1
    
    # In case the number of genomes in the current taxonomic level
    # is greater than or equals to the minimum number of genomes specified in input
    if how_may >= min_genomes:
        # Run kmtricks to build the kmers matrix
        kmtricks_matrix(os.path.join(tmp_level_dir, "genomes.fof"), 
                        os.path.join(tmp_level_dir, "matrix"), 
                        nproc,
                        os.path.join(tmp_level_dir, "kmers_matrix.txt")):
        
        # Extract boundaries from the kmtricks kmers matrix
        min_kmers, max_kmers = get_boundaries(os.path.join(tmp_level_dir, "kmers_matrix.txt"))

        # Get the full lineage from the level folder path
        lineage = list()
        kingdom_found = False
        for level in level_dir.split(os.sep):
            if level.startswith("k__"):
                kingdom_found = True
            if kingdom_found:
                lineage.append(level)
        lineage = "|".join(lineage)

        # Dump results to the boundaries table
        with open(output, "a+") as table:
            table.write("{}\t{}\t{}\n".format(lineage, min_kmers, max_kmers))

def boundaries(db_dir, kingdom, tmp_dir, output, min_genomes=3, logger=None, verbose=False, nproc=1):
    """
    Define boundaries for each of the taxonomic levels in the database
    Boundaries are defined as the minimum and maximum number of common kmers among all 
    the reference genomes under a specific taxonomic level

    :param db_dir:          Path to the database root folder
    :param kingdom:         Retrieve genomes that belong to a specific kingdom
    :param tmp_dir:         Path to the temporary folder
    :param output:          Path to the output table file with boundaries
    :param min_genomes:     Consider clusters with at least this number of genomes
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
        file.write("# --kingdom {}\n".format(kingdom))
        file.write("# --min-genomes {}\n".format(min_genomes))
        file.write("# {}\t{}\t{}\n".format("Lineage",
                                           "Min kmers",
                                           "Max kmers"))
    
    # Iterate over the taxonomic levels from species up to the phylum
    for level in ["species", "genus", "family", "order", "class", "phylum"]:
        printline("Defining {} boundaries".format(level))
        for level_dir in Path(os.path.join(db_dir, kingdom)).glob("**/{}__*".format(level[0])):
            if os.path.isdir(str(level_dir)):
                # Define boundaries for the current taxonomic level
                define_boundaries(str(level_dir), level, tmp_dir, output, 
                                  min_genomes=min_genomes, nproc=nproc)
    
    # Also define boundaries for the specified kingdom
    define_boundaries(os.path.join(db_dir, kingdom), "kingdom", tmp_dir, output, 
                      min_genomes=min_genomes, nproc=nproc)

    # Report the path to the output boundaries table file
    printline("Output table: {}".format(output))

def main():
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, verbose=args.verbose)

    # Check whether the database folder exists
    if not it_exists(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    
    # Check whether the output boundaries table alrady exists
    if os.path.exists(args.output):
        raise Exception("The output boundaries table already exists")

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    t0 = time.time()

    boundaries(args.db_dir, args.kingdom, args.tmp_dir, args.output, min_genomes=args.min_genomes, 
               logger=logger, verbose=args.verbose, nproc=args.nproc)

    if args.cleanup:
        # Remove the temporary folder
        println("Cleaning up temporary space".format(level), logger=logger, verbose=verbose)
        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()
    println("Total elapsed time {}s".format(int(t1 - t0)), 
            logger=logger, verbose=args.verbose)

if __name__ == "__main__":
    main()