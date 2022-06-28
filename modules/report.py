#!/usr/bin/env python3
"""
Create the report table for a specific database
"""

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 28, 2022"

import sys, os, time, errno
import argparse as ap
from pathlib import Path

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import it_exists
except:
    pass

# Define the module name
TOOL_ID = "report"

# Define the list of dependencies
DEPENDENCIES = list()

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",      # Database folder path
    "--output-file"  # Output file path
]

def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """
    
    p = ap.ArgumentParser(prog=TOOL_ID,
                          description="Build the database report table",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--db-dir",
                    type = os.path.abspath,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--output-file",
                    type = os.path.abspath,
                    required = True,
                    dest = "output_file",
                    help = "This is the path to the output table" )
    p.add_argument( "-v",
                    "--version",
                    action = "version",
                    version = "{} version {} ({})".format(TOOL_ID, __version__, __date__),
                    help = "Print the current {} version and exit".format(TOOL_ID) )
    return p.parse_args()

def report(db_dir: str, output_file: str) -> None:
    """
    Build the database report table

    :param db_dir:       Path to the database root folder
    :param output_file:  Path to the output table
    """

    # Initialise the report table
    with open(output_file, "w+") as output:
        # Write header line
        output.write("# {}\t{}\t{}\t{}\t{}\t{}\n".format("Lineage",                     # Taxonomic label
                                                         "MAGs",                        # Number of MAGs
                                                         "Reference genomes",           # Number of reference genomes
                                                         "Mean completeness",           # Mean QC completeness
                                                         "Mean contamination",          # Mean QC contamination
                                                         "Mean strain heterogeneity"))  # Mean QC strain heterogeneity
        
        # Search for all the species folders
        gen = Path(db_dir).glob("**/s__*")
        for species_dir in gen:
            if os.path.isdir(str(species_dir)):
                # Build the lineage from the folder path
                path_split = str(species_dir).split(os.sep)
                levels = list()
                collect = False
                for subpath in path_split:
                    if subpath.startswith("k__"):
                        collect = True
                    if collect:
                        levels.append(subpath)
                lineage = "|".join(levels)

                # Load the list of all reference genomes and MAGs
                mags_list = list()
                references_list = list()

                # In case the mags.txt file exists in the current species folder
                if it_exists(os.path.join(str(species_dir), "mags.txt"), path_type="file"):
                    with open(os.path.join(str(species_dir), "mags.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                mags_list.append(line)
                
                # Do the same with the genomes.txt file
                if it_exists(os.path.join(str(species_dir), "references.txt"), path_type="file"):
                    with open(os.path.join(str(species_dir), "references.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                references_list.append(line)

                # Keep track of all the completeness, contamination, and strain heterogeneity percentages
                # for all the genomes in the current species
                qc_stats = dict()

                # Also load the CheckM table
                if it_exists(os.path.join(str(species_dir), "checkm.tsv"), path_type="file"):
                    with open(os.path.join(str(species_dir), "checkm.tsv")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                line_split = line.split("\t")
                                qc_stats[line_split[0]] = {
                                    "completeness": float(line_split[-3]),
                                    "contamination": float(line_split[-2]),
                                    "strain_heterogeneity": float(line_split[-1])
                                }

                # Take track of the number of MAGs and reference genomes
                # Retrieve also the QC stats
                mags = 0
                references = 0
                completeness = list()
                contamination = list()
                strain_heterogeneity = list()

                # Check whether the reference genomes exist before counting
                for ref in references_list:
                    if os.path.join(os.path.join(str(species_dir), "genomes", "{}.fna.gz".format(ref))):
                        references += 1
                        
                        # Also retrieve its QC stats
                        if ref in qc_stats:
                            completeness.append(qc_stats[ref]["completeness"])
                            contamination.append(qc_stats[ref]["contamination"])
                            strain_heterogeneity.append(qc_stats[ref]["strain_heterogeneity"])

                # Do the same with the MAGs
                for mag in mags_list:
                    if os.path.join(os.path.join(str(species_dir), "genomes", "{}.fna.gz".format(mag))):
                        mags += 1

                        # Also retrieve its QC stats
                        if ref in qc_stats:
                            completeness.append(qc_stats[mag]["completeness"])
                            contamination.append(qc_stats[mag]["contamination"])
                            strain_heterogeneity.append(qc_stats[mag]["strain_heterogeneity"])

                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage,
                                                               mags,
                                                               references,
                                                               round(sum(completeness)/len(completeness), 2) if completeness else 0.0,
                                                               round(sum(contamination)/len(contamination), 2) if contamination else 0.0,
                                                               round(sum(strain_heterogeneity)/len(strain_heterogeneity), 2) if strain_heterogeneity else 0.0))

def main() -> None:
    # Load command line parameters
    args = read_params()

    # Check whether the input database folder path exists
    if not it_exists(args.db_dir, path_type="folder"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    
    # Check whether the output file and folder exist
    output_folder = os.path.dirname(args.output_file)
    if it_exists(args.output_file, path_type="file"):
        raise Exception("The output file already exists")
    if not it_exists(output_folder, path_type="folder"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.output_file)

    report(args.db_dir, args.output_file)

if __name__ == "__main__":
    main()