#!/usr/bin/env python

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 9, 2022"

import sys, os, time, errno
import argparse as ap
from pathlib import Path

# Define the module name
TOOL_ID = "report"

def read_params():
    p = ap.ArgumentParser(description="Build the database report table",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--db-dir",
                    type = str,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--output-file",
                    type = str,
                    required = True,
                    dest = "output_file",
                    help = "This is the path to the output table" )
    p.add_argument( "-v",
                    "--version",
                    action = "version",
                    version = "{} version {} ({})".format(TOOL_ID, __version__, __date__),
                    help = "Print the current {} version and exit".format(TOOL_ID) )
    return p.parse_args()

def report(db_dir, output_file):
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
                if os.path.exists(os.path.join(str(species_dir), "mags.txt")):
                    with open(os.path.join(str(species_dir), "mags.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                mags_list.append(line)
                
                # Do the same with the genomes.txt file
                if os.path.exists(os.path.join(str(species_dir), "genomes.txt")):
                    with open(os.path.join(str(species_dir), "genomes.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                references_list.append(line)

                # Keep track of all the completeness, contamination, and strain heterogeneity percentages
                # for all the genomes in the current species
                qc_stats = dict()

                # Also load the CheckM table
                if os.path.exists(os.path.join(str(species_dir), "checkm.tsv")):
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
                mags = 0
                references = 0
                completeness = list()
                contamination = list()
                strain_heterogeneity = list()

                # In case the genomes.fof file exists in the current species folder
                if os.path.exists(os.path.join(str(species_dir), "genomes.fof")):
                    with open(os.path.join(str(species_dir), "genomes.fof")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                line_split = line.split(" : ")
                                found = False
                                if line_split[0] in mags_list:
                                    mags += 1
                                    found = True
                                elif line_split[0] in references_list:
                                    references += 1
                                    found = True
                                if found and line_split[0] in qc_stats:
                                    completeness.append(qc_stats[line_split[0]["completeness"]])
                                    contamination.append(qc_stats[line_split[0]["contamination"]])
                                    strain_heterogeneity.append(qc_stats[line_split[0]["strain_heterogeneity"]])

                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage,
                                                               mags,
                                                               references,
                                                               round(sum(completeness)/len(completeness), 2),
                                                               round(sum(contamination)/len(contamination), 2),
                                                               round(sum(strain_heterogeneity)/len(strain_heterogeneity), 2)))

def main():
    # Load command line parameters
    args = read_params()

    # Check whether the input database folder path exists
    if not isinstance(args.db_dir, str):
        raise TypeError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    if not os.path.exists(args.db_dir) or not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    
    # In case the output file path is not specified
    if not args.output_file:
        raise TypeError(errno.ENOENT, os.strerror(errno.ENOENT), args.output_file)
    # Check whether its folder exists
    output_folder = os.path.dirname(args.output_file)
    if not os.path.exists(output_folder):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_folder)

    report(args.db_dir, args.output_file)

if __name__ == "__main__":
    main()