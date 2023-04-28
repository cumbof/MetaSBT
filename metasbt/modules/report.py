#!/usr/bin/env python3
"""
Create a database report table
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 13, 2023"

import argparse as ap
import datetime
import errno
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import get_bf_density  # type: ignore
except Exception:
    pass

# Define the module name
TOOL_ID = "report"

# Define the list of dependencies
DEPENDENCIES: List[str] = list()

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Database folder path
    "--diff",  # Path to a report table
    "--output-file",  # Output file path
]


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Create a database report table",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required=True,
        dest="db_dir",
        help="This is the database directory with the taxonomically organised sequence bloom trees",
    )
    p.add_argument(
        "--diff",
        type=os.path.abspath,
        dest="diff",
        help="Path to a report table that will be compared with the report produced on --db-dir",
    )
    p.add_argument(
        "--output-file",
        type=os.path.abspath,
        required=True,
        dest="output_file",
        help="This is the path to the output table",
    )
    p.add_argument(
        "--skip-metadata",
        action="store_true",
        default=False,
        dest="skip_metadata",
        help="Do not compare reports metadata in case of --diff",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def report(db_dir: str, output_file: str) -> None:
    """
    Build the database report table

    :param db_dir:       Path to the database root folder
    :param output_file:  Path to the output table
    """

    # Initialise the report table
    with open(output_file, "w+") as output:
        # Write metadata
        output.write("# {} v{} ({})\n".format(TOOL_ID, __version__, __date__))
        output.write("# Database: {}\n".format(db_dir))
        output.write("# Timestamp: {}\n".format(datetime.datetime.today().strftime("%Y%m%d")))

        # Write header line
        output.write(
            "# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                "Cluster",
                "Lineage",
                "MAGs",
                "Reference genomes",
                "Bloom Filter Density",
                "Mean completeness",
                "Mean contamination",
                "Mean strain heterogeneity",
            )
        )

        # Search for all the species folders
        gen = Path(db_dir).glob("**/s__*")

        for species_dir in gen:
            if os.path.isdir(str(species_dir)):
                # Check whether the metadata.tsv file exists
                # It must be there because it contains the cluster ID
                # If it doesn't exist, throw an exception
                metadata_filepath = os.path.join(str(species_dir), "metadata.tsv")

                if not os.path.isfile(metadata_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), metadata_filepath)

                # Retrieve the cluster ID
                cluster = None

                with open(metadata_filepath) as metadata:
                    for line in metadata:
                        line = line.strip()
                        if line:
                            if line.startswith("# Cluster ID:"):
                                cluster = line.split(" ")[-1]
                
                if not cluster:
                    raise Exception("Unable to retrieve the cluster ID in {}".format(metadata_filepath))
                
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
                if os.path.isfile(os.path.join(str(species_dir), "mags.txt")):
                    with open(os.path.join(str(species_dir), "mags.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                if os.path.isfile(os.path.join(str(species_dir), "strains", "filters", "{}.bf".format(line))):
                                    mags_list.append(line)

                # Do the same with the references.txt file
                if os.path.isfile(os.path.join(str(species_dir), "references.txt")):
                    with open(os.path.join(str(species_dir), "references.txt")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                if os.path.isfile(os.path.join(str(species_dir), "strains", "filters", "{}.bf".format(line))):
                                    references_list.append(line)
                
                if not mags_list and not references_list:
                    raise Exception("Unable to retrieve genomes in {}".format(str(species_dir)))

                # Keep track of all the completeness, contamination, and strain heterogeneity percentages
                # for all the genomes in the current species
                qc_stats = dict()

                # Also load the CheckM table
                if os.path.isfile(os.path.join(str(species_dir), "checkm.tsv")):
                    with open(os.path.join(str(species_dir), "checkm.tsv")) as file:
                        for line in file:
                            line = line.strip()
                            if line:
                                line_split = line.split("\t")
                                qc_stats[line_split[0]] = {
                                    "completeness": float(line_split[-3]),
                                    "contamination": float(line_split[-2]),
                                    "strain_heterogeneity": float(line_split[-1]),
                                }

                # Retrieve the QC stats
                completeness = list()
                contamination = list()
                strain_heterogeneity = list()

                for genome in references_list + mags_list:
                    if genome in qc_stats:
                        completeness.append(qc_stats[genome]["completeness"])
                        contamination.append(qc_stats[genome]["contamination"])
                        strain_heterogeneity.append(qc_stats[genome]["strain_heterogeneity"])

                output.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        cluster,
                        lineage,
                        len(mags_list),
                        len(references_list),
                        get_bf_density(os.path.join(str(species_dir), "{}.bf".format(os.path.basename(str(species_dir))))),
                        round(sum(completeness) / len(completeness), 2) if completeness else "na",
                        round(sum(contamination) / len(contamination), 2) if contamination else "na",
                        round(sum(strain_heterogeneity) / len(strain_heterogeneity), 2) if strain_heterogeneity else "na",
                    )
                )


def load_report_table(table_path: str) -> Tuple[str, str, datetime.datetime, Dict[str, Dict[str, Union[str, int, float]]]]:
    """
    Load a report table into a dictionary

    :param table_path:  Path to the report table
    :return:            Dictionary with the content of the report table indexed by cluster ID
    """

    database = None
    version = None
    timestamp = None
    report = dict()

    with open(table_path) as table:
        header = list()

        for line in table:
            line = line.strip()
            if line:

                if line.startswith("#"):
                    if line[1:].strip().startswith(TOOL_ID):
                        # Retrieve the version of this script used to build the report table
                        version = " ".join(line[1:].strip().split(" ")[1:])

                    elif line[1:].strip().startswith("Database"):
                        # Retrieve the database
                        database = line.split(":")[1].strip()

                    elif line[1:].strip().startswith("Timestamp"):
                        # Retrieve the time when the report table has been produced
                        timestamp = datetime.datetime.strptime(line.split(":")[1].strip(), "%Y-%m-%d %H:%M:%S.%f")

                    else:
                        # This is the header line with table columns
                        header = line[1:].strip().split("\t")

                else:
                    line_split = line.split("\t")

                    # Read the content of the report table and load it in a dictionary indexed by cluster ID
                    report[line_split[0]] = {
                        "lineage": line_split[1],
                        "mags": int(line_split[2]),
                        "references": int(line_split[3]),
                        "bf_density": float(line_split[4]),
                        "mean_completeness": float(line_split[5]),
                        "mean_contamination": float(line_split[6]),
                        "mean_strain_heterogeneity": float(line_split[7])
                    }

    return database, version, timestamp, report


def diff(new_table_path: str, old_table_path: str, skip_metadata: bool = False) -> None:
    """
    Compare two report tables and report differences

    :param new_table_path:  Path to the most recent report table
    :param old_table_path:  Path to an old report table
    """

    # Load new and old tables
    database_t1, version_t1, timestamp_t1, table_t1 = load_report_table(new_table_path)
    database_t2, version_t2, timestamp_t2, table_t2 = load_report_table(old_table_path)

    # Compare tables metadata
    if not skip_metadata:
        # Compare tool versions
        if version_t1 != version_t2:
            # This is just a warning
            print("Warning: the input report tables have been produced with two different versions of the \"{}\" tool".format(TOOL_ID))

        # Compare databases
        if database_t1 != database_t2:
            raise Exception(
                ("Unable to compare the input report tables!\n"
                 "You are trying to compare two reports built over two different databases")
            )

        # Compare timestamps
        # timestamp_t1 and timestamp_t2 are datetime.datetime objects
        if timestamp_t1 <= timestamp_t2:
            raise Exception(
                ("Unable to compare the input report tables!\n"
                 "The timestamp in the table you specified under --diff seems more recent than the newer table's timestamp")
            )

    # Compare lineages
    # Count how many references and MAGs are in the report table 1
    table_t1_references = 0
    table_t1_mags = 0
    table_t1_uncharacterized = 0

    print("Data in {}".format(os.path.basename(new_table_path)))

    for cluster in table_t1:
        table_t1_references += table_t1[cluster]["references"]
        table_t1_mags += table_t1[cluster]["mags"]

        if table_t1[cluster]["references"] == 0:
            table_t1_uncharacterized += 1

    print("\tReferences: {}".format(table_t1_references))
    print("\tMAGs: {}".format(table_t1_mags))
    print("\tUncharacterized: {}".format(table_t1_uncharacterized))

    # Do the same with the report table 2
    table_t2_references = 0
    table_t2_mags = 0
    table_t2_uncharacterized = 0

    print("Data in {}".format(os.path.basename(old_table_path)))

    for cluster in table_t2:
        table_t2_references += table_t2[cluster]["references"]
        table_t2_mags += table_t2[cluster]["mags"]

        if table_t2[cluster]["references"] == 0:
            table_t2_uncharacterized += 1

    print("\tReferences: {}".format(table_t2_references))
    print("\tMAGs: {}".format(table_t2_mags))
    print("\tUncharacterized: {}".format(table_t2_uncharacterized))

    # Intersect clusters
    intersection = set(table_t1.keys()).intersection(set(table_t2.keys()))

    # Report suppressed clusters
    suppressed = set(table_t2.keys()).difference(intersection)
    print("Suppressed clusters ({}):".format(len(suppressed)))

    for cluster in suppressed:
        print("\t{}".format(cluster))

    # Report new clusters
    new = set(table_t1.keys()).difference(intersection)
    print("New clusters ({}):".format(len(new)))

    for cluster in new:
        print("\t{}".format(cluster))

    # Count how many unknown clusters have been characterized
    characterized = list()

    for cluster in intersection:
        if table_t1[cluster]["references"] > 0 and table_t2[cluster]["references"] == 0:
            characterized.append(cluster)

    print("Characterized clusters ({}):".format(len(characterized)))

    for cluster in characterized:
        print("\t{}".format(cluster))


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Check whether the input database folder path exists
    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    # Check whether the output file already exists
    if os.path.isfile(args.output_file):
        raise Exception("The output file already exists")

    # Check whether the output folder exists
    output_folder = os.path.abspath(os.path.dirname(args.output_file))

    if not os.path.isdir(output_folder):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_folder)

    report(args.db_dir, args.output_file)

    # Check whether there is a report table that must be
    # compared to the output table
    if args.diff:
        # Check whether the report table file exists
        if not os.path.isfile(args.diff):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.diff)

        # Check for differences with the new report table
        diff(args.output_file, args.diff, skip_metadata=args.skip_metadata)


if __name__ == "__main__":
    main()
