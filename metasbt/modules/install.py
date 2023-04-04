#!/usr/bin/env python3
"""
Install pre-computed MetaSBT databases locally
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 1, 2023"

import argparse as ap
import errno
import math
import os
import shutil
import tarfile
import tempfile
import sys

from tabulate import tabulate

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import download, validate_url, run  # type: ignore  # isort: skip
except Exception:
    pass

# Define the module name
TOOL_ID = "install"

# Define the list of dependencies
DEPENDENCIES = [
    "wget",
]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--install-in",  # Install a database in a specific location
]


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Install pre-computed MetaSBT databases locally",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--database",
        type=str,
        help="Select a database",
    )
    p.add_argument(
        "--database-file",
        type=os.path.abspath,
        dest="database_file",
        help="Path to a local database tarball",
    )
    p.add_argument(
        "--database-url",
        type=str,
        dest="database_url",
        help="URL to a database tarball",
    )
    p.add_argument(
        "--database-version",
        type=str,
        dest="database_version",
        help="Select a database version. Must be used in conjunction with --database. Use the most updated one by default",
    )
    p.add_argument(
        "--hub",
        type=str,
        default="https://raw.githubusercontent.com/cumbof/MetaSBT-DBs/main/databases.tsv",
        help="Path to a local or remote table with a list of MetaSBT databases",
    )
    p.add_argument(
        "--install-in",
        type=os.path.abspath,
        dest="install_in",
        help="Install a database in a specific location",
    )
    p.add_argument(
        "--list-databases",
        action="store_true",
        default=False,
        dest="list_databases",
        help="List all the available MetaSBT databases and exit",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def convert_size(size: float, in_unit: str, out_unit: str) -> float:
    """
    Convert size from a unit to another (e.g., B to GB and viceversa)

    :param size:        Input size
    :param in_unit:     Unit of input size
    :param out_unit:    Convert the input size to this unit
    :return:            Size converted to a specific unit
    """

    # Define a list of units
    units = ["B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"]

    if in_unit.upper() not in units or out_unit.upper() not in units:
        raise ValueError("Cannot convert units")

    if size == 0.0:
        return 0.0

    if in_unit.upper() == out_unit.upper():
        return size

    p = math.pow(1024, abs(units.index(in_unit) - units.index(out_unit)))

    if units.index(in_unit) < units.index(out_unit):
        return size / p

    else:
        return size * p


def fix_paths(filepath: str, main_folder: str, replace_with: str) -> None:
    """
    Fix file paths

    :param filepath:        Path to the file with the list of paths
    :param main_folder:     Fix paths up to this folder
    :param replace_with:    Replace sub-paths with this
    """

    filepaths = list()

    with open(filepath) as infile:
        for line in infile:
            line = line.strip()
            if line:
                # In case of sbt files, count the number of stars before the file path
                stars = line.count("*")

                line_split = line.split(os.sep)
                up_to_idx = line_split.index(main_folder)
                
                # Define the new file path
                # Also prepend the stars in case of sbt files
                filepaths.append(
                    "{}{}".format(
                        "*"*stars,
                        os.path.join(replace_with, os.sep.join(line_split[up_to_idx:]))
                    )
                )

    # This overrides the original file
    with open(filepath, "w+") as infile:
        for line in filepaths:
            infile.write("{}\n".format(line))


def main() -> None:
    # Load command line parameters
    args = read_params()

    if not args.database_file and not args.database_url:
        databases = dict()

        with tempfile.TemporaryDirectory() as tmpdir:
            if os.path.isfile(args.hub):
                # Check whether the hub with the list of databases exists locally
                db_table_filepath = args.hub

            elif validate_url(args.hub):
                # Retrieve the list of pre-computed MetaSBT databases
                db_table_filepath = download(args.hub, tmpdir)

            else:
                raise Exception("Please enter a valid file path or URL with --hub")

            with open(db_table_filepath) as db_table:
                header = list()

                for line in db_table:
                    line = line.strip()
                    if line:
                        if line.startswith("#"):
                            header = line[1:].strip().split("\t")

                        else:
                            line_split = line.split("\t")

                            # Read databases information
                            if line_split[header.index("id")] not in databases:
                                databases[line_split[header.index("id")]] = list()

                            # id, version, references, mags, size, tarball, sha256, and info
                            db_version = {
                                value: line_split[pos] for pos, value in enumerate(header)
                            }

                            databases[line_split[header.index("id")]].append(db_version)

        if not databases:
            raise Exception("No databases available!")

        if args.list_databases:
            # Print the list of databases and exit
            table = [["ID", "Version", "References", "MAGs", "Size", "Info"]]

            for db in databases:
                for db_version in sorted(databases[db], key=lambda v: int(v["version"])):
                    table.append(
                        [
                            db,
                            db_version["version"],
                            db_version["references"],
                            db_version["mags"],
                            db_version["size"],
                            db_version["info"]
                        ]
                    )

            print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))

            sys.exit(os.EX_OK)

        if not args.database or args.database not in databases:
            raise Exception("No database available!")

        if not args.database_version:
            print("Warning: using the most updated version of {}".format(args.database))

            args.database_version = 0

            for db in databases[args.database]:
                # Search for the most updated version
                if db["version"] > args.database_version:
                    args.database_version = db["version"]

                    # Keep track of the link to the tarball
                    tarball_url = db["tarball"]

                    # Also keep track of the tarball size
                    tarball_size = db["size"]

        else:
            tarball_url = None
            tarball_size = None

            for db in databases[args.database]:
                # Search for a specific version
                if db["version"] == args.database_version:
                    tarball_url = db["tarball"]

                    tarball_size = db["size"]

                    break

        if not tarball_url:
            raise Exception("Unable to find version \"{}\" for database ID \"{}\"".format(args.database_version, args.database))

    else:
        if args.database_file:
            if not os.path.isfile(args.database_file):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.database_file)

            tarball_filepath = args.database_file

            # Get file size in bytes
            tarball_size = "{}B".format(os.path.getsize(args.database_file))

        elif args.database_url:
            if not validate_url(args.database_url):
                raise ValueError("Invalid URL")
            
            tarball_url = args.database_url

            # We cannot know the database size from a URL a priori
            # Assume you have enough space to maintain the tarball and print a warning
            tarball_size = "1B"

            print("Warning: downloading a database from a URL could lead to quickly running out of storage!")

        else:
            raise Exception("No database selected")

    if not args.install_in:
        # Use the current working directory in case of no --install-in
        args.install_in = os.getcwd()

    # Check whether there is enough space to store and extract the database on the local filesystem
    # Get free space on disk in bytes
    _, _, free_space_bytes = shutil.disk_usage(args.install_in)

    # Compare the database size with the free space on disk
    tarball_size_bytes = convert_size(float(tarball_size[:-2]), tarball_size[-2:], "B")

    # Databases are uncompressed tarballs, so they are approximately the same size as their extracted content
    # Consider twice the space in case we have to download and extract the tarball
    size_factor = 1 if args.database_file else 2

    if free_space_bytes < tarball_size_bytes * size_factor:
        raise Exception("Not enough space to install the database")

    if not args.database_file:
        # Download the specified database
        print("Downloading database \"{}\" version \"{}\"".format(args.database, args.database_version))
        print(tarball_url)
        tarball_filepath = download(tarball_url, args.install_in)

        if not os.path.isfile(tarball_filepath):
            raise Exception("Unable to retrieve tarball {}".format(tarball_url))

    # Extract the database in --install-in
    with tarfile.open(tarball_filepath, "r") as tf:
        print("Unpacking tarball into {}".format(args.install_in))
        tf.extractall(path=args.install_in)

    # Fix file paths
    print("Fixing absolute paths")

    # kingdom, phylum, class, order, family, genus, species
    taxonomic_levels_prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

    for subdir, _, _ in os.walk(os.path.join(args.install_in, os.path.splitext(os.path.basename(tarball_filepath))[0])):
        dirname = os.path.basename(subdir)
        if dirname[:3] in taxonomic_levels_prefixes:
            # Edit file with the list of paths to the bloom filter files under this taxonomic level
            fix_paths(
                os.path.join(subdir, "{}.txt".format(dirname)),
                os.path.splitext(os.path.basename(tarball_filepath))[0],
                args.install_in
            )

            # Gunzip the bloom filter representation of this taxonomic level
            bf_filepath = os.path.join(subdir, "{}.bf.gz".format(dirname))
            
            with open(os.path.splitext(bf_filepath)[0], "w+") as bf_file:
                run(["gzip", "-dc", bf_filepath], stdout=bf_file, stderr=bf_file)

            # Edit the HowDeSBT index file
            fix_paths(
                os.path.join(subdir, "index", "index.detbrief.sbt"),
                os.path.splitext(os.path.basename(tarball_filepath))[0],
                args.install_in
            )


if __name__ == "__main__":
    main()
