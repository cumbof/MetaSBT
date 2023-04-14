#!/usr/bin/env python3
"""
Create a MetaSBT database tarball
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 13, 2023"

import argparse as ap
import errno
import os
from pathlib import Path

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import run  # type: ignore
except Exception:
    pass

# Define the module name
TOOL_ID = "tar"

# Define the list of dependencies
DEPENDENCIES = ["tar"]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir",  # Install a database in a specific location
    "--out-file",  # Path to the output tarball file
]


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Create a MetaSBT database tarball",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        dest="db_dir",
        help="Path to the root folder of the MetaSBT database",
    )
    p.add_argument(
        "--out-file",
        type=os.path.abspath,
        dest="out_file",
        help="Path to the output tarball file",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the current "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def main() -> None:
    # Load command line parameters
    args = read_params()

    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    if os.path.isfile(args.out_file):
        raise Exception("The output tarball already exists")

    # kingdom, phylum, class, order, family, genus, species
    taxonomic_levels_prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

    # Take track of the paths to the compressed bloom filter representations of the taxonomic levels
    remove = list()

    for subdir, _, _ in os.walk(args.db_dir):
        dirname = os.path.basename(subdir)
        if dirname[:3] in taxonomic_levels_prefixes:
            # Compress the bloom filter representations of the taxonomic levels
            bf_filepath = os.path.join(subdir, "{}.bf".format(dirname))

            with open("{}.gz".format(bf_filepath), "w+") as bf_file:
                run(["gzip", "-c", bf_filepath], stdout=bf_file, stderr=bf_file)

            remove.append("{}.gz".format(bf_filepath))

            if dirname[0] == "s":
                # In case of species level, tag the genomes folders to be excluded from the tarball
                Path(os.path.join(subdir, "genomes", "exclude.tag")).touch()
                remove.append(os.path.join(subdir, "genomes", "exclude.tag"))

                # Do the same for the genomes folder under strains
                Path(os.path.join(subdir, "strains", "genomes", "exclude.tag")).touch()
                remove.append(os.path.join(subdir, "strains", "genomes", "exclude.tag"))

    # Also compress the database main bloom filter
    db_bf_filepath = os.path.join(args.db_dir, "{}.bf".format(os.path.basename(args.db_dir)))

    if os.path.isfile(db_bf_filepath):
        run(["gzip", db_bf_filepath], silence=True)

        remove.append("{}.gz".format(db_bf_filepath))

    # Create the tarball
    # Exclude the "genomes" folders and the uncompressed bloom filter files
    run(
        [
            "tar",
            "--exclude-tag-all='exclude.tag'",
            "--exclude='*.bf'",
            "-cf",
            args.out_file,
            args.db_dir
        ],
        silence=True
    )

    # Finally, remove unnecessary files
    for filepath in remove:
        os.unlink(filepath)


if __name__ == "__main__":
    main()
