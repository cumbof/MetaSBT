#!/usr/bin/env python3
"""Create a MetaSBT database tarball.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Mar 21, 2024"

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
    """Read and test the input arguments.

    Returns
    -------
    argparse.ArgumentParser
        The ArgumentParser object
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

    if args.db_dir.endswith(os.sep):
        # Trim the last char out of the string
        args.db_dir = args.db_dir[:-1]

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

            with open(os.path.join(subdir, ".tagignore"), "w+") as exclude_tag:
                exclude_tag.write("{}.bf\n".format(dirname))
                exclude_tag.write("howdesbt.log")

            remove.append("{}.gz".format(bf_filepath))
            remove.append(os.path.join(subdir, ".tagignore"))

            if dirname[0] == "s":
                if os.path.isdir(os.path.join(subdir, "strains")):
                    with open(os.path.join(subdir, "strains", ".tagignore"), "a+") as exclude_tag:
                        exclude_tag.write("strains.bf\n")
                        exclude_tag.write("tmp")

                    remove.append(os.path.join(subdir, "strains", ".tagignore"))

                    # Compress the bloom filter representation of the strains
                    strains_bf_filepath = os.path.join(subdir, "strains", "strains.bf")

                    if os.path.isfile(strains_bf_filepath):
                        with open("{}.gz".format(strains_bf_filepath), "w+") as bf_file:
                            run(["gzip", "-c", strains_bf_filepath], stdout=bf_file, stderr=bf_file)

                    if os.path.isdir(os.path.join(subdir, "strains", "genomes")):
                        with open(os.path.join(subdir, "strains", "genomes", ".tagignore"), "a+") as exclude_tag:
                            exclude_tag.write("*")

                        remove.append(os.path.join(subdir, "strains", "genomes", ".tagignore"))

                    if os.path.isdir(os.path.join(subdir, "strains", "filters")):
                        # Compress all the bloom filters under the filters folder at the strain level
                        for bf_filepath in Path(os.path.join(subdir, "strains", "filters")).glob("*.bf"):
                            with open("{}.gz".format(bf_filepath), "w+") as bf_file:
                                run(["gzip", "-c", bf_filepath], stdout=bf_file, stderr=bf_file)
                            
                            remove.append("{}.gz".format(bf_filepath))
                    
                        with open(os.path.join(subdir, "strains", "filters", ".tagignore"), "a+") as exclude_tag:
                            exclude_tag.write("*.bf")

                        remove.append(os.path.join(subdir, "strains", "filters", ".tagignore"))

                with open(os.path.join(subdir, "genomes", ".tagignore"), "a+") as exclude_tag:
                    exclude_tag.write("*")

                # Compress all the bloom filters under the filters folder at the species level
                for bf_filepath in Path(os.path.join(subdir, "filters")).glob("*.bf"):
                    with open("{}.gz".format(bf_filepath), "w+") as bf_file:
                        run(["gzip", "-c", bf_filepath], stdout=bf_file, stderr=bf_file)
                    
                    remove.append("{}.gz".format(bf_filepath))

                with open(os.path.join(subdir, "filters", ".tagignore"), "a+") as exclude_tag:
                    exclude_tag.write("*.bf")
                
                remove.append("{}.gz".format(strains_bf_filepath))
                remove.append(os.path.join(subdir, "genomes", ".tagignore"))
                remove.append(os.path.join(subdir, "filters", ".tagignore"))

    with open(os.path.join(args.db_dir, ".tagignore"), "w+") as exclude_tag:
        exclude_tag.write("*.log\n")
        exclude_tag.write("*.sh\n")

    # Also compress the database main bloom filter file
    db_bf_filepath = os.path.join(args.db_dir, "{}.bf".format(os.path.basename(args.db_dir)))

    if os.path.isfile(db_bf_filepath):
        with open("{}.gz".format(db_bf_filepath), "w+") as bf_file:
            run(["gzip", "-c", db_bf_filepath], stdout=bf_file, stderr=bf_file)

        with open(os.path.join(args.db_dir, ".tagignore"), "a+") as exclude_tag:
            exclude_tag.write("{}.bf".format(os.path.basename(args.db_dir)))

        remove.append("{}.gz".format(db_bf_filepath))

    remove.append(os.path.join(args.db_dir, ".tagignore"))

    # Create the tarball
    run(
        [
            "tar",
            "-cf",
            args.out_file,
            "--exclude-ignore=.tagignore",
            "-C",
            os.path.dirname(args.db_dir),
            os.path.basename(args.db_dir)
        ],
        silence=True
    )

    # Finally, remove unnecessary files
    for filepath in remove:
        os.unlink(filepath)


if __name__ == "__main__":
    main()
