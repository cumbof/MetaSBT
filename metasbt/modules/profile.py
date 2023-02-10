#!/usr/bin/env python3
"""
Query a specific sequence bloom tree with an input genome or file with a list of sequences, one per line.
In case of input genomes, results on single sequences are merged together
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Feb 8, 2023"

import argparse as ap
import errno
import os
import time
from functools import partial
from logging import Logger
from typing import Optional

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import get_file_info, init_logger, number, println, run  # type: ignore
except Exception:
    pass

# Define the module name
TOOL_ID = "profile"

# Define the list of dependencies
DEPENDENCIES = ["howdesbt"]

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--input-file",  # Path to the input file
    "--log",  # Path to the log file
    "--output-dir",  # Output folder path
    "--tree",  # Path to the tree definition file
]


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Query the sequence bloom trees at all the 7 taxonomic levels",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--expand",
        action="store_true",
        default=False,
        help="Expand the input query on all the taxonomic levels",
    )
    p.add_argument(
        "--input-file",
        type=os.path.abspath,
        required=True,
        dest="input_file",
        help="Path to the input query",
    )
    p.add_argument(
        "--input-id",
        type=str,
        dest="input_id",
        help="Unique identifier of the input query",
    )
    p.add_argument(
        "--input-type",
        type=str,
        required=True,
        dest="input_type",
        choices=["genome", "list"],
        help=(
            "Accepted input types are genomes and files with one nucleotide sequence per line. "
            "In case of genomes, if they contain multiple sequences, results are merged together"
        ),
    )
    p.add_argument("--log", type=os.path.abspath, help="Path to the log file")
    p.add_argument(
        "--output-dir",
        type=os.path.abspath,
        required=True,
        dest="output_dir",
        help="This is the output folder with queries results",
    )
    p.add_argument(
        "--output-prefix",
        type=str,
        dest="output_prefix",
        help="Prefix of the output files with query matches",
    )
    p.add_argument(
        "--stop-at",
        type=str,
        dest="stop_at",
        choices=["phylum", "class", "order", "family", "genus"],
        help=(
            "Stop expanding queries at a specific taxonomic level. "
            "Please note that this argument works in conjunction with --expand only. "
            "Available values: phylum, class, order, family, genus"
        ),
    )
    p.add_argument(
        "--threshold",
        type=number(float, minv=0.0, maxv=1.0),
        default=0.0,
        help=(
            "Fraction of query kmers that must be present in a leaf to be considered a match. "
            "This must be between 0.0 and 1.0"
        ),
    )
    p.add_argument(
        "--tree",
        type=os.path.abspath,
        required=True,
        help="This is the tree definition file",
    )
    p.add_argument("--verbose", action="store_true", default=False, help="Print results on screen")
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the current "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def profile_list(
    input_file: str,
    input_id: str,
    tree: str,
    output_dir: str,
    threshold: float = 0.0,
    output_prefix: Optional[str] = None,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> None:
    """
    Query a list of sequences against a specific tree

    :param input_file:      Path to the input file with query sequences
    :param input_id:        Unique identifier of the input file
    :param tree:            Path to the tree definition file
    :param output_dir:      Path to the output folder
    :param threhsold:       Query threshold
    :param output_prefix:   Prefix of the output files with profiles
    :param logger:          Logger object
    :param verbose:         Print messages on screen as alternative to the logger
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    printline("Input file: {}".format(input_file))
    printline("Input ID: {}".format(input_id))

    # Define the output file path
    output_file = os.path.join(output_dir, "{}__matches.txt".format(output_prefix))

    # Run HowDeSBT
    with open(output_file, "w+") as file:
        run(
            [
                "howdesbt",
                "query",
                "--sort",
                "--distinctkmers",
                "--tree={}".format(tree),
                "--threshold={}".format(threshold),
                input_file,
            ],
            stdout=file,
            stderr=file,
        )

    printline("Output: {}".format(output_file))


def profile_genome(
    input_file: str,
    input_id: str,
    tree: str,
    output_dir: str,
    threshold: float = 0.0,
    expand: bool = False,
    stop_at: Optional[str] = None,
    output_prefix: Optional[str] = None,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> None:
    """
    Query the input genome against a specific tree
    Also expand the query to the lower taxonomic levels if requested

    :param input_file:      Path to the input file with query genome
    :param input_id:        Unique identifier of the input file
    :param tree:            Path to the tree definition file
    :param output_dir:      Path to the output folder
    :param threhsold:       Query threshold
    :param expand:          Expand the query to the lower taxonomic levels
    :param stop_at:         Stop expanding the query at a specific taxonomic level
    :param output_prefix:   Prefix of the output files with profiles
    :param logger:          Logger object
    :param verbose:         Print messages on screen as alternative to the logger
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    printline("Input file: {}".format(input_file))
    printline("Input ID: {}".format(input_id))

    # Define the lineage
    lineage = dict()
    # Also define the closest genome
    closest_genome = ""
    closest_genome_common_kmers = 0
    closest_genome_total_kmers = 0
    closest_genome_score = 0.0

    while os.path.isfile(tree):
        printline("Querying {}".format(tree))

        # Retrieve the taxonomic level name from the tree file path
        id_dir = os.path.dirname(tree)
        level_dir = os.path.dirname(id_dir)
        level = os.path.basename(level_dir)

        # Define the output file path
        output_file = os.path.join(output_dir, "{}__{}__matches.txt".format(output_prefix, level))

        # There is no need to run HowDeSBT in case of clusters with only one node/genome
        # The result is the same of profiling the higher level
        # TODO

        # Run HowDeSBT
        with open(output_file, "w+") as file:
            run(
                [
                    "howdesbt",
                    "query",
                    "--sort",
                    "--distinctkmers",
                    "--tree={}".format(tree),
                    "--threshold={}".format(threshold),
                    input_file,
                ],
                stdout=file,
                stderr=file,
            )

        if os.path.isfile(output_file):
            # Take track of common and total kmers
            matches_kmers = dict()

            # Read the output file
            with open(output_file) as output:
                for line in output:
                    line = line.strip()
                    if line:
                        if not line.startswith("#") and not line.startswith("*"):
                            line_split = line.split(" ")
                            if line_split[0] not in matches_kmers:
                                matches_kmers[line_split[0]] = {"common": 0, "total": 0}
                            # Update the number of common and total kmers for the current sequence
                            hits = line_split[1].split("/")
                            matches_kmers[line_split[0]]["common"] += int(hits[0])
                            matches_kmers[line_split[0]]["total"] += int(hits[1])

            # Search for the best match
            # The best match is defined as the genome with more kmers in common with the input sequence
            best_match = ""
            best_kmers = 0
            for match in matches_kmers:
                if matches_kmers[match]["common"] > best_kmers:
                    best_kmers = matches_kmers[match]["common"]
                    best_match = match

            # Retrieve the score of the best match
            best_score = round(best_kmers / matches_kmers[best_match]["total"], 2)

            if level.startswith("s__"):
                # There is nothing to query after the species tree
                closest_genome = best_match
                closest_genome_common_kmers = best_kmers
                closest_genome_total_kmers = matches_kmers[best_match]["total"]
                closest_genome_score = best_score

                break

            else:
                # Update the lineage
                lineage[best_match] = {
                    "common": best_kmers,
                    "total": matches_kmers[best_match]["total"],
                    "score": best_score,
                }

                # Stop querying trees
                if not expand:
                    # In case --expand is not provided
                    break

                elif expand and stop_at:
                    if stop_at[0] == level[0]:
                        # In case both --expand and --stop_at are provided
                        # and the taxonomic level specified with --stop_at is reached
                        # Compare the first character of the current taxonomic level and the one specified with --stop_at
                        break

                # Keep querying in case of --expand
                tree = os.path.join(level_dir, best_match, "index", "index.detbrief.sbt")

    # Define the output profile path
    output_profile = os.path.join(output_dir, "{}__profiles.tsv".format(output_prefix))

    if not os.path.isfile(output_profile):
        # Write the header lines in case the output profile does not exist
        with open(output_profile, "w+") as output:
            output.write("# Tree: {}\n".format(tree))
            output.write(
                "# {}\t{}\t{}\t{}\t{}\n".format(
                    "Input ID",
                    "Level",
                    "Taxonomy",
                    "Common kmers",
                    "Score",
                )
            )

    # Reconstruct the lineage
    if lineage:
        # Print the closest lineage
        printline("Closest lineage: {}".format("|".join(lineage.keys())))

        # Print scores
        printline("Listing scores")
        for lin in lineage:
            printline("{}: {}".format(lin, lineage[lin]["score"]))

            # Expand the level ID to the full level name
            level_name = ""
            for tx_level in [
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]:
                if lin[0] == tx_level[0]:
                    level_name = tx_level
                    break

            # Report the characterisation to the output file
            with open(output_profile, "a+") as output:
                output.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        input_id,  # ID of the input query
                        level_name,  # Taxonomic level name
                        lin,  # Lineage
                        "{}/{}".format(
                            lineage[lin]["common"],  # Number of kmers in common
                            lineage[lin]["total"],  # Total number of kmers in query
                        ),
                        lineage[lin]["score"],  # Score
                    )
                )

    # Show the closest genome in case the input is a species tree
    # or the query has been expanded
    if closest_genome:
        # Print the closest genome and its score
        printline("Closest genome: {}".format(closest_genome))
        printline("Closest genome score: {}".format(closest_genome_score))

        # Report the closest genome to the output file
        with open(output_profile, "a+") as output:
            output.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    input_id,  # ID of the input query
                    "genome",  # Genome level
                    closest_genome,  # Closest genome in the database
                    "{}/{}".format(
                        closest_genome_common_kmers,  # Number of kmers in common
                        closest_genome_total_kmers,
                    ),  # Total number of kmers in query
                    closest_genome_score,  # Score
                )
            )

    printline("Output: {}".format(output_profile))


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    # Check whether the input files and folders exist on the file system
    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.input_file)

    if not os.path.isfile(args.tree):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.tree)

    # Create the output folder
    os.makedirs(args.output_dir, exist_ok=True)

    # In case --input_id is not provided
    if not args.input_id:
        # Use the input file name
        _, args.input_id, _, _ = get_file_info(args.input_id)

    # In case --output-prefix is not provided
    if not args.output_prefix:
        # Use the input_id
        args.output_prefix = args.input_id

    t0 = time.time()

    if args.input_type == "list":
        # In case the input file contains a list of sequences
        profile_list(
            args.input_file,
            args.input_id,
            args.tree,
            args.output_dir,
            threshold=args.threshold,
            output_prefix=args.output_prefix,
            logger=logger,
            verbose=args.verbose,
        )

    elif args.input_type == "genome":
        # In case of input genome
        if args.stop_at and not args.expand:
            # Raise an exception in case --expand is not provided
            raise Exception("The --stop-at argument must me used in conjunction with --expand only")

        profile_genome(
            args.input_file,
            args.input_id,
            args.tree,
            args.output_dir,
            threshold=args.threshold,
            expand=args.expand,
            stop_at=args.stop_at,
            output_prefix=args.output_prefix,
            logger=logger,
            verbose=args.verbose,
        )

    t1 = time.time()
    println(
        "Total elapsed time {}s".format(int(t1 - t0)),
        logger=logger,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
