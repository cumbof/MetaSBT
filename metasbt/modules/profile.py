#!/usr/bin/env python3
"""
Query a specific sequence bloom tree at all the seven taxonomic levels with an input genome 
or a file with a list of sequences, one per line. In case of an input genomes, results on 
single sequences are merged together
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 28, 2023"

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
        description=(
            "Query a specific sequence bloom tree at all the seven taxonomic levels with an input genome "
            "or a file with a list of sequences, one per line. In case of an input genomes, results on "
            "single sequences are merged together"
        ),
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
    p.add_argument(
        "--log",
        type=os.path.abspath,
        help="Path to the log file. Used to keep track of messages and errors printed on the stdout and stderr"
    )
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
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print messages and errors on the stdout"
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
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
    best_uncertainty: float = 25.0,
    logger: Optional[Logger] = None,
    verbose: bool = False,
) -> None:
    """
    Query the input genome against a specific tree
    Also expand the query to the lower taxonomic levels if requested

    :param input_file:          Path to the input file with query genome
    :param input_id:            Unique identifier of the input file
    :param tree:                Path to the tree definition file
    :param output_dir:          Path to the output folder
    :param threhsold:           Query threshold
    :param expand:              Expand the query to the lower taxonomic levels
    :param stop_at:             Stop expanding the query at a specific taxonomic level
    :param output_prefix:       Prefix of the output files with profiles
    :param best_uncertainty:    Consider this uncertainty percentage when selecting the best matches
    :param logger:              Logger object
    :param verbose:             Print messages on screen as alternative to the logger
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    printline("Input file: {}".format(input_file))
    printline("Input ID: {}".format(input_id))

    matches_dir = os.path.join(output_dir, input_id)
    os.makedirs(matches_dir, exist_ok=True)

    # Best matches and their scores
    best_matches = dict()

    # Path to the trees for keeping querying the database
    next_trees = [tree]

    while next_trees:
        curr_tree = next_trees.pop(0)
        
        printline("Querying {}".format(curr_tree))

        # Retrieve the taxonomic level name from the tree file path
        id_dir = os.path.dirname(curr_tree)
        level_dir = os.path.dirname(id_dir)
        level = os.path.basename(level_dir)

        levels = [level_dir.split(os.sep).index(level) for level in level_dir.split(os.sep) if level.startswith("k__")]
        curr_taxonomy = ""

        if levels:
            kingdom_index = levels[-1]
            curr_taxonomy = "|".join(level_dir.split(os.sep)[kingdom_index:])

        if level.startswith("s__") and expand:
            curr_tree = os.path.join(level_dir, "strains", "index", "index.detbrief.sbt")

        # Define the output file path
        output_file = os.path.join(matches_dir, "{}__{}__matches.txt".format(output_prefix, level))

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
                    "--tree={}".format(curr_tree),
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

                            node = line_split[0]
                            if level.startswith("s__"):
                                node = "t__{}".format(node)

                            if node not in matches_kmers:
                                matches_kmers[node] = {"common": 0, "total": 0}

                            # Update the number of common and total kmers for the current sequence
                            hits = line_split[1].split("/")
                            
                            matches_kmers[node]["common"] += int(hits[0])
                            matches_kmers[node]["total"] += int(hits[1])
                            matches_kmers[node]["score"] = round(matches_kmers[node]["common"] / matches_kmers[node]["total"], 2)
            
            # Search for the best match
            # The best match is defined as the genome with the highest score
            best_match = sorted(matches_kmers.keys(), key=lambda match: matches_kmers[match]["score"])[-1]
            best_score = matches_kmers[best_match]["score"]

            # Define a threshold on the best score
            score_threshold = float((best_score * best_uncertainty) / 100.0)

            # Get the best matches
            curr_best_matches = {"{}{}".format("{}|".format(curr_taxonomy) if curr_taxonomy else "", match): matches_kmers[match]
                                 for match in matches_kmers if matches_kmers[match]["score"] >= (best_score - score_threshold)}

            # Stop querying trees
            if not expand:
                # In case --expand is not provided
                break

            elif expand and stop_at:
                if stop_at[0] == level[0] and level != "strains":
                    # In case both --expand and --stop_at are provided
                    # and the taxonomic level specified with --stop_at is reached
                    # Compare the first character of the current taxonomic level and the one specified with --stop_at
                    break

            selected_best_matches = dict()

            same_level_matches = {
                match: best_matches[match] for match in best_matches if match.split("|")[-1].startswith(best_match.split("|")[-1][:3])
            }

            if same_level_matches:
                same_level_best = sorted(same_level_matches.keys(), key=lambda match: same_level_matches[match]["score"])[-1]
                same_level_best_score = same_level_matches[same_level_best]["score"]
                same_level_best_score_threshold = float((same_level_best_score * best_uncertainty) / 100.0)

            # Keep querying in case of --expand
            for match in curr_best_matches:
                process_match = False

                if not same_level_matches:
                    selected_best_matches[match] = curr_best_matches[match]
                    process_match = True

                elif curr_best_matches[match]["score"] >= (same_level_best_score - same_level_best_score_threshold):
                    selected_best_matches[match] = curr_best_matches[match]
                    process_match = True

                    if same_level_matches:
                        curr_best_score_threshold = float((curr_best_matches[match]["score"] * best_uncertainty) / 100.0)

                        if same_level_best_score < (curr_best_matches[match]["score"] - curr_best_score_threshold):
                            for same_level_match in list(best_matches.keys()):
                                if same_level_match.split("|")[-1].startswith(best_match.split("|")[-1][:3]):
                                    del best_matches[same_level_match]
                                    del same_level_matches[same_level_match]

                if process_match and not level.startswith("s__"):
                    # Convert match to level path
                    match_split = match.split("|")

                    if level in match_split:
                        match_split = match_split[match_split.index(level)+1:]

                    match_partial = os.sep.join(match_split)

                    # Add the path to the tree definition file to the tree bucket
                    next_trees.append(
                        os.path.join(
                            level_dir,
                            match_partial,
                            "index",
                            "index.detbrief.sbt"
                        )
                    )

            best_matches.update(selected_best_matches)

            if same_level_matches:
                best_matches.update(same_level_matches)

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

    if best_matches:
        for level in [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "t__",
        ]:
            best_level_matches = list(
                filter(
                    lambda match: match.split("|")[-1].startswith("{}__".format(level[0])),
                    best_matches.keys()
                )
            )

            best_level_match = sorted(best_level_matches, key=lambda match: best_matches[match]["score"])[-1]

            equally_best = list(
                filter(
                    lambda match: best_matches[match]["score"] == best_matches[best_level_match]["score"],
                    best_level_matches
                )
            )

            best_common = sorted(equally_best, key=lambda match: best_matches[match]["common"])[-1]

            printline(
                "Closest {}: {} (score {})".format(
                    level if not level.startswith("t__") else "genome",
                    best_common,
                    best_matches[best_common]["score"],
                )
            )

            # Report the characterisation to the output file
            with open(output_profile, "a+") as output:
                output.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        input_id,  # ID of the input query
                        level if not level.startswith("t") else "genome",  # Taxonomic level name
                        best_common,  # Lineage
                        "{}/{}".format(
                            best_matches[best_common]["common"],  # Number of kmers in common
                            best_matches[best_common]["total"],  # Total number of kmers in query
                        ),
                        best_matches[best_common]["score"],  # Score
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
