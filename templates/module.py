#!/usr/bin/env python3
"""
Describe your module here and locate this file under the modules folder
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Feb 8, 2023"

import argparse as ap
import os
import time
from typing import List

# Local modules are not available when the main controller
# tries to load them for accessing their variables
try:
    # Load utility functions
    from utils import init_logger, println  # type: ignore
except Exception:
    pass

# Define the module name
TOOL_ID = "template"

# Define the list of dependencies
DEPENDENCIES: List[str] = list()

# Define the list of input files and folders
FILES_AND_FOLDERS: List[str] = list()


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Describe your module here",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )

    # TODO Add other parameters

    p.add_argument("--log", type=os.path.abspath, help="Path to the log file")
    p.add_argument("--verbose", action="store_true", default=False, help="Print results on screen")
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

    # Initialise the logger
    logger = init_logger(filepath=args.log, toolid=TOOL_ID, verbose=args.verbose)

    t0 = time.time()

    # TODO Add your code here

    t1 = time.time()
    println(
        "Total elapsed time {}s".format(int(t1 - t0)),
        logger=logger,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
