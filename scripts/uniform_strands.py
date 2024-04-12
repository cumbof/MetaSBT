#!/usr/bin/env python3
"""Uniform the strands of the input sequences.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 11, 2024"

import argparse as ap
import errno
import multiprocessing as mp
import os
import sys

import metasbt  # type: ignore

TOOL_ID = "uniform_strands"

# Define the list of dependencies
DEPENDENCIES = [
    "metasbt",
]

# Assume MetaSBT is installed
# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(metasbt.__file__), "modules")

# Define the path to utils.py
UTILS_FILEPATH = os.path.join(MODULES_DIR, "utils.py")

# Check whether utils.py exists
if not os.path.isfile(UTILS_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), UTILS_FILEPATH)

# Add the modules dir to the system path
sys.path.append(MODULES_DIR)

# Finally import utils
import utils


def read_params():
    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Uniform the strands of the input sequences.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--filepath",
        type=os.path.abspath,
        required=True,
        help="File with the list of paths to the input FASTA files",
    )
    p.add_argument(
        "--out-dir",
        type=os.path.abspath,
        required=True,
        dest="out_dir",
        help="Path to the output folder",
    )
    p.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Process the input files in parallel",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def main() -> None:
    args = read_params()

    if not os.path.isfile(args.filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.filepath)
    
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    if args.nproc < 1:
        args.nproc = 1

    if args.nproc > os.cpu_count():
        args.nproc = os.cpu_count()

    filepaths = [line.strip() for line in open(args.filepath).readlines() if line.strip() and not line.startswith("#")]

    if not filepaths:
        raise Exception("No input files available!")

    with mp.Pool(processes=args.nproc) as pool:
        jobs = [pool.apply_async(utils.strand, args=(str(filepath), args.out_dir,)) for filepath in filepaths]

        for job in jobs:
            job.get()


if __name__ == "__main__":
    main()
