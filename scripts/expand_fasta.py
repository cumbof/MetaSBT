#!/usr/bin/env python3
"""Expand a FASTA file. One file for each read.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 4, 2024"

import argparse as ap
import errno
import os

TOOL_ID = "expand_fasta"


def read_params():
    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Expand a FASTA file. One file for each read.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--filepath",
        type=os.path.abspath,
        required=True,
        help="Path to the input FASTA file",
    )
    p.add_argument(
        "--out-dir",
        type=os.path.abspath,
        required=True,
        dest="out_dir",
        help="Path to the output folder",
    )
    p.add_argument(
        "--out-file-prefix",
        type=str,
        default="Read",
        dest="out_file_prefix",
        help="Output file prefix followed by an incremental numerical ID",
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

    out_file_prefix = "{}_".format(args.out_file_prefix) if args.out_file_prefix else ""

    with open(args.filepath) as infile:
        outfile = None
        filecount = 0

        for line in infile:
            if line.startswith(">"):
                if outfile:
                    # Close the output file buffer
                    outfile.close()

                # Just an incremental numerical ID
                filecount += 1

                # Open a new output file buffer
                outfile = open(os.path.join(args.out_dir, "{}{}.fna".format(args.out_file_prefix, filecount)))

            # Write the line to the output file
            outfile.write(line)

        # Close the last output file buffer
        outfile.close()

if __name__ == "__main__":
    main()
