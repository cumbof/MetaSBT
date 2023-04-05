#!/usr/bin/env python3
"""
Build minimal bloom filter sketches with cluster-specific marker kmers 
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 3, 2023"

import argparse as ap
import errno
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List


def read_params():
    p = ap.ArgumentParser(
        prog="bf_sketch",
        description="Build minimal bloom filter sketches with cluster-specific marker kmers",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required=True,
        dest="db_dir",
        help="Path to the root folder of a MetaSBT database",
    )
    p.add_argument(
        "--level",
        type=str,
        choices=["phylum", "class", "order", "family", "genus", "species"],
        help=(
            "Build a sketch for each of the clusters under a particular taxonomic level "
            "(e.g., \"--level genus\" will build a sketch for every genus in the database)"
        )
    )
    p.add_argument(
        "--out-dir",
        type=os.path.abspath,
        required=True,
        dest="out_dir",
        help="Path to the output folder",
    )
    p.add_argument(
        "--taxonomy",
        type=str,
        help=(
            "Build the sketch for this specific full or partial taxonomy. "
            "It must start with a kingdom up to a species. "
            "Note that it does not necessarily end with a species. "
            "Look at the following examples as valid taxonomic labels: "
            "k__Viruses|p__Nucleocytoviricota|c__Pokkesviricetes|o__Chitovirales|f__Poxviridae|g__Orthopoxvirus|s__Monkeypox_virus, "
            "k__Viruses|p__Nucleocytoviricota|c__Pokkesviricetes|o__Chitovirales|f__Poxviridae, "
            "k__Viruses|p__Nucleocytoviricota"
        ),
    )
    return p.parse_args()


def build_sketch(bf_filepath: str, bf_bucket: List[str], out_dir: str) -> str:
    """
    Build a bloom filter sketch with marker kmers

    :param bf_filepath: Path to the input bloom filter file under a specific taxonomic level
    :param bf_bucket:   List with paths to the bloom filter files under the same taxonomic level
                        of the input bf_filepath
    :param out_dir:     Path to the output folder
    :return:            Path to the output file
    """

    # Remove the input bf_filepath from the bucket
    bf_files_bucket = list(set(bf_bucket).difference({bf_filepath}))

    bf_filename = os.path.splitext(os.path.basename(bf_filepath))[0]

    # Path to the bloom filter as the result of the bitwise OR operator
    # on all the bloom filters in the input bf_bucket
    sum_bf_filepath = os.path.join(out_dir, "{}__sum.bf".format(bf_filename))

    # Path to the bloom filter as the result of the NOT operator
    # on the result of the OR operator on the bloom filters in the bf_bucket
    not_bf_filepath = os.path.join(out_dir, "{}__not.bf".format(bf_filename))

    # Path to the output bloom filter as the result of the bitwise AND operator
    # between bf_filepath and not_bf_filepath
    diff_bf_filepath = os.path.join(out_dir, "{}.bf".format(bf_filename))

    # Given a input node A and a bucket with nodes B, C, and D
    # the following pipeline produces a node X as the result
    # of the following operations:
    # X = A - (B + C + D)
    # X = A and (not (B or C or D))
    try:
        with tempfile.NamedTemporaryFile() as bflist:
            with open(bflist.name, "wt") as bfilters:
                for bf_path in bf_files_bucket:
                    # Check whether the bf file exists
                    if not os.path.isfile(bf_path):
                        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), bf_path)

                    # Write the list of bloom filter file paths in bf_bucket
                    bfilters.write("{}\n".format(bf_path))

            # Bitwise OR on every bloom filter in the bucket
            # Called (sum)
            subprocess.check_call(
                [
                    "howdesbt",
                    "bfoperate",
                    "--list={}".format(bflist.name),
                    "--or",
                    "--out={}".format(sum_bf_filepath)
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

        # Comlement NOT of the (sum) bloom filter
        # Called (not)
        subprocess.check_call(
            [
                "howdesbt",
                "bfoperate",
                sum_bf_filepath,
                "--not",
                "--out={}".format(not_bf_filepath)
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        # Output bloom filter as the bitwise AND between
        # the input bloom filter and (not)
        subprocess.check_call(
            [
                "howdesbt",
                "bfoperate",
                bf_filepath,
                not_bf_filepath,
                "--and",
                "--out={}".format(diff_bf_filepath)
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

    except subprocess.CalledProcessError as e:
        raise Exception("An error has occurred while running howdesbt").with_traceback(e.__traceback__)

    # Delete intermediate bloom filters
    if os.path.isfile(sum_bf_filepath):
        os.unlink(sum_bf_filepath)

    if os.path.isfile(not_bf_filepath):
        os.unlink(not_bf_filepath)

    return diff_bf_filepath


def main() -> None:
    args = read_params()

    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    if args.taxonomy:
        tax_folder = os.path.join(args.db_dir, args.taxonomy.replace("|", os.sep))

        if not os.path.isdir(tax_folder):
            raise ValueError("The specified taxonomy does not exist in the database")
        
        # Get the cluster level
        bf_level_id = os.path.basename(tax_folder)[0]
    
    elif args.level:
        # Get the cluster level
        bf_level_id = args.level[0]
    
    # Search for cluster at the same taxonomic level
    clusters = Path(args.db_dir).glob("**/{}__*".format(bf_level_id))

    bf_bucket = [
        os.path.join(str(cluster), "{}.bf".format(os.path.basename(str(cluster)))) \
            for cluster in clusters if os.path.isdir(str(cluster))
    ]

    if not bf_bucket:
        raise Exception("No bloom filter found in the database")

    if args.taxonomy:
        # Build a bloom filter sketch for a specific cluster
        bf_filepath = os.path.join(tax_folder, "{}.bf".format(os.path.basename(tax_folder)))

        # Finally build the sketch
        ouf_file = build_sketch(bf_filepath, bf_bucket, args.out_dir)

        print(ouf_file)

    elif args.level:
        # Build a bloom filter sketch for each cluster under a specific taxonomic level
        for bf_filepath in bf_bucket:
            out_file = build_sketch(bf_filepath, bf_bucket, args.out_dir)

            print(out_file)


if __name__ == "__main__":
    main()
