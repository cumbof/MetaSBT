#!/usr/bin/env python3
"""Average-linkage hierarchical clustering based on ANI distances with Fastcluster.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.2"
__date__ = "Sep 12, 2024"

import argparse as ap
import errno
import math
import multiprocessing as mp
import os
import subprocess
import sys
from collections import Counter
from typing import List, Tuple

import fastcluster
import scipy.cluster.hierarchy as hier

import metasbt  # type: ignore

TOOL_ID = "cluster"

# Define the list of dependencies
DEPENDENCIES = [
    "gzip",
    "howdesbt",
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
        description="Average-linkage hierarchical clustering based on ANI distances with Fastcluster.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--filepath",
        type=os.path.abspath,
        required=True,
        help=(
            "Path to a two-columns TSV file with the list of input genome file paths (can be Gzip compressed) "
            "and their full taxonomic labels"
        ),
    )
    p.add_argument(
        "--kmer-size",
        type=int,
        default=21,
        dest="kmer_size",
        help="Kmer size for building Bloom Filter sketches",
    )
    p.add_argument(
        "--min-occurrences",
        type=int,
        default=21,
        dest="min_occurrences",
        help="Minimum number of kmer occurrences",
    )
    p.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Make it parallel",
    )
    p.add_argument(
        "--out-file",
        type=str,
        required=True,
        dest="out_file",
        help="Path to the output file with clusters assignments",
    )
    p.add_argument(
        "--sketch-size",
        type=int,
        default=10000,
        dest="sketch_size",
        help="Sketch size for building Bloom Filters",
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="ANI distance threshold for defining clusters with Fastcluster",
    )
    p.add_argument(
        "--tmp-dir",
        type=os.path.abspath,
        required=True,
        dest="tmp_dir",
        help="Path to the temporary folder",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def bf_distances(
    bf_sketches: List[os.path.abspath],
    kmer_size: int,
    tmp_dir: os.path.abspath,
    nproc: int=1,
) -> Tuple[List[str], List[float]]:
    """Compute the ANI distances between sketches, and build the condensed distance matrix.

    Parameters
    ----------
    bf_sketches : list
        List with paths to the BF sketch files.
    kmer_size : int
        The length of the k-mers.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    nproc : int, default 1
        Make the computation of BF distances parallel.

    Raises
    ------
    Exception
        In case an unexpected error occurred while computing BF distances.

    Returns
    -------
    tuple
        A tuple with the list of sequence names and the condensed distance matrix as list.
    """

    dist_folder = os.path.join(tmp_dir, "distances")

    if not os.path.isdir(dist_folder):
        os.makedirs(dist_folder, exist_ok=True)

    input_names = [os.path.splitext(os.path.basename(filepath))[0] for filepath in bf_sketches]

    condensed_distance_matrix = list()

    dists = dict()

    with mp.Pool(processes=nproc) as pool:
        jobs = [
            pool.apply_async(
                utils.ani_distance,
                args=(
                    bf_sketch_1,
                    kmer_size,
                    bf_sketches[bf_sketch_pos+1:],
                    dist_folder,
                )
            ) for bf_sketch_pos, bf_sketch_1 in enumerate(bf_sketches)
        ]

        for job in jobs:
            source, target_dists = job.get()

            dists[source] = target_dists

    for bf_sketch in bf_sketches:
        condensed_distance_matrix.extend(dists[bf_sketch])

    return input_names, condensed_distance_matrix


def bf_sketch(
    filepaths: List[os.path.abspath],
    tmp_dir: os.path.abspath,
    kmer_size: int=21,
    min_occurrences: int=2,
    sketch_size: int=10000,
    nproc: int=1
) -> List[os.path.abspath]:
    """Build Bloom Filter sketches.

    Parameters
    ----------
    filepaths : list
        List with paths to the input files.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    kmer_size : int, default 21
        Length of kmers for building BF sketches.
    min_occurrences : int, default 2
        Minimum number of kmer occurrences.
    sketch_size : int, default 10000
        Size of the BF sketch.
    nproc : int, default 1
        Make it parallel.

    Raises
    ------
    FileNotFoundError
        If a sketch file does not exist.
    Exception
        In case an unexpected error occurred while building a BF sketch file.

    Returns
    -------
    list
        A list of paths to the BF sketches.
    """

    bf_folder = os.path.join(tmp_dir, "sketches")

    if not os.path.isdir(bf_folder):
        os.makedirs(bf_folder, exist_ok=True)

    bf_filepaths = list()

    # TODO: Make it parallel
    for filepath in filepaths:
        # Assume these files are all unzipped, in fasta format, and with .fna extension
        bf_filepath = os.path.join(bf_folder, "{}.bf".format(os.path.splitext(os.path.basename(filepath))[0]))

        if not os.path.isfile(bf_filepath):
            try:
                # Build the BF sketch
                subprocess.check_call(
                    [
                        "howdesbt",
                        "makebf",
                        "--k={}".format(kmer_size),
                        "--min={}".format(min_occurrences),
                        "--bits={}".format(sketch_size),
                        "--hashes=1",
                        "--seed=0,0",
                        filepath,
                        "--out={}".format(bf_filepath),
                        "--threads={}".format(nproc),
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )

                if not os.path.isfile(bf_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), bf_filepath)

            except subprocess.CalledProcessError as e:
                raise Exception("An error has occurred while running BF sketch on {}".format(filepath)).with_traceback(e.__traceback__)

        bf_filepaths.append(bf_filepath)

    return bf_filepaths


def unzip_files(filepaths: List[os.path.abspath], tmp_dir: os.path.abspath) -> List[os.path.abspath]:
    """Unzip Gzip compressed files.

    Parameters
    ----------
    filepaths : list
        List with paths to the input files.
    tmp_dir : os.path.abspath
        Path to the temporary folder.

    Raises
    ------
    NotImplementedError
        If the one or more input files are not in fasta format (.fa, .fasta, or .fna).
    FileNotFoundError
        If the unzipped file does not exist.
    Exception
        In case of an unexpected error while running gzip.

    Returns
    -------
    list
        A list of paths to the unzipped files.
    """

    unzip_folder = os.path.join(tmp_dir, "unzip")

    if not os.path.isdir(unzip_folder):
        os.makedirs(unzip_folder, exist_ok=True)

    unzipped_filepaths = list()

    zipped_filepaths = list()

    for filepath in filepaths:
        if filepath.endswith(".gz"):
            # Only Gzip compressed files are supported here
            zipped_filepaths.append(filepath)

        elif filepath.endswith(".fa") or filepath.endswith(".fasta") or filepath.endswith(".fna"):
            # In case it is not Gzip compressed, consider it uncompressed
            # Only fasta files are supported here
            # Always use .fna
            if filepath.endswith(".fa") or filepath.endswith(".fasta"):
                fna_filepath = "{}.fna".format(os.path.splitext(filepath)[0])

                if not os.path.isfile(fna_filepath):
                    os.symlink(filepath, fna_filepath)

                unzipped_filepaths.append(fna_filepath)

            else:
                unzipped_filepaths.append(filepath)

        else:
            raise NotImplementedError("File format not supported\n{}".format(filepath))

    if zipped_filepaths:
        for filepath in zipped_filepaths:
            try:
                unzipped_filepath = os.path.join(unzip_folder, os.path.splitext(os.path.basename(filepath))[0])

                if not os.path.isfile(unzipped_filepath):
                    with open(unzipped_filepath, "w+") as outfile:
                        # Unzip file
                        subprocess.check_call(
                            [
                                "gzip",
                                "-dc",
                                filepath,
                            ],
                            stdout=outfile,
                            stderr=outfile
                        )

                    if not os.path.isfile(unzipped_filepath):
                        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), unzipped_filepath)

                unzipped_filepaths.append(unzipped_filepath)    

            except subprocess.CalledProcessError as e:
                raise Exception("An error has occurred while running gzip on {}".format(filepath)).with_traceback(e.__traceback__)

    return unzipped_filepaths


def main() -> None:
    args = read_params()

    if not os.path.isfile(args.filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.filepath)

    if os.path.isfile(args.out_file):
        raise Exception("The output file already exists!")

    if not os.path.isdir(args.tmp_dir):
        os.makedirs(args.tmp_dir, exist_ok=True)

    # Load the list of input file paths
    print("Loading the list of input files")

    if not os.path.isfile(args.filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.filepath)

    # Load the list of input genome file paths and their full taxonomic labels
    taxa = {
        line.strip().split("\t")[0]: line.strip().split("\t")[1]
            for line in open(args.filepath).readlines() if line.strip() and not line.strip().startswith("#")
    }

    input_dict = dict()

    for filepath in taxa:
        filename = os.path.splitext(os.path.basename(filepath))[0]

        if filepath.lower().endswith(".gz"):
            filename = os.path.splitext(filename)[0]

        input_dict[filename] = filepath

    input_filepaths = list(input_dict.values())

    # Unzip input files
    print("Unzipping file if Gzip compressed")
    input_filepaths = unzip_files(input_filepaths, args.tmp_dir)

    # Compute BF sketches
    print("Computing BF sketches")
    bf_sketches = bf_sketch(
        input_filepaths,
        args.tmp_dir,
        kmer_size=args.kmer_size,
        min_occurrences=args.min_occurrences,
        sketch_size=args.sketch_size,
        nproc=args.nproc,
    )

    # Compute the BF distances pair-wise
    print("Computing BF ANI distances and building the condensed similarity matrix")
    input_names, condensed_matrix = bf_distances(bf_sketches, args.kmer_size, args.tmp_dir, args.nproc)

    # Run fastcluster
    print("Computing a average-linkage hierarchical clustering")
    dendro = fastcluster.linkage(condensed_matrix, method="average")

    # Define clusters with threshold
    print("Defining clusters with threshold {}".format(args.threshold))
    clusters = hier.fcluster(dendro, args.threshold, criterion="distance")

    clusters_dict = dict()

    # Define a mapping cluster - input files
    for filename, cluster in zip(input_names, clusters):
        if cluster not in clusters_dict:
            clusters_dict[cluster] = list()

        clusters_dict[cluster].append(input_dict[filename])

    taxa_clusters = dict()

    for cluster in clusters_dict:
        taxonomies = sorted([taxa[filepath] for filepath in clusters_dict[cluster]])

        counter = Counter(taxonomies)

        # Get the most occurrent taxonomic label in cluster
        assignment = counter.most_common(1)[0][0]

        if assignment not in taxa_clusters:
            taxa_clusters[assignment] = list()

        taxa_clusters[assignment].append(cluster)

    clusters_taxa_fixed = dict()

    # Fix cluster names if a taxonomy is assigned to more than one cluster
    for taxonomy in taxa_clusters:
        if len(taxa_clusters[taxonomy]) == 1:
            clusters_taxa_fixed[taxa_clusters[taxonomy][0]] = taxonomy

        else:
            for i, cluster in enumerate(taxa_clusters[taxonomy]):
                clusters_taxa_fixed[cluster] = "{}__clade_{}".format(taxonomy, i)

    # Write assignments
    with open(args.out_file, "w+") as outfile:
        for cluster in clusters_dict:
            for filepath in clusters_dict[cluster]:
                outfile.write("{}\t{}\n".format(filepath, clusters_taxa_fixed[cluster]))


if __name__ == "__main__":
    main()
