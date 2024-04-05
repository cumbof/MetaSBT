#!/usr/bin/env python3
"""Average-linkage hierarchical clustering based on MASH genomic distances with Fastcluster.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 4, 2024"

import argparse as ap
import errno
import os
import subprocess
from typing import List, Tuple

import fastcluster
import scipy.cluster.hierarchy as hier

TOOL_ID = "fastcluster"

# Define the list of dependencies
DEPENDENCIES = [
    "gzip",
    "mash",
]


def read_params():
    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Average-linkage hierarchical clustering based on MASH genomic distances with Fastcluster.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--filepath",
        type=os.path.abspath,
        required=True,
        help="Path to the input file with the list of paths to the input genomes (can be Gzip compressed)",
    )
    p.add_argument(
        "--kmer-size",
        type=int,
        default=21,
        dest="kmer_size",
        help="Kmer size for building MASH sketches",
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
        help="Sketch size for building MASH sketches",
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="MASH distance threshold for defining clusters with Fastcluster",
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


def mash_distance(
    mash_sketches: List[os.path.abspath],
    mash_pastes: List[os.path.abspath],
    tmp_dir: os.path.abspath,
    nproc: int=1,
) -> Tuple[List[str], List[float]]:
    """Compute the MASH distances between sketches and pastes, and build the condensed distance matrix.

    Parameters
    ----------
    mash_sketches : list
        List with paths to the MASH sketch files.
    mash_pastes : list
        List with paths to the MASH paste files.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    nproc : int, default 1
        Make the computation of MASH distances parallel.

    Raises
    ------
    Exception
        In case an unexpected error occurred while computing MASH dist.

    Returns
    -------
    tuple
        A tuple with the list of sequence names and the condensed distance matrix as list.
    """

    mash_folder = os.path.join(tmp_dir, "distances")

    if not os.path.isdir(mash_folder):
        os.makedirs(mash_folder, exist_ok=True)

    input_names = [os.path.splitext(os.path.basename(filepath))[0] for filepath in mash_sketches]

    condensed_distance_matrix = list()

    for mash_sketch_pos, mash_sketch in enumerate(mash_sketches):
        distance_filepath = os.path.join(mash_folder, "{}.tsv".format(os.path.splitext(os.path.basename(mash_sketch))[0]))

        if not os.path.isfile(distance_filepath):
            # TODO: A distance file cannot be recovered in this way.
            #       It should take track of the last mash paste file first.
            for mash_paste in mash_pastes:
                try:
                    with open(distance_filepath, "a+") as dist_file:
                        subprocess.check_call(
                            [
                                "mash",
                                "dist",
                                "-p",
                                str(nproc),
                                "-t",
                                mash_sketch,
                                mash_paste,
                            ],
                            stdout=dist_file,
                            stderr=subprocess.DEVNULL
                        )

                except subprocess.CalledProcessError as e:
                    raise Exception("An error has occurred while running MASH dist on {}".format(mash_sketch)).with_traceback(e.__traceback__) 

        distances_dict = {
            os.path.splitext(os.path.basename(line.strip().split("\t")[0]))[0]: float(line.strip().split()[1])
            for line in open(distance_filepath).readlines() if line.strip()
        }

        # Sort distances according to input_names
        distances = [distances_dict[input_name] for input_name in input_names]

        # Populate the condensed distance matrix
        condensed_distance_matrix.extend(distances[mash_sketch_pos+1:])

    return input_names, condensed_distance_matrix


def mash_paste(filepaths: List[os.path.abspath], tmp_dir: os.path.abspath) -> os.path.abspath:
    """Paste MASH sketches.

    Parameters
    ----------
    filepaths : list
        List of paths to the MASH sketch files..
    tmp_dir : os.path.abspath
        Path to the temporary folder.

    Raises
    ------
    FileNotFoundError
        If a paste file does not exist.
    Exception
        In case an unexpected error occurred while building a MASH paste file.

    Returns
    -------
    list
        A list of paths to the MASH pastes.
    """

    mash_folder = os.path.join(tmp_dir, "pastes")

    if not os.path.isdir(mash_folder):
        os.makedirs(mash_folder, exist_ok=True)

    mash_pastes = list()

    # Split the list of input file paths into N chunks of up to 50000 elements each
    filepaths_chunks = [filepaths[i:i+50000] for i in range(0, len(filepaths), 50000)]

    for chunk_count, filepaths_chunk in enumerate(filepaths_chunks):
        sketches_filepath = os.path.join(mash_folder, "paste_{}.txt".format(chunk_count))

        if not os.path.isfile(sketches_filepath):
            # Write a file with the list of sketches file paths to paste
            with open(sketches_filepath, "w+") as sketches_paste_list:
                for sketch_filepath in filepaths_chunk:
                    sketches_paste_list.write("{}\n".format(sketch_filepath))

        mash_paste_filepath = os.path.join(
            os.path.dirname(sketches_filepath),
            "{}.msh".format(os.path.splitext(os.path.basename(sketches_filepath))[0])
        )

        if not os.path.isfile(mash_paste_filepath):
            try:
                # Paste the MASH sketches
                subprocess.check_call(
                    [
                        "mash",
                        "paste",
                        "-l",
                        os.path.splitext(mash_paste_filepath)[0],
                        sketches_filepath,
                    ],
                    stderr=subprocess.DEVNULL
                )

                if not os.path.isfile(mash_paste_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mash_paste_filepath)

            except subprocess.CalledProcessError as e:
                raise Exception("An error has occurred while running MASH paste").with_traceback(e.__traceback__)

        mash_pastes.append(mash_paste_filepath)

    return mash_pastes


def mash_sketch(
    filepaths: List[os.path.abspath],
    tmp_dir: os.path.abspath,
    kmer_size: int=21,
    sketch_size: int=10000,
) -> List[os.path.abspath]:
    """Build MASH sketches.

    Parameters
    ----------
    filepaths : list
        List with paths to the input files.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    kmer_size : int, default 21
        Length of kmers for building MASH sketches.
    sketch_size : int, default 10000
        Size of the MASH sketch.

    Raises
    ------
    FileNotFoundError
        If a sketch file does not exist.
    Exception
        In case an unexpected error occurred while building a MASH sketch file.

    Returns
    -------
    list
        A list of paths to the MASH sketches.
    """

    mash_folder = os.path.join(tmp_dir, "sketches")

    if not os.path.isdir(mash_folder):
        os.makedirs(mash_folder, exist_ok=True)

    mash_filepaths = list()

    for filepath in filepaths:
        # Assume these files are all unzipped, in fasta format, and with .fna extension
        mash_filepath = os.path.join(mash_folder, "{}.msh".format(os.path.splitext(os.path.basename(filepath))[0]))

        if not os.path.isfile(mash_filepath):
            try:
                # Build the MASH sketch
                subprocess.check_call(
                    [
                        "mash",
                        "sketch",
                        "-k",
                        str(kmer_size),
                        "-s",
                        str(sketch_size),
                        "-o",
                        filepath,
                        mash_filepath,
                    ],
                    stderr=subprocess.DEVNULL
                )

                if not os.path.isfile(mash_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mash_filepath)

            except subprocess.CalledProcessError as e:
                raise Exception("An error has occurred while running MASH sketch on {}".format(filepath)).with_traceback(e.__traceback__)

        mash_filepaths.append(mash_filepath)

    return mash_filepaths


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
    input_filepaths = [line.strip() for line in open(args.filepath).readlines() if line.strip()]

    # Unzip input files
    print("Unzipping file if Gzip compressed")
    input_filepaths = unzip_files(input_filepaths, args.tmp_dir)

    # Compute MASH sketches
    print("Computing MASH sketches")
    mash_sketches = mash_sketch(input_filepaths, args.tmp_dir, kmer_size=args.kmer_size, sketch_size=args.sketch_size)

    # Paste all the MASH sketches
    print("Pasting MASH sketches")
    mash_pastes = mash_paste(mash_sketches, args.tmp_dir)

    # Compute the MASH distances pair-wise
    print("Computing MASH distances")
    mash_distances = mash_distance(mash_sketches, mash_pastes, args.tmp_dir)

    # Build a condensed matrix with distances
    print("Computing a condensed distance matrix")
    input_names, condensed_matrix = condense_matrix(mash_distances)

    # Run fastcluster
    print("Computing a average-linkage hierarchical clustering")
    dendro = fastcluster.linkage(condensed_matrix, method="average")

    # Define clusters with threshold
    print("Defining clusters with threshold {}".format(args.threshold))
    clusters = hier.fcluster(dendro, args.threshold, criterion="distance")

    # Write assignments
    with open(args.out_file, "w+") as outfile:
        outfile.write("# {}\t{}\n".format("Sequence", "Cluster"))

        for name, cluster in zip(input_names, clusters):
            outfile.write("{}\t{}\n".format(name, cluster))


if __name__ == "__main__":
    main()
