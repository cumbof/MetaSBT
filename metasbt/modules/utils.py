"""Utility functions.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.7"
__date__ = "Sep 24, 2024"

import argparse as ap
import errno
import gzip
import logging
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import tempfile
from ast import literal_eval
from collections.abc import Callable
from logging import Logger
from logging.config import dictConfig
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, TextIO, Tuple, Union

import fastcluster  # type: ignore
import numpy as np  # type: ignore
import scipy.cluster.hierarchy as hier  # type: ignore
import tqdm  # type: ignore
from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore

# Define the list of dependencies
# This is never used but helps to keep track of the external
# software dependencies required by the functions implemented here
DEPENDENCIES = [
    "checkm2",
    "checkv",
    "eukcc",
    "howdesbt",
    "kitsune",
    "ntcard",
    "wget",
]

# Define the list of supported extensions for compressed files
# .rrr compression is used by HowDeSBT only
# Everything else can be Gzip compressed only
COMPRESSED_FILES = [
    ".gz",
    ".rrr",
]

# Define the list of supported extensions for uncompressed files
UNCOMPRESSED_FILES = [
    ".fa",
    ".fna",
    ".fasta",
    ".bf",
    ".txt",
    ".tsv",
]


def bfaction(
    genomes: List[os.path.abspath],
    tmpdir: os.path.abspath,
    kmer_len: int,
    min_occurrences: int=2,
    filter_size: Optional[int]=None,
    nproc: int=1,
    action: str="bfdistance",
    mode: str="theta",
) -> Union[Dict[str, Dict[str, float]], Dict[str, int]]:
    """Wrapper around HowDeSBT bfdistance and bfoperate sub-commands.

    Parameters
    ----------
    genomes : list
        List with paths to the genome or bloom filter files.
    tmpdir : os.path.abspath
        Path to the temporary folder with bloom filters.
    kmer_len : int
        The kmers length.
    min_occurrences : int, default 2
        Exclude kmers with a number of occurrences lower than this parameter.
    filter_size : int, optional
        The bloom filter size.
    nproc : int, default 1
        Make it parallel.
    action : {"bfoperate", "bfdistance"}, default "bfdistance"
        The HowDeSBT sub-command.
    mode : {"and", "or", "xor", "eq", "not", "hamming", "intersect", "union", "theta"}, default "theta"
        The sub-command mode.
        Available values when `action` is "bfoperate" are "and", "or", "xor", "eq", and "not".
        Available values when `action` is "bfdistance" are "hamming", "intersect", "union", and "theta".

    Raises
    ------
    FileNotFoundError
        If a path to a genome file does not exist.
    ValueError
        - If the provided `action` is not supported;
        - If the provided `mode` is not supported based on the provided `action`.
    Exception
        If the number of input genomes is <=2.

    Returns
    -------
    dict
        A dictionary with the results of the bfoperate or bfdistance commands indexed by genome name.

    Notes
    -----
    Please refer to the HowDeSBT documentation for additional information about the bfoperate and bfdistance sub-commands:
    https://github.com/medvedevgroup/HowDeSBT
    """

    # Define supported actions
    actions = ["bfoperate", "bfdistance"]

    if action not in actions:
        raise ValueError('Unsupported action "{}"!'.format(action))

    mode = mode.lower()

    # Define supported modes
    bfoperate_modes = ["and", "or", "xor", "eq", "not"]
    bfdistance_modes = ["hamming", "intersect", "union", "theta"]

    if (action == "bfoperate" and mode not in bfoperate_modes) or (action == "bfdistance" and mode not in bfdistance_modes):
        raise ValueError('Unsupported mode "{}" for action "{}"!'.format(mode, action))

    # Check whether the input genomes exist
    for filepath in genomes:
        if not os.path.isfile(filepath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

    # Keep going in case of 2 or more input genomes
    if len(genomes) < 2:
        raise Exception("The number of input genomes must be >2!")

    # Check whether the temporary folder exists, otherwise create it
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir, exist_ok=True)

    if not filter_size:
        # Estimate the bloom filter size with ntCard
        filter_size = estimate_bf_size(
            genomes,
            kmer_len=kmer_len,
            min_occurrences=min_occurrences,
            prefix="genomes",
            tmp_dir=tmpdir,
            nproc=nproc
        )

    # Take track of all the bloom filter file paths
    bf_files = list()

    howdesbt_log_filepath = os.path.join(tmpdir, "howdesbt.log")
    howdesbt_log = open(howdesbt_log_filepath, "w+")

    # Build the bloom filters
    for genome_path in genomes:
        # Retrieve genome file info
        _, genome_name, extension, compression = get_file_info(genome_path)

        # Define the uncompressed genome path
        genome_file = os.path.join(tmpdir, "{}{}".format(genome_name, extension))

        if not os.path.exists(genome_file):
            if not compression:
                try:
                    # Make a symbolic link in case of an uncompressed file
                    os.symlink(genome_path, genome_file)

                except FileExistsError:
                    # This may occur in a multiprocessing scenario, where different processes
                    # call this function on the same `genome_path`
                    pass

            else:
                # Uncompress the genome file
                with open(genome_file, "w+") as file:
                    run(["gzip", "-dc", genome_path], stdout=file, stderr=file)

        # Define the bloom filter file path
        bf_filepath = os.path.join(tmpdir, "{}.bf".format(genome_name))

        if not os.path.exists(bf_filepath):
            # Build the bloom filter representation of the genome
            run(
                [
                    "howdesbt",
                    "makebf",
                    "--k={}".format(kmer_len),
                    "--min={}".format(min_occurrences),
                    "--bits={}".format(filter_size),
                    "--hashes=1",
                    "--seed=0,0",
                    genome_file,
                    "--out={}".format(bf_filepath),
                    "--threads={}".format(nproc),
                ],
                stdout=howdesbt_log,
                stderr=howdesbt_log,
            )

        if os.path.isfile(bf_filepath):
            bf_files.append(bf_filepath)

    dist = dict()

    if action == "bfdistance":
        bfdistance_tmpdir = os.path.join(tmpdir, "distances")

        if not os.path.isdir(bfdistance_tmpdir):
            os.makedirs(bfdistance_tmpdir, exist_ok=True)

        processes = nproc if nproc <= len(bf_files) else len(bf_files)

        if processes == 1:
            for filepath in bf_files:
                source, target_dists = bfdistance(filepath, bf_files, mode, bfdistance_tmpdir)

                dist[source] = target_dists

        else:
            with mp.Pool(processes=processes) as pool:
                jobs = [pool.apply_async(bfdistance, args=(filepath, bf_files, mode, bfdistance_tmpdir,)) for filepath in bf_files]

                for job in jobs:
                    source, target_dists = job.get()

                    dist[source] = target_dists

    elif action == "bfoperate":
        with tempfile.NamedTemporaryFile() as bflist, tempfile.NamedTemporaryFile() as bfaction_out:
            # Dump the list of bloom filter file paths
            with open(bflist.name, "wt") as bflist_file:
                for filepath in bf_files:
                    bflist_file.write("{}\n".format(filepath))

            with open(bfaction_out.name, "wt") as bfaction_out_file:
                run(
                    [
                        "howdesbt",
                        "bfoperate",
                        "--list={}".format(bflist.name),
                        "--noout",
                        "--{}".format(mode),
                        "--report:counts",
                    ],
                    stdout=bfaction_out_file,
                    stderr=howdesbt_log,
                )

                # Retrieve the output of howdesbt bfoperate
                with open(bfaction_out.name) as bfaction_out_file:
                    for line in bfaction_out_file:
                        line = line.strip()
                        if line:
                            # TODO: this may break in case the absolute paths to the bloom filters contain one or more spaces;
                            #       need to find a better solution
                            line_split = line.split(" ")

                            key = "result"

                            if line_split[0] != key:
                                # Get genome name
                                # Set check_exist=True is not required, but it makes sure that the absolute file path to the
                                # bloom filters has not been altered because of the previous .split(" ")
                                _, key, _, _ = get_file_info(line_split[0], check_supported=False, check_exists=False)

                            # Get active bits
                            dist[key] = int(line_split[-3])

    # Close the log
    howdesbt_log.close()

    return dist


def bfdistance(
    source: os.path.abspath,
    targets: List[os.path.abspath],
    mode: str,
    tmpdir: os.path.abspath,
) -> Tuple[str, Dict[str, float]]:
    """Wrapper around `howdesbt bfdistance` for computing the distance between
    a bloom filter file and a set of bloom filters.

    Parameters
    ----------
    source : os.path.abspath
        Path to the input bloom filter file.
    targets : list
        List with paths to a series of bloom filter files.
    mode : str
        The `bfdistance` mode.
    tmpdir : os.path.abspath
        Path to the temporary foder with the output of `howdesbt bfdistance`.

    Raises
    ------
    FileNotFoundError
        If the input source file does not exist.
    Exception
        If there are not target files in the input list or the output distance file already exists.

    Returns
    -------
    tuple
        A tuple with the source input file path and a dictionary with the distances to the target files.
    """

    if not os.path.isfile(source):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), source)

    if not targets:
        raise Exception("There should be at least one target file!")

    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir, exist_ok=True)

    distance_filepath = os.path.join(tmpdir, "{}__{}.txt".format(os.path.splitext(os.path.basename(source))[0], mode))

    with open(distance_filepath, "a+") as distance_file:
        for target in targets:
            run(
                [
                    "howdesbt",
                    "bfdistance",
                    source,
                    target,
                    "--show:{}".format(mode),
                ],
                stdout=distance_file,
                stderr=subprocess.DEVNULL,
            )

    dist = dict()

    with open(distance_filepath) as distance_file:
        for line in distance_file:
            line = line.strip()
            if line:
                # This is require to replace consecutive space instances with a single space
                # TODO: this may break in case the absolute paths to the bloom filters contain one or more spaces;
                #       need to find a better solution
                line_split = " ".join(line.split()).split(" ")

                # Get genome names
                # Set check_exist=True is not required, but it makes sure that the absolute file path to the
                # bloom filters has not been altered because of the previous .split(" ")
                _, target_name, _, _ = get_file_info(line_split[1].split(":")[0], check_supported=False, check_exists=True)

                # Remove non informative fields
                if line_split[-1] == "({})".format("intersection" if mode == "intersect" else mode):
                    line_split = " ".join(line_split[:-1]).strip().split(" ")

                # Get distance
                dist[target_name] = float(line_split[-1])

    _, source_name, _, _ = get_file_info(source, check_supported=False, check_exists=False)

    return source_name, dist


def build_sh(argv: List[str], module: str, outfolder: os.path.abspath) -> None:
    """Build a sh script with the command line used to launch a module.

    Parameters
    ----------
    argv : list
        List of arguments.
    module : str
        The module name.
    outfolder : os.path.abspath
        Path to the output folder.
    """

    with open(os.path.join(outfolder, "{}.sh".format(module)), "w+") as sh:
        sh.write("#!/bin/bash\n\n")

        # Add metasbt
        argv.insert(0, "metasbt")

        # Replace the path to the python script with the module ID
        argv[1] = module

        # Finally build the command line
        sh.write("{}\n".format(" ".join([os.path.abspath(v) if os.path.exists(v) else v for v in argv])))


def checkm2(
    genomes_paths: Optional[List[os.path.abspath]]=None,
    tmp_dir: Optional[os.path.abspath]=None,
    file_extension: str="fna.gz",
    nproc: int=1,
) -> Dict[str, Dict[str, str]]:
    """Run CheckM2 on a set of genomes. Organize genomes in chunks with 1000 genomes at most.

    Parameters
    ----------
    genomes_paths : list, optional
        List of paths to the input genomes.
    tmp_dir : os.path.abspath, optional
        Path to the temporary folder.
    file_extension : str, default "fna.gz"
        Assume all genomes have the same file extension.
    nproc : int, default 1
        Make the execution of CheckM2 parallel.

    Returns
    -------
    dict
        A dictionary with the CheckM2 stats indexed by the genome names.

    Notes
    -----
    Please refer to the official documentation for additional information about CheckM2:
    https://github.com/chklovski/CheckM2
    """

    # Define the output dictionary
    output_dict = dict()

    # Check whether there is at least one genome path in list
    if genomes_paths:
        run_tmp_dir = os.path.join(tmp_dir, "tmp")

        # Organize genomes
        counter = 0
        run_id = 1

        os.makedirs(os.path.join(run_tmp_dir, "bins_{}".format(run_id)), exist_ok=True)
        os.makedirs(os.path.join(run_tmp_dir, "tmp_{}".format(run_id)), exist_ok=True)

        # Iterate over the list of paths to the genome files
        for genome_path in genomes_paths:
            # Reorganize genomes in chunks with 1000 genomes at most
            if counter % 1000 > 0:
                counter = 0
                run_id += 1
                os.makedirs(os.path.join(run_tmp_dir, "bins_{}".format(run_id)), exist_ok=True)
                os.makedirs(os.path.join(run_tmp_dir, "tmp_{}".format(run_id)), exist_ok=True)

            # Symlink genome files to the bins folder of the current chunk
            os.symlink(
                genome_path,
                os.path.join(run_tmp_dir, "bins_{}".format(run_id), os.path.basename(genome_path)),
            )

        # Iterate over the genomes chunk folders
        for bins_folder in Path(run_tmp_dir).glob("bins_*"):
            if os.path.isdir(str(bins_folder)):
                # Retrieve the run ID from the file path
                run_id = int(os.path.splitext(os.path.basename(str(bins_folder)))[0].split("_")[-1])

                # Create the run folder
                run_dir = os.path.join(tmp_dir, "run_{}".format(run_id))
                os.makedirs(run_dir, exist_ok=True)

                # Define the output table path for the current run
                table_path = os.path.join(tmp_dir, "run_{}.tsv".format(run_id))

                try:
                    with open(table_path, "w+") as table:
                        # Run CheckM2
                        run(
                            [
                                "checkm2",
                                "predict",
                                "--input",
                                bins_folder,
                                "--output-directory",
                                run_dir,
                                "--specific",
                                "--extension",
                                file_extension,
                                "--tmpdir",
                                os.path.join(run_tmp_dir, "tmp_{}".format(run_id)),
                                "--threads",
                                nproc,
                                "--stdout",
                            ],
                            stdout=table,
                            stderr=subprocess.DEVNULL,
                        )

                    if os.path.isfile(table_path):
                        with open(table_path) as table:
                            header = table.readline().strip()

                            for line in table:
                                line = line.strip()

                                if line:
                                    line_split = line.split("\t")

                                    genome = line_split[header.index("Name")]

                                    if "." in file_extension:
                                        # Fix the genome name
                                        genome = ".".join(genome.split(".")[:-1])

                                    output_dict[genome] = dict()

                                    for idx, value in enumerate(line_split):
                                        output_dict[genome][header[idx]] = value

                except Exception:
                    pass

    return output_dict


def checkv(
    genomes_paths: Optional[List[os.path.abspath]]=None,
    tmp_dir: Optional[os.path.abspath]=None,
    file_extension: str="fna.gz",
    nproc: int=1,
) -> Dict[str, Dict[str, str]]:
    """Run CheckV on a set of genomes.

    Parameters
    ----------
    genomes_paths : list, optional
        List of paths to the input genomes.
    tmp_dir : os.path.abspath, optional
        Path to the temporary folder.
    file_extension : str, default "fna.gz"
        Assume all genomes have the same file extension.
    nproc : int, default 1
        Make the execution of CheckV parallel.

    Returns
    -------
    dict
        A dictionary with the CheckV stats indexed by the genome names.

    Notes
    -----
    Please refer to the official documentation for additional information about CheckV:
    https://bitbucket.org/berkeleylab/checkv
    """

    # Define the output dictionary
    output_dict = dict()

    # Check whether there is at least one genome path in list
    if genomes_paths:
        # Merge input files
        merged_filepath = os.path.join(tmp_dir, "merged.fna")

        # Mapping contig ID to file name
        contig_to_filename = dict()

        with open(merged_filepath, "w+") as outfile:
            for genome_path in genomes_paths:
                # This should not make any memory issue since input files are pretty small here
                if file_extension.lower().endswith(".gz"):
                    with gzip.open(genome_path, "rt") as infile:
                        # First line is the contig ID
                        contig_id = infile.readline().strip()[1:]

                        # Take track of the mapping between the contig IDs and file names
                        contig_to_filename[contig_id] = get_file_info(genome_path, check_supported=False, check_exists=False)[1]

                        outfile.write(infile.read())

                else:
                    with open(genome_path) as infile:
                        # First line is the contig ID
                        contig_id = infile.readline().strip()[1:]

                        # Take track of the mapping between the contig IDs and file names
                        contig_to_filename[contig_id] = get_file_info(genome_path)[1]

                        outfile.write(infile.read())

        # Run CheckV
        run(
            [
                "checkv",
                "end_to_end",
                merged_filepath,
                tmp_dir,
                "-t",
                nproc,
            ],
            silence=True,
        )

        output_filepath = os.path.join(tmp_dir, "quality_summary.tsv")

        if os.path.isfile(output_filepath):
            # Read the content of the output table and add a new "Name" column 
            # with the input file names base on the mapping contig_to_filename
            with open(os.path.join(tmp_dir, "run_1.tsv")) as run_table:
                with open(output_filepath) as output_table:
                    header = output_table.readline().strip().split("\t")

                    # Add a new column "Name" with the name of the input files
                    header.insert("Name", 0)

                    # Write the new header line
                    run_table.write("{}\n".format("\t".join(header)))

                    for line in output_table:
                        line = line.strip()
                        if line:
                            line_split = line.split("\t")

                            # Retrieve the contig ID
                            contig_id = line_split[header.index("contig_id")+1]

                            # Add the file name
                            line_split.insert(contig_to_filename[contig_id], 0)

                            # Write the new line
                            run_table.write("{}\n".format("\t".join(line_split)))

                            # Update the output dictionary with quality stats
                            output_dict[contig_to_filename[contig_id]] = dict()

                            for idx, value in enumerate(line_split):
                                output_dict[contig_to_filename[contig_id]][header[idx]] = value

    return output_dict


def ani_distance(
    bf_filepath: os.path.abspath,
    kmer_size: int,
    target_bf_filepaths: List[os.path.abspath],
    tmp_dir: os.path.abspath
) -> Tuple[os.path.abspath, List[float]]:
    """Compute the ANI distances between one sketch and a list of sketches.

    Parameters
    ----------
    bf_filepath : os.path.abspath
        Path to the source BF.
    kmer_size : int
        The length of the k-mers.
    target_bf_filepaths : list
        List with paths to the target BFs.
    tmp_dir : os.path.abspath
        Path to the temporary folder.

    Returns
    -------
    tuple
        A tuple with the path to the source BF file and a list with ANI distances.
    """

    distances = list()

    if target_bf_filepaths:
        # Compute the number of kmer in common between the sketches
        _, intersects = bfdistance(bf_filepath, target_bf_filepaths, "intersect", tmp_dir)

        # Also compute the number of total kmers between sketches
        _, unions = bfdistance(bf_filepath, target_bf_filepaths, "union", tmp_dir)

        for bf_filepath_2 in target_bf_filepaths:
            bf_filepath_2_name = os.path.splitext(os.path.basename(bf_filepath_2))[0]

            # Compute the Jaccard similarity
            jaccard_similarity = round(intersects[bf_filepath_2_name]/unions[bf_filepath_2_name], 5)

            if jaccard_similarity == 0.0:
                # This is to avoid ValueError: math domain error
                jaccard_similarity = sys.float_info.min

            # https://doi.org/10.1093/bioinformatics/btae452 Formula #8
            # Same distance measure implemented in MASH
            ani_distance = 1 - (1 + (1/kmer_size) * math.log((2*jaccard_similarity) / (1+jaccard_similarity)))

            distances.append(ani_distance)

    return bf_filepath, distances


def cluster(
    genomes_list: List[os.path.abspath],
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    manifest_filepath: os.path.abspath,
    profiles_dir: os.path.abspath,
    tmpdir: os.path.abspath,
    outpath: os.path.abspath,
    cluster_prefix: str="MSBT",
    min_occurrences: int=2,
    nproc: int=1,
) -> Tuple[Dict[str, List[os.path.abspath]], List[os.path.abspath]]:
    """Define new clusters with unassigned genomes.

    Parameters
    ----------
    genomes_list : list
        List with paths to the unassigned genomes.
    boundaries : dict
        Dictionary with boundaries as produced by the boundaries module.
    manifest_filepath : os.path.abspath
        Path to the manifest file.
    profiles_dir : oa.path.abspath
        Path to the temporary folder with the genomes profiles defined by the profile module.
    tmpdir : os.path.abspath
        Path to the temporary folder for building the bloom filters.
    outpath : os.path.abspath
        Path to the output file with the new assignments.
    cluster_prefix : str, default "MSBT"
        Prefix of clusters.
    min_occurrences : int, default 2
        Exclude kmers with a number of occurrences lower than this parameter.
    nproc : int, default 1
        Make the HowDeSBT bfdistance sub-command parallel

    Raises
    ------
    FileExistsError
        If the output file with assignments already exists.
    FileNotFoundError
        - If the manifest file does not exist;
        - If the temporary folder with the genomes profiles does not exist.

    Returns
    -------
    tuple
        A tuple with the assignments as a dictionary mapping the taxonomic label with the genome file paths 
        in addition to the list of paths to the unassigned genomes.
    """

    # Check whether the output file already exists
    if os.path.isfile(outpath):
        raise FileExistsError(errno.ENOENT, os.strerror(errno.ENOENT), outpath)

    # Also check whether the input files already exist
    # Otherwise, raise an exception
    if not os.path.isfile(manifest_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    if not os.path.isdir(profiles_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), profiles_dir)

    # Load the database manifest file
    manifest = load_manifest(manifest_filepath)

    # Retrieve the kmer length
    kmer_len = manifest["kmer_len"]

    # Also, retrieve the clusters counter
    clusters_counter = manifest["clusters_counter"]

    # Build up a condensed ANI distance matrix
    condensed_distance_matrix = list()

    dists = dict()

    # Define the path to the bloom filters
    bfs = list()

    for filepath in genomes_list:
        genome_name = get_file_info(filepath, check_supported=False, check_exists=False)[1]

        # Filters should all already be under this temporary folder
        bf_path = os.path.join(tmpdir, "{}.bf".format(genome_name))

        bfs.append(bf_path)

    with mp.Pool(processes=nproc) as pool:
        jobs = [
            pool.apply_async(
                ani_distance,
                args=(
                    filepath,
                    kmer_len,
                    bfs[filepath_pos+1:],
                    os.path.join(tmpdir, "distances"),
                )
            ) for filepath_pos, filepath in enumerate(bfs)
        ]

        for job in jobs:
            source_bf, target_dists = job.get()

            dists[source_bf] = target_dists

    for filepath in genomes_list:
        genome_name = get_file_info(filepath, check_supported=False, check_exists=False)[1]

        bf_filepath = os.path.join(tmpdir, "{}.bf".format(genome_name))

        condensed_distance_matrix.extend(dists[bf_filepath])

    dendro = fastcluster.linkage(condensed_distance_matrix, method="average")

    profiles = dict()

    # Define the list of taxonomic levels
    levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    for filepath in genomes_list:
        genome_name = get_file_info(filepath, check_supported=False, check_exists=False)[1]

        # Retrieve the genome profile
        profile = os.path.join(profiles_dir, "{}__profiles.tsv".format(genome_name))

        # Check whether the profile exists
        if os.path.isfile(profile):
            profiles[filepath] = dict()

            with open(profile) as file:
                for line in file:
                    line = line.strip()

                    if line:
                        if not line.startswith("#"):
                            line_split = line.split("\t")

                            if line_split[1] in levels:
                                profiles[filepath][line_split[1]] = {
                                    "taxonomy": line_split[2],
                                    "score": float(line_split[4])
                                }

                            if line_split[1] == "genome":
                                profiles[filepath][line_split[1]] = {
                                    "taxonomy": "|".join(line_split[2].split("|")[:-1]),
                                    "score": float(line_split[4])
                                }

    assignments = dict()

    unassigned = list()

    for level_pos, level in enumerate(reversed(levels)):
        # Start from the species level up to the kingdom
        processed = set()

        for filepath_pos, filepath in enumerate(genomes_list):
            if filepath in profiles and filepath not in processed:
                # Retrieve the closest cluster under the current taxonomic level
                level_taxonomy = profiles[filepath][level]["taxonomy"]

                level_taxonomy_split = level_taxonomy.split("|")

                # Retrieve the closeness score
                level_score = profiles[filepath][level]["score"]

                # Retrive the closest cluster boundaries
                min_score, _, level_centroid = get_level_boundaries(boundaries, level_taxonomy)

                # `level_centroid` could be `None` in case it does not contain enough clusters/genomes to compute its boundaries
                # In that case we can still use the closest cluster/genome
                if level_centroid is not None and level_taxonomy_split[-1][3:] != level_centroid:
                    # Retrieve the input genome name
                    genome_name = get_file_info(filepath, check_supported=False, check_exists=False)[1]

                    # Retrieve the database root folder from the manifest file path
                    db_dir = os.path.dirname(manifest_filepath)

                    if level == "species":
                        centroid_bf = os.path.join(db_dir, level_taxonomy.replace("|", os.sep), "filters", "{}.bf".format(level_centroid))

                    else:
                        centroid_bf = os.path.join(db_dir, level_taxonomy.replace("|", os.sep), level_centroid, "{}.bf".format(level_centroid))

                    # Compute the distance between the input genome and the centroid
                    _, distances = ani_distance(
                        os.path.join(tmpdir, "{}.bf".format(genome_name)),
                        kmer_len,
                        [centroid_bf],
                        tmpdir,
                    )

                    # `distances` contains a single value here
                    # Convert the distance back to a similarity measure
                    # Replace `level_score` with the similarity between the input genome and the centroid of the closest cluster
                    level_score = abs(distances[0] - 1)

                if level_score >= min_score:
                    # The current genome must be assigned to the closest cluster
                    # We should cut the dendrogram using the closest cluster boundaries to check for other assignments
                    # `min_score` is a similarity, thus, it must be converted to a distance measure
                    clusters = hier.fcluster(dendro, 1 - min_score, criterion="distance")

                    # Search for the cluster id assigned to the current genome
                    cluster_id = clusters[filepath_pos]

                    # Collect all the genomes that have been assigned to the same cluster
                    for filepath2, assigned_cluster in zip(genomes_list, clusters):
                        if assigned_cluster == cluster_id:
                            # All genomes that fall under this cluster are assigned to an already defined cluster in MetaSBT
                            assignments[filepath2] = level_taxonomy

                            # Mark the assigned genomes as processed
                            processed.add(filepath2)

            if level == "kingdom" and filepath not in processed:
                # This genome has not been characterized, not even at the kingdom level
                # We should simply report these genomes, we are not going to define new kingdoms
                unassigned.append(filepath)

    # Get full and partial assignments first
    characterized = dict()

    partially_characterized = dict()

    for filepath in assignments:
        if "|s__" not in assignments[filepath]:
            if assignments[filepath] not in partially_characterized:
                partially_characterized[assignments[filepath]] = set()

            # These genomes have not been fully characterized and must be further processed
            partially_characterized[assignments[filepath]].add(filepath)

        else:
            if assignments[filepath] not in characterized:
                characterized[assignments[filepath]] = set()

            characterized[assignments[filepath]].add(filepath)

    # We can now finish characterizing the partially characterized genomes
    # Genomes in `partially_characterized` are all characterized at the kingdom level at least
    while len(partially_characterized) > 0:
        for taxonomic_assignment in list(partially_characterized.keys()):
            # We previously defined the assignments starting from the species level up to the kingdom
            # We should now start from the kingdom level down to the species level to finish characterizing genomes
            # If a genome has not been characterized, it means that it is not close enough to its closest cluster according to its boundaries
            # We can now create new clusters and estimate their boundaries, then cut the dendrogram to see how many genomes fall in them
            tmp_level_taxonomy = "{}|{}__tmp".format(taxonomic_assignment, levels[taxonomic_assignment.count("|")+1][0])

            # Estimate the boundaries for the temporary cluster
            # This is a totally new cluster. There is no need to search for a centroid here
            min_score, _, _ = get_level_boundaries(boundaries, tmp_level_taxonomy)

            # Cluster genomes according to the temporary cluster's boundaries
            clusters = hier.fcluster(dendro, 1 - min_score, criterion="distance")

            # Keep track of the mapping between fcluster ids and MetaSBT cluster ids
            clusters_map = dict()

            for filepath, assigned_cluster in zip(genomes_list, clusters):
                # Consider genomes under the current partially characterized taxonomy only
                if filepath in partially_characterized[taxonomic_assignment]:
                    if assigned_cluster not in clusters_map:
                        # Create a new MSBT cluster
                        # Increment the cluster counter
                        clusters_counter += 1

                        clusters_map[assigned_cluster] = clusters_counter

                    # Define the assigned taxonomic label
                    assigned_taxonomy = "{}|{}__{}{}".format(
                        taxonomic_assignment,
                        levels[taxonomic_assignment.count("|")+1][0],
                        cluster_prefix,
                        clusters_map[assigned_cluster]
                    )

                    if "|s__" not in assigned_taxonomy:
                        # This is not fully characterized yet
                        if assigned_taxonomy not in partially_characterized:
                            partially_characterized[assigned_taxonomy] = set()

                        partially_characterized[assigned_taxonomy].add(filepath)

                    else:
                        # This is fully characterized now
                        if assigned_taxonomy not in characterized:
                            characterized[assigned_taxonomy] = set()

                        characterized[assigned_taxonomy].add(filepath)

            # The partial taxonomy must be removed from the bucket of partially characterized labels
            del partially_characterized[taxonomic_assignment]

    # Update the manifest with the new clusters counter
    if clusters_counter > manifest["clusters_counter"]:
        # Load the manifest file
        with open(manifest_filepath) as manifest_file:
            manifest_lines = manifest_file.readlines()

        # Update the --clusters-counter info
        with open(manifest_filepath, "w+") as manifest_file:
            for line in manifest_lines:
                line = line.strip()

                if line:
                    line_split = line.split(" ")

                    if line_split[0] == "--clusters-counter":
                        line_split[-1] = str(clusters_counter)

                    manifest_file.write("{}\n".format(" ".join(line_split)))

    # Dump the new assignments to the output file
    with open(outpath, "w+") as out:
        # Add header line
        out.write("# Genome\tAssignment\tCluster ID\n")

        for taxonomy in sorted(characterized.keys()):
            for genome_path in sorted(characterized[taxonomy]):
                # Get genome name
                _, genome, _, _ = get_file_info(genome_path)

                cluster_id = taxonomy.split("|")[-1][3:]

                out.write("{}\t{}\t{}\n".format(genome, taxonomy, cluster_id))

    return characterized, unassigned


def dereplicate_genomes(
    genomes: List[os.path.abspath],
    tax_id: str,
    tmp_dir: os.path.abspath,
    kmer_len: int,
    min_occurrences: int=2,
    filter_size: Optional[int]=None,
    nproc: int=1,
    similarity: float=0.99,
) -> List[os.path.abspath]:
    """Dereplicate genomes based on ANI similarity.

    Parameters
    ----------
    genomes : list
        List with paths to the genome files.
    tax_id : str
        The NCBI tax ID.
    tmp_dir : os.path.abspath
        Path to the temporary folder.
    kmer_len : int
        The length of the kmers.
    min_occurrences : int, default 2
        Minimum number of kmer occurrences.
    filter_size : int, optional
        The size of the bloom filters.
    nproc : int, default 1
        Make it parallel.
    similarity : float, default 0.99
        The similarity threshold on the ANI similarity.

    Returns
    -------
    list
        A list with paths to the genome files that passed the dereplication process.
    """

    # Define the temporary folders
    os.makedirs(tmp_dir, exist_ok=True)

    # Compute the intersection of kmers between all the input genomes
    bfdistance_intersect = bfaction(
        genomes,
        tmp_dir,
        kmer_len,
        min_occurrences=min_occurrences,
        filter_size=filter_size,
        nproc=nproc,
        action="bfdistance",
        mode="intersect"
    )

    # Compute the union of kmers between all the input genomes
    bfdistance_union = bfaction(
        genomes,
        tmp_dir,
        kmer_len,
        min_occurrences=min_occurrences,
        filter_size=filter_size,
        nproc=nproc,
        action="bfdistance",
        mode="union"
    )

    filtered_genomes_filepath = os.path.join(tmp_dir, "filtered.txt")
 
    filtered_genomes = set()

    # Pair-wise comparison of input genomes
    for i, genome1_path in enumerate(genomes):
        if genome1_path not in filtered_genomes:
            # Get genome1 file name
            _, genome1, _, _ = get_file_info(genome1_path)

            for genome2_path in genomes[i+1:]:
                # Get genome2 file name
                _, genome2, _, _ = get_file_info(genome2_path)

                # Compute the Jaccard similarity
                jaccard_similarity = bfdistance_intersect[genome1][genome2]/bfdistance_union[genome1][genome2]

                if jaccard_similarity == 0.0:
                    # This is to avoid ValueError: math domain error
                    jaccard_similarity = sys.float_info.min

                # https://doi.org/10.1093/bioinformatics/btae452 Formula #8
                # Same similarity measure implemented in MASH
                ani_similarity = round(1 + (1/kmer_len) * math.log((2*jaccard_similarity) / (1+jaccard_similarity)), 5)

                if ani_similarity >= similarity:
                    # Keep track of the excluded genomes in the filtered file
                    with open(filtered_genomes_filepath, "a+") as f:
                        f.write("{}\n".format(genome2_path))

                    filtered_genomes.add(genome2_path)

    # Redefine the list of genomes by removing the filtered ones
    genomes = list(set(genomes).difference(set(filtered_genomes)))

    return genomes


def download(
    url: Optional[str]=None,
    urls: Optional[List[str]]=None,
    folder: os.path.abspath=os.getcwd(),
    retries: int=10,
    raise_exception: bool=True
) -> Optional[Union[os.path.abspath, List[os.path.abspath]]]:
    """Download a file from URL to the specified folder.

    Parameters
    ----------
    url : str, optional
        The URL to the source file.
    urls : list, optional
        A list of URLs to the source files.
    folder : os.path.abspath, default os.getcwd
        Path to the output folder.
    retries : int, default 10
        Try downloading again in case of errors.
    raise_exception : bool, default True
        Raise an exception in case of errors if True.

    Raises
    ------
    ValueError
        If no URLs are provided.
    Exception
        If an error occurs while downloading the files.
    FileNotFoundError
        If pass the download process but the file does not exist.

    Returns
    -------
    A path or a list of paths to the downloaded files.

    Examples
    --------
    >>> import os
    >>> from utils import download
    >>> url = "https://raw.githubusercontent.com/cumbof/MetaSBT-DBs/main/databases.tsv"
    >>> filepath = download(url=url)
    >>> os.path.isfile(filepath)
    True

    Download the databases table from the MetaSBT-DBs repository and check whether the 
    file exists on the file system. By default, it is saved into the home directory.
    """

    if not url and not urls:
        raise ValueError("No URLs provided")

    # Check whether the destination folder path exists
    if not os.path.isdir(folder):
        os.makedirs(folder, exist_ok=True)

    try:
        if url:
            # Download file from URL to the destination folder
            run(
                ["wget", "-N", url, "-P", folder],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                retries=retries,
            )

            filepath = os.path.join(folder, url.split(os.sep)[-1])

            if not os.path.isfile(filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

            return filepath

        elif urls:
            filepaths = list()

            with tempfile.NamedTemporaryFile() as tmpfile:
                # Dump the list of bloom filter file paths
                with open(tmpfile.name, "wt") as tmpfile_list:
                    for url in urls:
                        tmpfile_list.write("{}\n".format(url))

                        filepaths.append(os.path.join(folder, url.split(os.sep)[-1]))

                # Download a list of files from URL
                run(
                    ["wget", "-N", "-i", tmpfile.name, "-P", folder],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    retries=retries,
                )

            for filepath in filepaths:
                if not os.path.isfile(filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

            return filepaths

    except Exception as e:
        if raise_exception:
            raise Exception(
                "An error has occurred while trying to download {}".format(url)
            ).with_traceback(e.__traceback__)

        # This file does not seem really important after all
        return None


def estimate_bf_size(
    genomes: List[os.path.abspath],
    kmer_len: int=21,
    min_occurrences: int=2,
    prefix: str="genomes",
    tmp_dir: os.path.abspath=os.getcwd(),
    nproc: int=1
) -> int:
    """Estimate the bloom filter size with ntCard.

    Parameters
    ----------
    genomes : list
        List of paths to the genome files.
    kmer_len : int, default 21
        Length of the kmers.
    min_occurrences : int, default 2
        Exclude kmers with a number of occurrences lower than this parameter.
    prefix : str, default "genomes"
        Prefix of the output histogram file.
    tmp_dir : os.path.abspath, default os.getcwd
        Path to the temporary folder.
    nproc : int, default 1
        Make it parallel.

    Raises
    ------
    Exception
        If it is not able to estimate the bloom filter size.

    Returns
    -------
    int
        The estimated bloom filter size.
    """

    os.makedirs(tmp_dir, exist_ok=True)

    with tempfile.NamedTemporaryFile() as genomes_file:
        # Dump the list of genome file paths
        with open(genomes_file.name, "wt") as gfile:
            for filepath in genomes:
                gfile.write("{}\n".format(filepath))

        # Estimate the bloom filter size with ntCard
        run(
            [
                "ntcard",
                "--kmer={}".format(kmer_len),
                "--threads={}".format(nproc),
                "--pref={}".format(os.path.join(tmp_dir, prefix)),
                "@{}".format(genomes_file.name),
            ],
            silence=True,
        )

    # Total number of kmers in the reads
    F1 = 0

    # Number of distinct kmers
    F0 = 0

    # List with the number of kmers occurring less than min_occurrences
    fs = list()

    # Read the ntCard output hist file
    hist_filepath = os.path.join(tmp_dir, "{}_k{}.hist".format(prefix, kmer_len))

    with open(hist_filepath) as histfile:
        for line in histfile:
            line = line.strip()
            if line:
                line_split = line.split()

                if line_split[0] == "F1":
                    F1 = int(line_split[-1])

                elif line_split[0] == "F0":
                    F0 = int(line_split[-1])

                elif isinstance(line_split[0], int):
                    if int(line_split[0]) < min_occurrences:
                        fs.append(int(line_split[-1]))

                    else:
                        break

    if F0 == 0:
        # This could happen in case of a single very small genome
        # Use F1 as the bloom filter size in this case
        if F1 == 0:
            raise Exception("Unable to estimate the bloom filter size: {}".format(hist_filepath))

        return F1

    # Estimate the bloom filter size
    return F0 - sum(fs)


def eukcc(
    genomes_paths: Optional[List[os.path.abspath]]=None,
    tmp_dir: Optional[os.path.abspath]=None,
    file_extension: str="fna.gz",
    nproc: int=1,
) -> Dict[str, Dict[str, str]]:
    """Run EukCC on a set of genomes.

    Parameters
    ----------
    genomes_paths : list, optional
        List of paths to the input genomes.
    tmp_dir : os.path.abspath, optional
        Path to the temporary folder.
    file_extension : str, default "fna.gz"
        Assume all genomes have the same file extension.
    nproc : int, default 1
        Make the execution of EukCC parallel

    Returns
    -------
    dict
        A dictionary with the EukCC stats indexed by the genome names.

    Notes
    -----
    Please refer to the official documentation for additional information about EukCC:
    https://github.com/Finn-Lab/EukCC
    """

    # Define the output dictionary
    output_dict = dict()

    # Check whether there is at least one genome path in list
    if genomes_paths:
        bins_folder = os.path.join(tmp_dir, "bins")

        os.makedirs(bins_folder, exist_ok=True)

        for genome_path in genomes_paths:
            os.symlink(genome_path, os.path.join(bins_folder, os.path.basename(genome_path)))

        # TODO
        # Use file_extension to change the default file type as argument to EukCC
        # Does it work with Gzip compressed sequences?

        # Run EukCC
        run(
            [
                "eukcc",
                "folder",
                "--out",
                tmp_dir,
                "--threads",
                nproc,
                bins_folder,
            ],
            silence=True,
        )

        output_filepath = os.path.join(tmp_dir, "eukcc.tsv")

        if os.path.isfile(output_filepath):
            # Read the content of the output table and add a new "Name" column 
            with open(os.path.join(tmp_dir, "run_1.tsv")) as run_table:
                with open(output_filepath) as output_table:
                    header = output_table.readline().strip().split("\t")

                    # Add a new column "Name" with the name of the input files
                    header.insert("Name", 0)

                    # Write the new header line
                    run_table.write("{}\n".format("\t".join(header)))

                    for line in output_table:
                        line = line.strip()
                        if line:
                            line_split = line.split("\t")

                            # Retrieve the file name
                            filename = get_file_info(line_split[header.index("fasta")+1], check_supported=False, check_exists=False)[1]

                            # Add the file name
                            line_split.insert(filename, 0)

                            # Write the new line
                            run_table.write("{}\n".format("\t".join(line_split)))

                            # Update the output dictionary with quality stats
                            output_dict[filename] = dict()

                            for idx, value in enumerate(line_split):
                                output_dict[filename][header[idx]] = value

    return output_dict


def filter_quality(
    quality_dict: Dict[str, Dict[str, str]], completeness: float=0.0, contamination: float=100.0
) -> List[str]:
    """Filter genomes according to completeness and contamination criteria.

    Parameters
    ----------
    quality_dict : dict
        Dictionary with quality information indexed by genome names.
    completeness : float, default 0.0
        Minimum allowed completeness.
    contamination : float, default 100.0
        Maximum allowed contamination.

    Returns
    -------
    list
        A list of genomes that passed the quality-control criteria.
    """

    # Define the list of genomes that passed the quality control
    genomes = list()

    for genome in quality_dict:
        # Check whether the current genome respect both the completeness and contamination criteria
        comp = None
        cont = None

        if "Completeness_Specific" in quality_dict[genome]:
            # CheckM2
            comp = float(quality_dict[genome]["Completeness_Specific"])
            cont = float(quality_dict[genome]["Contamination"])

        elif "completeness" in quality_dict[genome]:
            # CheckV and EukCC
            comp = float(quality_dict[genome]["completeness"])
            cont = float(quality_dict[genome]["contamination"])

        if isinstance(comp, float) and isinstance(cont, float):
            if comp >= completeness and cont <= contamination:
                genomes.append(genome)

    return genomes


def get_bf_density(filepath: os.path.abspath) -> float:
    """Retrieve the bloom filter density.

    Parameters
    ----------
    filepath : os.path.abspath
        Path to the bloom filter file.

    Raises
    ------
    Exception
        If an error occurs while retrieving the bloom filter density with HowDeSBT.

    Returns
    -------
    float
        The density of the bloom filter.
    """

    density = 0.0

    with tempfile.NamedTemporaryFile() as dumpbf:
        # Retrieve bloom filter density
        run(
            [
                "howdesbt",
                "dumpbf",
                filepath,
                "--show:density",
            ],
            stdout=dumpbf,
            stderr=dumpbf,
        )

        try:
            # Get the result
            density = float(open(dumpbf.name, "rt").readline().strip().split(" ")[-1])

        except Exception as ex:
            raise Exception("An error has occurred while retrieving bloom filter density:\n{}".format(filepath)).with_traceback(
                ex.__traceback__
            )

    return density


def get_boundaries(
    bfs: List[os.path.abspath],
    tmpdir: os.path.abspath,
    kmer_len: int,
    filter_size: Optional[int]=None,
    min_occurrences: int=2,
    nproc: int=1
) -> Tuple[int, int, str]:
    """Return kmers boundaries for a specific set of genomes defined as the minimum and
    maximum ANI similarities between all the genomes in the current taxonomic level.

    Parameters
    ----------
    bfs : list
        List with paths to the bloom filter representations of the genomes.
    tmpdir : os.path.abspath
        Path to the temporary folder.
    kmer_len : int
        The length of the kmers.
    filter_size : int, optional
        The bloom filter size.
    min_occurrences : int, default 2
        Minimum number of kmer occurrences.
    nproc : int, default 1
        Make it parallel.

    Returns
    -------
    tuple
        A tuple with the minimum and maximum ANI similarities between all the input genomes
        in addition to the centroid defined as the genome that maximize the similarity with
        all the other genomes
    """

    # Define the clusters' boundaries
    minv = np.Inf
    maxv = 0

    # Compute the number of kmers in common between all pairs of genomes
    bfdistance_intersect = bfaction(
        bfs,
        tmpdir,
        kmer_len,
        min_occurrences=min_occurrences,
        filter_size=filter_size,
        nproc=nproc,
        action="bfdistance",
        mode="intersect"
    )

    # Also, compute the total number of kmers as the union of kmers in every pair of genomes
    bfdistance_union = bfaction(
        bfs,
        tmpdir,
        kmer_len,
        min_occurrences=min_occurrences,
        filter_size=filter_size,
        nproc=nproc,
        action="bfdistance",
        mode="union"
    )

    similarities = dict()

    # Pair-wise comparison of input bloom filters
    for i, bf1_path in enumerate(bfs):
        # Get bf1 file name
        _, bf1, _, _ = get_file_info(bf1_path)

        if bf1 not in similarities:
            similarities[bf1] = list()

        for bf2_path in bfs[i+1:]:
            # Get bf2 file name
            _, bf2, _, _ = get_file_info(bf2_path)

            if bf2 not in similarities:
                similarities[bf2] = list()

            try:
                if bfdistance_intersect[bf1][bf2] == 0:
                    # This could be only due to a wrong classification and must be reported
                    with open(os.path.join(tmpdir, "zero_common_kmers.tsv"), "a+") as zck:
                        zck.write("{}\t{}\n".format(bf1, bf2))

                    # Pass to the next comparison
                    continue

                # Compute the Jaccard similarity
                jaccard_similarity = bfdistance_intersect[bf1][bf2]/bfdistance_union[bf1][bf2]

                if jaccard_similarity == 0.0:
                    # This is to avoid ValueError: math domain error
                    jaccard_similarity = sys.float_info.min

                # https://doi.org/10.1093/bioinformatics/btae452 Formula #8
                # Same similarity measure implemented in MASH
                ani_similarity = round(1 + (1/kmer_len) * math.log((2*jaccard_similarity) / (1+jaccard_similarity)), 5)

                similarities[bf1].append(ani_similarity)
                similarities[bf2].append(ani_similarity)

                # Update the minimum and maximum scores
                if ani_similarity > maxv:
                    maxv = ani_similarity

                if ani_similarity < minv:
                    minv = ani_similarity

            except KeyError:
                # Only in case one of the two current genomes is not into the bfdistance_intersect dictionary
                # This should never happen
                pass

    # Search for the centroid as the genome that maximize the similarity with all the other genomes
    centroid = None

    average_similarity = 0.0

    for bf in similarities:
        bf_average_similarity = sum(similarities[bf])/len(similarities[bf])

        if bf_average_similarity > average_similarity:
            average_similarity = bf_average_similarity

            centroid = bf

    return minv, maxv, centroid


def get_file_info(filepath: os.path.abspath, check_supported: bool=True, check_exists: bool=True) -> Tuple[str, str, str, str]:
    """Get file path, name, extension, and compression.

    Parameters
    ----------
    filepath : os.path.abspath
        Path to the input file.
    check_supported : bool, default True
        Check if the file extension is supported if True.
    check_exists : bool, default True
        Check if the file exists if True.

    Raises
    ------
    FileNotFoundError
        If `check_exists` is True and the file does not exist.
    Exception
        If `check_supported` is True and the file extension is not recognized.

    Returns
    -------
    tuple
        A tuple with the file's folder path, name, extension, and compression.
    """

    if check_exists and not os.path.isfile(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

    # Trim the folder path out
    basename = os.path.basename(filepath)

    # Split the basename
    filename, extension = os.path.splitext(basename)

    # Take track of the extension in case of compression
    compression = None

    # Check whether it is compressed
    if extension in COMPRESSED_FILES:
        compression = extension
        filename, extension = os.path.splitext(filename)

    # Check whether the input file is supported
    if check_supported and extension not in UNCOMPRESSED_FILES:
        raise Exception("Unrecognized input file")

    # Retrieve the absolute path to the file folder
    absdir = os.path.abspath(os.path.dirname(filepath))

    return absdir, filename, extension, compression


def get_level_boundaries(boundaries: Dict[str, Dict[str, Union[int, float]]], taxonomy: str) -> Tuple[int, int]:
    """Retrieve boundaries for a given taxonomic label.

    Parameters
    ----------
    boundaries : dict
        Dictionary with boundaries as produced by the boundaries module.
    taxonomy : str
        The taxonomic label.

    Raises
    ------
    Exception
        If it is not able to retrieve the boundaries.

    Returns
    -------
    tuple
        A tuple with the minimum and maximum scores as the specified cluster's boundaries, plus the centroid.
        The centroid could be `None` in case of temporary boundaries.
    """

    if taxonomy in boundaries:
        return boundaries[taxonomy]["min_ani"], boundaries[taxonomy]["max_ani"], boundaries[taxonomy]["centroid"]

    taxonomy_levels = taxonomy.split("|")[:-1]

    # Keep track of the min and max 
    min_bounds = list()
    max_bounds = list()

    factor = 1

    while not min_bounds and not max_bounds:
        for tax in boundaries:
            tax_split = tax.split("|")

            if len(tax_split) == len(taxonomy_levels) + factor and "|".join(tax_split[:-factor]) == "|".join(taxonomy_levels):
                min_bounds.append(boundaries[tax]["min_ani"])
                max_bounds.append(boundaries[tax]["max_ani"])

        if not min_bounds and not max_bounds and len(taxonomy_levels) == 1:
            raise Exception("Unable to retrieve boundaries for {}".format(taxonomy))

        factor += 1

        taxonomy_levels = taxonomy_levels[:-1]

    min_bound = round(sum(min_bounds)/len(min_bounds), 5)
    max_bound = round(sum(max_bounds)/len(max_bounds), 5)

    return min_bound, max_bound, None


def howdesbt(
    level_dir: os.path.abspath,
    extension: str="fna.gz",
    kmer_len: int=21,
    min_occurrences: int=2,
    filter_size: int=10000,
    nproc: int=1,
    flat_structure: bool=False,
) -> None:
    """Run HowDeSBT on a specific taxonomic level. Genomes must be in the "genomes" folder under `level_dir`.

    Parameters
    ----------
    level_dir : os.path.abspath
        Path to the taxonomic level folder.
    extension : str, default "fna.gz"
        Input file extension.
    kmer_len : int, default 21
        The length of the kmers.
    min_occurrences : int, default 2
        Exclude kmers with a number of occurrences lower than this parameter.
    filter_size : int, default 10000
        The size of the bloom filters.
    nproc : int, default 1
        Make it parallel.
    flat_structure : bool, default False
        Genomes are not taxonomically organized if False.

    Raises
    ------
    FileNotFoundError
        If the provided `level_dir` does not exist.

    Notes
    -----
    Please refer to the official documentation for additional information about HowDeSBT:
    https://github.com/medvedevgroup/HowDeSBT
    """

    # Check whether the input folder exists
    if not os.path.isdir(level_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), level_dir)

    # Extract the level name from the level folder path
    level_name = os.path.basename(level_dir)

    # Define the index folder
    index_dir = os.path.join(level_dir, "index")
    if os.path.isdir(index_dir):
        # Remove old index folder if any
        shutil.rmtree(index_dir, ignore_errors=True)

    # Define the path to the file with the list of genome under the current taxonomic level
    level_list = os.path.join(level_dir, "{}.txt".format(level_name))
    if os.path.isfile(level_list):
        os.unlink(level_list)

    # Define the path to the bloom filter representation of the current taxonomic level
    level_filter = os.path.join(level_dir, "{}.bf".format(level_name))
    if os.path.isfile(level_filter):
        os.unlink(level_filter)

    # Define the log file
    howdesbt_log_filepath = os.path.join(level_dir, "howdesbt.log")
    howdesbt_log = open(howdesbt_log_filepath, "w+")

    # Take track of how many genomes under the specific taxonomic levels
    how_many = 0

    if os.path.basename(level_dir).startswith("s__") or flat_structure:
        # Search for all the genomes under the current taxonomic level
        genomes_folder = os.path.join(level_dir, "genomes")
        if os.path.isdir(genomes_folder):
            # Create the filters folder
            filters_dir = os.path.join(os.path.dirname(genomes_folder), "filters")
            os.makedirs(filters_dir, exist_ok=True)

            # Iterate over the genome files
            for genome_path in Path(genomes_folder).glob("*.{}".format(extension)):
                # Retrieve genome file info
                _, genome_name, genome_extension, genome_compression = get_file_info(genome_path)

                # Define the path to the bloom filter representation of the genome
                bf_filepath = os.path.join(filters_dir, "{}.bf".format(genome_name))

                if not os.path.isfile(bf_filepath):
                    # Define the uncompressed genome path
                    genome_file = os.path.join(genomes_folder, "{}{}".format(genome_name, genome_extension))

                    if genome_compression:
                        # Uncompress the genome file
                        with open(genome_file, "w+") as file:
                            run(["gzip", "-dc", genome_path], stdout=file, stderr=file)

                    # Build the bloom filter file from the current genome
                    run(
                        [
                            "howdesbt",
                            "makebf",
                            "--k={}".format(kmer_len),
                            "--min={}".format(min_occurrences),
                            "--bits={}".format(filter_size),
                            "--hashes=1",
                            "--seed=0,0",
                            genome_file,
                            "--out={}".format(bf_filepath),
                            "--threads={}".format(nproc),
                        ],
                        stdout=howdesbt_log,
                        stderr=howdesbt_log,
                    )

                    if genome_compression:
                        # Get rid of the uncompressed genome file
                        os.unlink(genome_file)

        filters_folder = os.path.join(level_dir, "filters")

        if os.path.isdir(filters_folder):
            # Take track of the bloom filter files
            with open(level_list, "w+") as level_list_file:
                for bf_filepath in Path(filters_folder).glob("*.bf"):
                    level_list_file.write("{}\n".format(bf_filepath))

                    # Increment the genomes counter
                    how_many += 1

    else:
        # Find all the other taxonomic levels
        # Define the new list of bloom filters
        for level in os.listdir(level_dir):
            if os.path.isdir(os.path.join(level_dir, level)):
                # Defile the path to the bloom filter file
                bf_filepath = os.path.join(level_dir, level, "{}.bf".format(level))

                if os.path.isfile(bf_filepath):
                    with open(level_list, "a+") as level_list_file:
                        level_list_file.write("{}\n".format(bf_filepath))

                    # Increment the genomes counter
                    how_many += 1

    # Build the index folder
    os.makedirs(index_dir, exist_ok=True)

    # Move to the index folder
    # This will force howdesbt to build the compressed nodes into the index folder
    os.chdir(index_dir)

    if how_many > 1:
        # Create the tree topology file
        run(
            [
                "howdesbt",
                "cluster",
                "--list={}".format(level_list),
                "--bits={}".format(filter_size),
                "--tree={}".format(os.path.join(index_dir, "union.sbt")),
                "--nodename={}".format(os.path.join(index_dir, "node{number}")),
                "--keepallnodes",
            ],
            stdout=howdesbt_log,
            stderr=howdesbt_log,
        )

    else:
        # With only one bloom filter it does not make sense to cluster genomes
        bf_filepath = [line.strip() for line in open(level_list).readlines() if line.strip()][0]

        # There is only one line which refers to the only bloom filter file
        shutil.copy(bf_filepath, os.path.join(level_dir, "{}.bf".format(level_name)))

        # Manually define the union.sbt file with the single node
        with open(os.path.join(index_dir, "union.sbt"), "w+") as union:
            union.write("{}\n".format(bf_filepath))

    # Build all the bloom filter files
    run(
        [
            "howdesbt",
            "build",
            "--howde",
            "--tree={}".format(os.path.join(index_dir, "union.sbt")),
            "--outtree={}".format(os.path.join(index_dir, "index.detbrief.sbt")),
        ],
        stdout=howdesbt_log,
        stderr=howdesbt_log,
    )

    # Remove the union.sbt file
    os.unlink(os.path.join(index_dir, "union.sbt"))

    # Fix node paths in the final index.detbrief.sbt file
    with open(os.path.join(index_dir, "index.full.detbrief.sbt"), "w+") as file1:
        with open(os.path.join(index_dir, "index.detbrief.sbt")) as file2:
            for line in file2:
                line = line.strip()
                if line:
                    # Define the depth of the node in the tree
                    stars = line.count("*")
                    # Remove the stars to retrieve the node name
                    node_name = line[stars:]
                    # Define the absolute path to the node bloom filter file
                    node_path = os.path.join(index_dir, node_name)
                    # Define the new node in the tree
                    file1.write("{}{}\n".format("*" * stars, node_path))

    # Get rid of the old tree
    os.unlink(os.path.join(index_dir, "index.detbrief.sbt"))

    # Rename the new tree
    shutil.move(
        os.path.join(index_dir, "index.full.detbrief.sbt"),
        os.path.join(index_dir, "index.detbrief.sbt"),
    )

    if how_many > 1:
        # Build the bloom filter representation of the current taxonomic level
        bf_filepath = os.path.join(level_dir, "{}.bf".format(level_name))

        # Merge all the leaves together by applying the OR logic operator on the bloom filter files
        # The resulting bloom filter is the representative one, which is the same as the root node of the tree
        run(
            [
                "howdesbt",
                "bfoperate",
                "--list={}".format(level_list),
                "--or",
                "--out={}".format(bf_filepath),
            ],
            stdout=howdesbt_log,
            stderr=howdesbt_log,
        )

    # Close the log file handler
    howdesbt_log.close()


def init_logger(filepath: Optional[os.path.abspath]=None, toolid: Optional[str]=None, verbose: bool=True) -> Optional[Logger]:
    """Define a logger to print on console, on file, or both.

    Parameters
    ----------
    filepath : os.path.abspath, optional
        Path to the log file.
    toolid : str, optional
        The ID of the module.
    verbose : bool, default True
        Print messages on the stdout.

    Raises
    ------
    FileNotFoundError
        If the folder of the provided `filepath` does not exist.

    Returns
    -------
    logging.Logger
        The logging.Logger object.
    """

    # Define the logger config
    # TODO configure other logging levels (i.e., NOTSET, DEBUG, INFO, WARN, ERROR, and CRITICAL)
    logging_config: Dict[str, Any] = dict(
        version=1,
        formatters={
            "verbose": {
                "format": "[%(toolid)s][%(levelname)s][%(asctime)s] %(message)s",
                "datefmt": "%d/%b/%Y %H:%M:%S",
            }
        },
        handlers={
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "verbose",
                "stream": sys.stdout,
            },
            "file": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "INFO",
                "formatter": "verbose",
                "filename": os.devnull,
                "maxBytes": 52428800,
                "backupCount": 7,
            },
        },
        loggers={
            "console": {"handlers": ["console"], "level": logging.INFO},
            "file": {"handlers": ["file"], "level": logging.INFO},
            "full": {"handlers": ["console", "file"], "level": logging.INFO},
        },
    )

    # In case of log file
    if filepath:
        # Check whether its folder exists
        log_dir = os.path.dirname(filepath)
        if not os.path.isdir(log_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), log_dir)

        # Update the log file path in the config dictionary
        logging_config["handlers"]["file"]["filename"] = filepath

    # Load the logging config
    dictConfig(logging_config)

    # Get the record factory
    factory = logging.getLogRecordFactory()

    # Customise the record_factory function to add the toolid attribute
    def record_factory(*args, **kwargs):
        record = factory(*args, **kwargs)
        record.toolid = toolid
        return record

    # Register the new record factory
    logging.setLogRecordFactory(record_factory)

    # Define the logger type
    logtype = None

    if filepath and verbose:
        # Full logger will print on screen and on the log file
        logtype = "full"

    elif filepath and not verbose:
        # File logger will only print on the log file
        logtype = "file"

    elif not filepath and verbose:
        # Console logger will only print message on the screen
        logtype = "console"

    if logtype:
        # Define and return the logger object
        logger = logging.getLogger(logtype)
        return logger

    # In case no file path and verbose have been specified
    return None


def integrity_check(filepath: os.path.abspath) -> bool:
    """This is for Gzipped files only.

    Parameters
    ----------
    filepath : os.path.abspath
        Path to the Gzipped file.

    Raises
    ------
    Exception
        If the input file does not have a ".gz" extension.

    Returns
    -------
    bool
        True if the input file passes the integrity check.
    """

    # This checks whether the input file exists and its extension and compression are supported
    _, _, _, compression = get_file_info(filepath, check_supported=True, check_exists=True)

    if compression != ".gz":
        # Limit the compression to Gzipped files only
        raise Exception("Unsupported file type")

    try:
        # It always throws an Exception in case of a return value > 0
        run(["gzip", "-t", filepath], silence=True)

    except Exception:
        return False

    return True


def load_boundaries(boundaries_filepath: os.path.abspath) -> Dict[str, Dict[str, Union[int, float]]]:
    """Load the table produced by the boundaries module.

    Parameters
    ----------
    boundaries_filepath : os.path.abspath
        Path to the boundaries table.

    Returns
    -------
    dict
        A dictionary with the content of the boundaries table indexed by taxa.
    """

    if not os.path.isfile(boundaries_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), boundaries_filepath)

    boundaries = dict()

    with open(boundaries_filepath) as table:
        for line in table:
            line = line.strip()
            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")

                    # Indexed by taxonomic labels
                    boundaries[line_split[0]] = {
                        "nodes": int(line_split[1]),
                        "min_ani": float(line_split[2]),
                        "max_ani": float(line_split[3]),
                        "centroid": line_split[4],
                    }

    return boundaries


def load_input_table(filepath: os.path.abspath, input_extension: str="fna.gz") -> Dict[str, List[os.path.abspath]]:
    """Load the input table with the list of paths to the genome files and eventually their taxonomic labels.

    Parameters
    ----------
    filepath : os.path.abspath
        Path to the input file.
    input_extension : str, default "fna.gz"
        The extension of the genome files listed in the input table.

    Raises
    ------
    FileNotFoundError
        If the provided `filepath` does not exist.
    Exception
        - If the input table contains more than 2 columns;
        - If a taxonomic label contains more or less than 7 taxonomic levels;
        - If a genome file has a different extension compared to the provided one;
        - If a genome is reported twice.

    Returns
    -------
    dict
        A dictionary with the list of genome paths indexed by their taxonomic labels.
    """

    if not os.path.isfile(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

    genome2taxonomy: Dict[str, str] = dict()

    with open(filepath) as input_file:
        for line in input_file:
            line = line.strip()

            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")

                    taxonomy = "NA"

                    if len(line_split) > 2:
                        raise Exception("Malformed input file! It must contain two columns at most")

                    elif len(line_split) == 2:
                        taxonomy = line_split[1]

                        if len(taxonomy.split("|")) != 7:
                            # Taxonomic labels must have 7 levels
                            raise Exception(
                                "Invalid taxonomic label! Please note that taxonomies must have 7 levels:\n{}".format(
                                    line_split[1]
                                )
                            )

                    # This automatically check whether extension and compression are supported
                    dirpath, genome_name, extension, compression = get_file_info(line_split[0])

                    if not line_split[0].endswith(".{}".format(input_extension)):
                        raise Exception(
                            "Unexpected input file extension! "
                            "File: {}; Expected extension: {}".format(line_split[0], input_extension)
                        )

                    genome_path = os.path.join(dirpath, "{}{}{}".format(
                        genome_name, extension, compression if compression else ""
                    ))

                    if genome_path in genome2taxonomy:
                        if genome2taxonomy[genome_path] != taxonomy:
                            raise Exception(
                                "Genome \"{}\" appears twice in the input file with two different taxonomic labels:\n{}\n{}".format(
                                    genome_name, genome2taxonomy[genome_path], taxonomy
                                )
                            )

                    genome2taxonomy[genome_path] = taxonomy

    taxonomy2genomes: Dict[str, List[str]] = dict()

    if genome2taxonomy:
        for genome_path in genome2taxonomy:
            if genome2taxonomy[genome_path] not in taxonomy2genomes:
                taxonomy2genomes[genome2taxonomy[genome_path]] = list()

            taxonomy2genomes[genome2taxonomy[genome_path]].append(genome_path)

    return taxonomy2genomes


def load_manifest(manifest_filepath: os.path.abspath) -> Dict[str, Union[str, int, float]]:
    """Load the manifest file.

    Parameters
    ----------
    manifest_filepath : os.path.abspath
        Path to the manifest file.

    Raises
    ------
    FileNotFoundError
        If the input manifest file does not exist.

    Returns
    -------
    dict
        A dictionary with manifest data.
    """

    if not os.path.isfile(manifest_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    manifest = dict()

    with open(manifest_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                line_split = line.split(" ")

                # e.g., key: --kmer-len > kmer_len
                key = line_split[0][2:].replace("-", "_")

                try:
                    # Try to cast values to the appropriate type
                    manifest[key] = literal_eval(line_split[1])

                except Exception:
                    # Otherwise, maintain value as string
                    manifest[key] = line_split[1]

    return manifest


def number(
    typev: type,
    minv: Optional[Union[int, float]]=None,
    maxv: Optional[Union[int, float]]=None,
) -> Callable:
    """Take full control of input numeric types by defining custom intervals.

    Parameters
    ----------
    typev : type
        A numerical type, int or float.
    minv : {int, float}, optional
        The lower bound of the numerical interval.
    maxv : {int, float}, optional
        The upper bound of the numerical interval.

    Raises
    ------
    argparse.ArgumentTypeError
        - If the type of the provided input is not numerical.
        - If the provided number falls outside the interval.

    Returns
    -------
    collections.abc.Collable
        The collections.abc.Callable object.
    """

    def type_func(value: Union[int, float]) -> Union[int, float]:
        """
        Test data type and ranges on the input value
        """

        try:
            value = typev(value)

            if minv and value < minv:
                raise ap.ArgumentTypeError("Minimum value is {}".format(minv))

            if maxv and value > maxv:
                raise ap.ArgumentTypeError("Maximum value is {}".format(maxv))

            return value

        except Exception as e:
            raise ap.ArgumentTypeError("Input value must be {}".format(typev)).with_traceback(e.__traceback__)

    return type_func


def optimal_k(
    genomes: List[os.path.abspath],
    kl: int,
    tmpdir: os.path.abspath,
    closely_related: bool=False,
    nproc: int=1,
    threads: int=1
) -> int:
    """Given a set of genomes, try to define the best k-mer length with kitsune.

    Parameters
    ----------
    genomes : list
        List of paths to the genome files (Gzip compressed or not).
    kl : int
        kitsune tests different kmer lengths, starting from k=4 up to this parameter.
    tmpdir : os.path.abspath
        Path to the temporary folder.
    closely_related : bool, default False
        For closely related genomes, set this parameter to True.
    nproc : int, default 1
        Maximum number of processes.
    threads : int, default 1
        Maximum number of threads.

    Raises
    ------
    ValueError
        if the provided `kl` is lower than 4.
    Exception
        - If the size of `genomes` is <2;
        - If the actual number of genomes is <2;
        - If an error occurs while running kitsune kopt.

    Returns
    -------
    int
        The optimal kmer length.

    Notes
    -----
    Please refer to the official documentation for additional information about kitsune:
    https://github.com/natapol/kitsune
    """

    if len(genomes) < 2:
        raise Exception("Not enough genomes")

    if kl < 4:
        raise ValueError("Initial k-mer length is too small")

    # Check whether the destination folder path exists
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir, exist_ok=True)

    # Take track of the genome file paths in the tmp folder
    genomes_paths = list()

    for genome_path in genomes:
        _, genome_name, extension, compression = get_file_info(genome_path, check_supported=True, check_exists=True)

        # Define the uncompressed genome path
        genome_file = os.path.join(tmpdir, "{}{}".format(genome_name, extension))

        if not compression:
            # Make a symbolic link in case of an uncompressed file
            os.symlink(genome_path, genome_file)

        else:
            # Uncompress the genome file
            # It can always be Gzip compressed here
            with open(genome_file, "w+") as file:
                run(["gzip", "-dc", genome_path], stdout=file, stderr=file)

        genomes_paths.append(genome_file)

    if len(genomes_paths) < 2:
        raise Exception("Not enough genomes. Something went wrong while processing your genomes")

    with tempfile.NamedTemporaryFile() as inputlist, tempfile.NamedTemporaryFile() as outres:
        # Dump the list of bloom filter file paths
        with open(inputlist.name, "wt") as inputlist_file:
            for filepath in genomes_paths:
                inputlist_file.write("{}\n".format(filepath))

        # Run kitsune
        # This may take a while and a considerable amount of computational resources
        run(
            [
                "kitsune",
                "kopt",
                "--filenames",
                inputlist.name,
                "--k-max",
                kl,
                "--canonical",
                "--fast",
                "--nproc",
                nproc,
                "--threads",
                threads,
                "--in-memory",
                "--output",
                outres.name,
                "--closely_related" if closely_related else ""
            ],
            silence=True
        )

        # Get kitsune output message
        out_content = open(outres.name).read().strip()

        try:
            # Try to retrieve the optimal k
            return int(out_content.split(" ")[-1])

        except Exception as ex:
            raise Exception("An error has occurred while running kitsune kopt:\n{}".format(out_content)).with_traceback(
                ex.__traceback__
            )


def println(message: str, logger: Optional[Logger]=None, verbose: bool=True) -> None:
    """Send messages to the logger. It will print messages on screen, send messages to the log file, or both.

    Parameters
    ----------
    message : str
        The custom message.
    logger : logging.Logger, optional
        The logging.Logger object.
    verbose : bool, default True
        Print messages on the stdout if True and `logger` is None.
    """

    if logger:
        # Redirect messages to the logger
        logger.info(message)

    elif verbose:
        # In case the logger is not defined
        # Redirect messages to the screen
        print(message)


def quality(
    method: str="CheckM2",
    args: Optional[Dict[str, Any]]=None,
) -> Dict[str, Dict[str, str]]:
    """Wrapper around CheckM2, CheckV, and EukCC for assessing the quality of the input genomes based on their kingdom.

    Parameters
    ----------
    method : {"CheckM2", "CheckV", "EukCC"}, default "CheckM2"
        The quality-control method.
    args : dict, optional
        A dictionary with method-specific arguments.

    Raises
    ------
    ValueError
        If the provided `method` is not supported.

    Returns
    -------
    dict
        A dictionary with the quality stats indexed by genome names.
    """

    # Define the set of supported methods
    supported_methods = {
        "checkm2": checkm2,
        "checkv": checkv,
        "eukcc": eukcc
    }

    if method.lower() not in supported_methods:
        raise ValueError("Unsupported method \"{}\"".format(method))

    # Unpack the args dictionary and pass arguments to the selected function
    return supported_methods[method.lower()](**args)


def run(
    cmdline: List[Union[str, int, float]],
    stdout: Union[int, TextIO]=sys.stdout,
    stderr: Union[int, TextIO]=sys.stderr,
    silence: bool=False,
    extended_error: bool=False,
    retries: int=1,
) -> None:
    """Wrapper for the subprocess.check_call function.

    Parameters
    ----------
    cmdline : list
        The command line list of arguments and values.
    stdout : TextIO, default sys.stdout
        The standard output.
    stderr : TextIO, default sys.stderr
        The standard error.
    silence : bool, default False
        Redirect the stdout and stderr to /dev/null.
    extended_error : bool, default False
        Raise errors with traceback in case of unexpected exceptions.
    retries : int, default 1
        Try running the process again in case of errors.

    Raises
    ------
    Exception
        - If the provided command line is empty;
        - If the provided command line failed to run.
    """

    # Check whether ther is something to run
    if cmdline:
        while retries > 0:
            try:
                # Cast everything to string in cmdline
                cmdline = [str(cmd) for cmd in cmdline]

                # In case of silence
                if silence:
                    # Redirect the stdout and stderr to /dev/null
                    stdout = subprocess.DEVNULL
                    stderr = subprocess.DEVNULL

                # Run a specific command line and redirect the stdout and stderr
                # to those specified in input
                subprocess.check_call(cmdline, stdout=stdout, stderr=stderr)

                # At this point, the execution of the command did not raise any exception
                # Set retries to 0
                retries = 0

            except subprocess.CalledProcessError as e:
                if retries == 1:
                    # Define the error message
                    error_message = "\nAn error has occurred while running the following command:\n{}\n\n".format(
                        " ".join(cmdline)
                    )

                    if extended_error:
                        # Extend the error message
                        error_message += (
                            "If you think this is a bug and need support, please open an Issue or a new Discussion on the official GitHub repository.\n"
                            "We would be happy to answer your questions and help you troubleshoot any kind of issue with our framework.\n"
                        )

                    raise Exception(error_message).with_traceback(e.__traceback__)

                # Try again
                retries -= 1

    else:
        # There is nothing to run
        raise Exception("Empty command line!")


def strand(filepath: os.path.abspath, tmpdir: os.path.abspath) -> os.path.abspath:
    """Restructure a fasta file evaluating its records and their reverse-complement. 
    Sequences are evaluated by translating them and searching for ORFs.
    The best one is selected according to the highest number of ORFs.

    Parameters
    ----------
    filepath : os.path.abspath
        Path to the input sequence file (can be Gzip compressed).
    tmpdir : os.path.abspath
        Path to the temporary folder used to uncompress the input file (if Gzip compressed) 
        and store the new sequence file with the evaluated records.

    Returns
    -------
    str
        Path to the new fasta file in the temporary folder (always uncompressed).
    """

    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir, exist_ok=True)

    # Call get_file_info to check whether the file exists and it is supported
    # Also check if it is compressed
    _, filename, extension, compression = get_file_info(filepath, check_supported=True, check_exists=True)

    if compression:
        uncompressed_filepath = os.path.join(tmpdir, "{}{}~".format(filename))

        # Uncompress the input file
        # This file will be removed at the end of this function
        with open(uncompressed_filepath, "w+") as file:
            run(["gzip", "-dc", filepath], stdout=file, stderr=file)

        filepath = uncompressed_filepath

    # A regular expression to search for ORFs
    regex_orf = re.compile(r'M[^*]{25,}?\*')

    # Keep track of the fasta file records
    records = list()

    # Iterate over all the record in the input fasta file
    # Input files are always in fasta format
    for record in SeqIO.parse(filepath, "fasta"):
        sequence = record.seq.upper()

        if not re.match("^[ACTUGN]*$", str(sequence)):
            # Not a valid sequence
            # Go ahead with the next record
            continue

        # Replace Us with Ts
        sequence = sequence.replace("U", "T")

        positive_strand = sequence
        longest_CDS = 0

        # Search for the positive strand among the original sequence
        # and its reverse complement
        strands = [sequence, sequence.reverse_complement()]

        for strand in strands:
            for frame in range(3):
                protein_sequence = ""

                for fragment in range(frame, len(strand), 3):
                    codon = strand[fragment:fragment+3]

                    if len(codon) == 3:
                        # Translate to aminoacid
                        # Return a Bio.Seq object
                        protein_sequence += codon.translate()

                # Search for ORFs
                matches = regex_orf.findall(str(protein_sequence))
                all_ORFs = "".join([match for match in matches if match])

                if len(all_ORFs)/float(len(strand)) > longest_CDS:
                    longest_CDS = len(all_ORFs)/float(len(strand))
                    positive_strand = strand

        strand_record = SeqRecord(positive_strand, id=record.id, name=record.name,
                                  description=record.description)

        records.append(strand_record)

    if compression:
        # Get rid of the uncompressed file
        os.unlink(filepath)

    # Output sequence file
    out_filepath = os.path.join(tmpdir, "{}{}".format(filename, extension))

    with open(out_filepath, "w+") as strand_file:
        for record in records:
            strand_file.write(strand_record.format("fasta"))

    return out_filepath


def validate_url(url: str) -> bool:
    """Validate a URL.

    Parameters
    ----------
    url : str
        The input URL that must be validated.

    Returns
    -------
    bool
        True if the input URL is valid, False otherwise.
    """

    regex = re.compile(
        r'^(?:http|ftp)s?://'  # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain
        r'localhost|'  # localhost
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # or ip
        r'(?::\d+)?'  # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE
    )

    return re.match(regex, url) is not None
