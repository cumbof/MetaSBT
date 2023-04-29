"""
Utility functions
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 28, 2023"

import argparse as ap
import errno
import logging
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
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union

import numpy as np  # type: ignore

# Define the list of dependencies
# This is never used but helps to keep track of the external
# software dependencies required by the functions implemented here
DEPENDENCIES = [
    "checkm",
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
    genomes: List[str],
    tmpdir: str,
    kmer_len: int,
    min_occurrences: int = 2,
    filter_size: Optional[int] = None,
    nproc: int = 1,
    action: str = "bfdistance",
    mode: str = "theta",
) -> Union[Dict[str, Dict[str, float]], Dict[str, int]]:
    """
    bfdistance and bfoperate wrapper

    :param genomes:             List with paths to the genome or bloom filter files
    :param tmpdir:              Path to the temporary folder with bloom filters
    :param kmer_len:            Kmer length
    :param min_occurrences:     Exclude kmers with a number of occurrences less than this param
    :param filter_size:         Bloom filter size
    :param nproc:               Make it parallel
    :param action:              "bfoperate" or "bfdistance"
    :param mode:                bfoperate modes: "and", "or", "xor", "eq", and "not"
                                bfdistance modes: "hamming", "intersect", "union", and "theta"
    :return:                    Dictionary with the result of bfdistance or bfoperate
    """

    # Define supported actions
    actions = ["bfoperate", "bfdistance"]

    if action not in actions:
        raise Exception('Unsupported action "{}"!'.format(action))

    mode = mode.lower()

    # Define supported modes
    bfoperate_modes = ["and", "or", "xor", "eq", "not"]
    bfdistance_modes = ["hamming", "intersect", "union", "theta"]

    if (action == "bfoperate" and mode not in bfoperate_modes) or (
        action == "bfdistance" and mode not in bfdistance_modes
    ):
        raise Exception('Unsupported mode "{}" for action "{}"!'.format(mode, action))

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
                # Make a symbolic link in case of an uncompressed file
                os.symlink(genome_path, genome_file)

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

    with tempfile.NamedTemporaryFile() as bflist, tempfile.NamedTemporaryFile() as bfaction_out:
        # Dump the list of bloom filter file paths
        with open(bflist.name, "wt") as bflist_file:
            for filepath in bf_files:
                bflist_file.write("{}\n".format(filepath))

        with open(bfaction_out.name, "wt") as bfaction_out_file:
            if action == "bfdistance":
                run(
                    [
                        "howdesbt",
                        "bfdistance",
                        "--list={}".format(bflist.name),
                        "--show:{}".format(mode),
                    ],
                    stdout=bfaction_out_file,
                    stderr=howdesbt_log,
                )

                # Retrieve the output of howdesbt bfdistance
                with open(bfaction_out.name) as bfaction_out_file:
                    for line in bfaction_out_file:
                        line = line.strip()
                        if line:
                            # This is require to replace consecutive space instances with a single space
                            line_split = " ".join(line.split()).split(" ")

                            # Get genome names
                            _, genome1, _, _ = get_file_info(line_split[0].split(":")[0], check_exists=False)
                            _, genome2, _, _ = get_file_info(line_split[1].split(":")[0], check_exists=False)

                            # Remove non informative fields
                            if line_split[-1] == "({})".format("intersection" if mode == "intersect" else mode):
                                line_split = " ".join(line_split[:-1]).strip().split(" ")

                            if genome1 not in dist:
                                dist[genome1] = dict()

                            # Get distance
                            dist[genome1][genome2] = float(line_split[-1])

            elif action == "bfoperate":
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
                            line_split = line.split(" ")

                            key = "result"
                            if line_split[0] != key:
                                # Get genome name
                                _, key, _, _ = get_file_info(line_split[0], check_exists=False)

                            # Get active bits
                            dist[key] = int(line_split[-3])

    # Close the log
    howdesbt_log.close()

    return dist


def build_sh(argv: List[str], module: str, outfolder: str) -> None:
    """
    Build a sh script with the command line used to launch a module

    :param argv:        List of arguments
    :param module:      Module ID
    :param outfolder:   Output folder path
    """

    with open(os.path.join(outfolder, "{}.sh".format(module)), "w+") as sh:
        sh.write("#!/bin/bash\n\n")

        # Add metasbt
        argv.insert(0, "metasbt")

        # Replace the path to the python script with the module ID
        argv[1] = module

        # Finally build the command line
        sh.write("{}\n".format(" ".join([os.path.abspath(v) if os.path.exists(v) else v  for v in argv])))


def checkm(
    genomes_paths: List[str],
    tmp_dir: str,
    file_extension: str = "fna.gz",
    nproc: int = 1,
    pplacer_threads: int = 1,
) -> List[str]:
    """
    Run CheckM on a set of genomes
    Organise genomes in chunks with 1000 genomes at most

    :param genomes_paths:       List of paths to the input genomes
    :param tmp_dir:             Path to the temporary folder
    :param file_extension:      Assume all genomes have the same file extension
    :param nproc:               Make the execution CheckM parallel
    :param pplacer_threads:     Maximum number of threads for pplacer
    :return:                    Return the list of paths to the CheckM output tables
    """

    # Define the output list of paths to the CheckM tables
    output_tables = list()

    # Check whether there is at least one genome path in list
    if genomes_paths:
        run_tmp_dir = os.path.join(tmp_dir, "tmp")

        # Organise genomes
        counter = 0
        run_id = 1
        os.makedirs(os.path.join(run_tmp_dir, "bins_{}".format(run_id)), exist_ok=True)

        # Iterate over the list of paths to the genome files
        for genome_path in genomes_paths:
            # Reorganise genomes in chunks with 1000 genomes at most
            if counter % 1000 > 0:
                counter = 0
                run_id += 1
                os.makedirs(os.path.join(run_tmp_dir, "bins_{}".format(run_id)), exist_ok=True)

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
                    # Run CheckM
                    # TODO update to CheckM2
                    run(
                        [
                            "checkm",
                            "lineage_wf",
                            "-t",
                            nproc,
                            "-x",
                            file_extension,
                            "--pplacer_threads",
                            pplacer_threads,
                            "--tab_table",
                            "-f",
                            table_path,
                            bins_folder,
                            run_dir,
                        ],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )

                    # Add the output table path to the output list
                    output_tables.append(table_path)

                except Exception:
                    pass

    return output_tables


def cluster(
    genomes_list: List[str],
    boundaries: Dict[str, Dict[str, Union[int, float]]],
    manifest_filepath: str,
    profiles_dir: str,
    tmpdir: str,
    outpath: str,
    cluster_prefix: str = "MSBT",
    min_occurrences: int = 2,
    nproc: int = 1,
) -> Dict[str, str]:
    """
    Define new clusters with the unassigned MAGs

    :param genomes_list:        List with paths to the unassigned genomes
    :param boundaries:          Boundaries table produced by the boundaries module
    :param manifest_filepath:   Path to the manifest file
    :param profiles_dir:        Path to the temporary folder with the genomes profiles defined by the profile module
    :param tmpdir:              Path to the temporary folder for building bloom filters
    :param outpath:             Path to the output file with the new assignments
    :param cluster_prefix:      Prefix of clusters numerical identifiers
    :param min_occurrences:     Exclude kmers with a number of occurrences less than this param
    :param nproc:               Make bfdistance parallel
    :return:                    Return the assignments as a dictionary <genome_path, taxonomy>
                                Also return the list of paths to the unassigned genomes
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

    # Retrieve the kmer length, filter size, and clusters counter from the manifest file
    manifest = load_manifest(manifest_filepath)

    kmer_len = manifest["kmer_len"]
    clusters_counter_manifest = manifest["clusters_counter"]

    # Estimate the proper bloom filter size for the set of unassigned genomes
    filter_size = estimate_bf_size(
        genomes_list,
        kmer_len=kmer_len,
        min_occurrences=min_occurrences,
        prefix="genomes",
        tmp_dir=tmpdir,
        nproc=nproc
    )

    # Retrieve the list of input genomes
    genomes = [get_file_info(genome_path)[1] for genome_path in genomes_list]

    # Start counting new clusters
    clusters_counter = clusters_counter_manifest

    # Define the list of taxonomic levels for sorting profiles
    levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    level_ids = [lv[0] for lv in levels]

    # Keep track of the already assigned genomes
    assigned_taxa: Dict[str, List[str]] = dict()
    assigned_genomes: List[str] = list()

    # Keep track of those genomes that MetaSBT is not able to assign
    unassigned: List[str] = list()

    # Compute pair-wise distance between genomes as the number of common kmers
    # This could take a while
    bfdistance_intersect = bfaction(
        genomes_list, tmpdir, kmer_len, filter_size=filter_size, nproc=nproc, action="bfdistance", mode="intersect"
    )

    # Iterate over genomes
    for i in range(len(genomes_list)):
        if genomes[i] not in assigned_genomes:
            # Retrieve the genome profile
            profile = os.path.join(profiles_dir, "{}__profiles.tsv".format(genomes[i]))

            # Check whether the profile exists
            if os.path.isfile(profile):
                # Load levels and scores
                level2match = dict()

                with open(profile) as file:
                    for line in file:
                        line = line.strip()
                        if line:
                            if not line.startswith("#"):
                                line_split = line.split("\t")
                                if line_split[1] in levels:
                                    # Key: taxonomic level
                                    # Value: kmers in common with the taxonomic level
                                    level2match[line_split[1]] = {
                                        "taxonomy": line_split[2],
                                        "common_kmers": int(line_split[3].split("/")[0])
                                    }

                                if line_split[1] == "genome":
                                    # Override the species level with strains info
                                    level2match["species"] = {
                                        "taxonomy": "|".join(line_split[2].split("|")[:-1]),
                                        "common_kmers": int(line_split[3].split("/")[0])
                                    }

                assigned_taxonomy = None

                last_known_level_mink = 0

                # From the species up to the kingdom level
                for level in reversed(levels):
                    # Get level boundaries
                    mink, _ = get_level_boundaries(boundaries, level2match[level]["taxonomy"])
                
                    if level2match[level]["common_kmers"] >= mink and mink > 0:
                        assigned_taxonomy = level2match[level]["taxonomy"]

                        last_known_level_mink = mink

                        break

                if not assigned_taxonomy:
                    # Unable to assign a taxonomic label to the current genome
                    unassigned.append(genomes_list[i])
                
                else:
                    assignment = assigned_taxonomy.split("|")

                    # Fill the assignment with missing levels
                    assigned_levels = len(assignment)

                    for pos in range(assigned_levels, len(levels)):
                        clusters_counter += 1

                        # Create new clusters
                        assignment.append("{}__{}{}".format(level_ids[pos], cluster_prefix, clusters_counter))

                    # Compose the assigned (partial) label
                    assigned_taxonomy = "|".join(assignment)

                    # Assigne current genome to the taxonomy
                    if assigned_taxonomy not in assigned_taxa:
                        assigned_taxa[assigned_taxonomy] = list()
                    
                    assigned_taxa[assigned_taxonomy].append(genomes_list[i])

                    # Mark current genome as assigned
                    assigned_genomes.append(genomes[i])

                    # Check whether other input genomes look pretty close to the current genome by computing
                    # the number of kmers in common between the current genome and all the other input genomes
                    for j in range(i + 1, len(genomes_list)):
                        # Kmers in common have been already computed
                        # It returns a float by default
                        common = int(bfdistance_intersect[genomes[i]][genomes[j]])

                        if common >= last_known_level_mink:
                            # Set the second genome as assigned
                            assigned_genomes.append(genomes[j])
                            # Also assign these genomes to the same taxonomy assigned to the current genome
                            assigned_taxa[assigned_taxonomy].append(genomes_list[j])

    # Update the manifest with the new clusters counter
    if clusters_counter > clusters_counter_manifest:
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

    # Mapping genome -> taxonomy, cluster
    assignment = dict()

    # Dumpt the new assignments to the output file
    with open(outpath, "w+") as out:
        # Add header line
        out.write("# Genome\tAssignment\tCluster ID\n")

        for taxonomy in sorted(assigned_taxa.keys()):
            for genome_path in sorted(assigned_taxa[taxonomy]):
                # Get genome name
                _, genome, _, _ = get_file_info(genome_path)

                cluster_id = taxonomy.split("|")[-1][3:]

                out.write("{}\t{}\t{}\n".format(genome, taxonomy, cluster_id))

                # Take track of mapping genome - taxonomy
                assignment[genome_path] = {
                    "taxonomy": taxonomy,
                    "cluster": cluster_id
                }

    return assignment, unassigned


def dereplicate_genomes(
    genomes: list,
    tax_id: str,
    tmp_dir: str,
    kmer_len: int,
    filter_size: Optional[int] = None,
    nproc: int = 1,
    similarity: float = 1.0,
) -> List[str]:
    """
    Dereplicate genomes

    :param genomes:         List of genome file paths
    :param tax_id:          NCBI tax ID
    :param tmp_dir:         Path to the temporary folder
    :param kmer_len:        Length of the kmers
    :param filter_size:     Size of the bloom filters
    :param nproc:           Make it parallel
    :param similarity:      Similarity threshold on the theta distance
                            Theta between two genomes A and B is defined as N/D, where
                            N is the number of 1s in common between A and B, and
                            D is the number of 1s in A
    :return:                List of genome file paths for genomes that passed the dereplication
    """

    # Define the HowDeSBT temporary folder
    howdesbt_tmp_dir = os.path.join(tmp_dir, "howdesbt", tax_id)
    os.makedirs(howdesbt_tmp_dir, exist_ok=True)

    filtered_genomes_filepath = os.path.join(howdesbt_tmp_dir, "filtered.txt")
    filtered_genomes = list()

    # Compute the theta distance between all the input genomes
    bfdistance_theta = bfaction(
        genomes, howdesbt_tmp_dir, kmer_len, filter_size=filter_size, nproc=nproc, action="bfdistance", mode="theta"
    )

    # Pair-wise comparison of input genomes
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            # Get genome file names
            _, genome1, _, _ = get_file_info(genomes[i])
            _, genome2, _, _ = get_file_info(genomes[j])

            excluded = None

            if bfdistance_theta[genome1][genome2] >= similarity:
                filtered_genomes.append(genomes[i])
                excluded = genomes[i]
            
            if bfdistance_theta[genome2][genome1] >= similarity:
                filtered_genomes.append(genomes[j])
                excluded = genomes[j]

            if excluded:
                # Also take note if the excluded genomes in the filtered file
                with open(filtered_genomes_filepath, "a+") as f:
                    f.write("{}\n".format(excluded))

                break

    # Redefine the list of genomes by removing the filtered ones
    genomes = list(set(genomes).difference(set(filtered_genomes)))

    return genomes


def download(
    url: Optional[str] = None,
    urls: Optional[List[str]] = None,
    folder: str = os.getcwd(),
    retries: int = 10,
    raise_exception: bool = True
) -> Optional[Union[str, List[str]]]:
    """
    Download a file from URL to the specified folder

    :param url:             Source file URL
    :param urls:            List with source file URLs
    :param folder:          Target destination folder path
    :param retries:         Try downloading again in case of errors
    :param raise_exception: Raise an exception in case of error
    :return:                Path or list of paths to the downloaded files
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

        elif urls:
            with tempfile.NamedTemporaryFile() as tmpfile:
                # Dump the list of bloom filter file paths
                with open(tmpfile.name, "wt") as tmpfile_list:
                    for url in urls:
                        tmpfile_list.write("{}\n".format(url))

                # Download a list of files from URL
                run(
                    ["wget", "-N", "-i", tmpfile.name, "-P", folder],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    retries=retries,
                )

    except Exception as e:
        if raise_exception:
            raise Exception(
                "An error has occurred while trying to download {}".format(url)
            ).with_traceback(e.__traceback__)

        # This file does not seem really important after all
        return None

    return os.path.join(folder, url.split(os.sep)[-1])


def estimate_bf_size(
    genomes: str,
    kmer_len: int = 21,
    min_occurrences: int = 2,
    prefix: str = "genomes",
    tmp_dir: str = os.getcwd(),
    nproc: int = 1
) -> int:
    """
    Estimate the bloom filter size with ntCard

    :param genomes:         List of paths to the genome files
    :param kmer_len:        Length of the kmers
    :param min_occurrences: Exclude kmers with a number of occurrences less than this param
    :param prefix:          Prefix of the output histogram file
    :param tmp_dir:         Path to the temporary folder
    :param nproc:           Make it parallel
    :return:                The estimated bloom filter size
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


def filter_checkm_tables(
    checkm_tables: List[str], completeness: float = 0.0, contamination: float = 100.0
) -> List[str]:
    """
    Filter genomes according to completeness and contamination criteria

    :param checkm_tables:   List of paths to the CheckM output tables
    :param completeness:    Minimum allowed completeness
    :param contamination:   Maximum allowed contamination
    :return:                The list of genomes that passed the quality-control criteria
    """

    # Define the list of genomes that passed the quality control
    genomes = list()

    # Iterate over the CheckM output tables
    for filepath in checkm_tables:
        if os.path.isfile(filepath):
            with open(filepath) as table:
                line_count = 0
                for line in table:
                    line = line.strip()
                    if line:
                        # Always skip the first header line
                        if line_count > 0:
                            line_split = line.split("\t")

                            # Check whether the current genome respect both the completeness and contamination criteria
                            if float(line_split[-3]) >= completeness and float(line_split[-2]) <= contamination:
                                genomes.append(line_split[0])

                        line_count += 1

    return genomes


def get_bf_density(filepath: str) -> float:
    """
    Retrieve the bloom filter density

    :param filepath:    Path to the bloom filter file
    :return:            Density of the bloom filter
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
    bfs: List[str], tmpdir: str, kmer_len: int, filter_size: Optional[int] = None, nproc: int = 1
) -> Tuple[int, int, int]:
    """
    Return kmers boundaries for a specific set of genomes defined as the minimum and
    maximum number of common kmers among all the genomes in the current taxonomic level

    :param genomes:     List with paths to the bloom filter representations of the genomes
    :param tmpdir:      Path to the temporary folder
    :param kmer_len:    Kmer length
    :param filter_size: Bloom filter size
    :param nproc:       Make it parallel
    :return:            Return a tuple with boundaries
                        Total number, minimum, and maximum amount of kmers in common among the input genomes
    """

    # Search for the minimum and maximum number of common kmers among all the input genomes
    kmers = 0
    minv = np.Inf
    maxv = 0

    # Compute the number of kmers in common between all pairs of genomes
    bfdistance_intersect = bfaction(
        bfs, tmpdir, kmer_len, filter_size=filter_size, nproc=nproc, action="bfdistance", mode="intersect"
    )

    # Iterate over the bloom filters
    for i in range(len(bfs)):
        for j in range(i + 1, len(bfs)):
            # Get genome file names
            _, genome1, _, _ = get_file_info(bfs[i])
            _, genome2, _, _ = get_file_info(bfs[j])

            # Result is under key "result"
            # It returns a float by default
            common = int(bfdistance_intersect[genome1][genome2])

            if common == 0:
                # This could be only due to a wrong classification and must be reported
                with open(os.path.join(tmpdir, "zero_common_kmers.tsv"), "a+") as zck:
                    zck.write("{}\t{}\n".format(genome1, genome2))

                # Pass to the next comparison
                continue

            # Update the minimum and maximum number of common kmers
            if common > maxv:
                maxv = common
            if common < minv:
                minv = common

    # Use bfoperate --or (union) to retrieve the total number of kmers
    bfoperate_or = bfaction(bfs, tmpdir, kmer_len, filter_size=filter_size, nproc=nproc, action="bfoperate", mode="or")

    # Result is under the key "result"
    kmers = bfoperate_or["result"]

    return kmers, minv, maxv


def get_file_info(filepath: str, check_supported: bool = True, check_exists: bool = True) -> Tuple[str, str, str, str]:
    """
    Get file path, name, extension, and compression

    :param filepath:    Path to the input file
    :return:            File path, name, extension, and compression
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
    """
    Retrieve boundaries for a given taxonomic label

    :param boundaries:  Boundaries table produced by the boundaries module
    :param taxonomy:    Taxonomic label
    :return:            Taxonomy-specific boundaries
    """

    minv = 0
    maxv = 0

    # Keep track of the min and max common kmers
    min_bounds = list()
    max_bounds = list()

    # Try searching for boundaries again with a redefined taxonomic label
    retry = False

    while minv == 0 and maxv == 0:
        taxonomic_boundaries = dict()

        if not retry:
            if taxonomy in boundaries:
                # Exact search of the current taxonomic label in boundaries
                taxonomic_boundaries[taxonomy] = boundaries[taxonomy]

        else:
            for tax in boundaries:
                if tax.startswith("{}|".format(taxonomy)):
                    # Expand the search to all the taxonomies with a common prefix
                    taxonomic_boundaries[tax] = boundaries[tax]

        if taxonomic_boundaries:
            # In case the current taxonomy is in the boundaries file
            for tax in taxonomic_boundaries:
                min_bounds.append(taxonomic_boundaries[tax]["min_kmers"])
                max_bounds.append(taxonomic_boundaries[tax]["max_kmers"])

            minv = int(sum(min_bounds) / len(min_bounds))
            maxv = int(sum(max_bounds) / len(max_bounds))

        else:
            # Split the taxonomic label into levels
            taxonomic_levels = taxonomy.split("|")

            if len(taxonomic_levels) == 1:
                # Get out of the loop of there are no other levels available
                break

            # Redefine the taxonomic label
            taxonomy = "|".join(taxonomic_levels[:-1])

            # Retry
            retry = True

    return minv, maxv


def howdesbt(
    level_dir: str,
    extension: str = "fna.gz",
    kmer_len: int = 21,
    min_occurrences: int = 2,
    filter_size: int = 10000,
    nproc: int = 1,
    flat_structure: bool = False,
) -> None:
    """
    Run HowDeSBT on a specific taxonomic level
    Genomes must be in the "genomes" folder under level_dir

    :param level_dir:       Path to the taxonomic level folder
    :param extension:       Input file extension
    :param kmer_len:        Length of the kmers
    :param min_occurrences: Exclude kmers with a number of occurrences less than this param
    :param filter_size:     Size of the bloom filters
    :param nproc:           Make it parallel
    :param flat_structure:  Genomes are not taxonomically organized
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


def init_logger(filepath: Optional[str] = None, toolid: Optional[str] = None, verbose: bool = True) -> Optional[Logger]:
    """
    Define a logger to print on console, on file, or both

    :param filepath:    Path to the log file
    :param verbose:     Print on screen
    :return:            Logger object or None
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


def integrity_check(filepath) -> bool:
    """
    This is for Gzipped files only

    :param filepath:    Path to the Gzipped file
    :return:            True if it passes the integrity check
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


def load_boundaries(boundaries_filepath: str) -> Dict[str, Dict[str, Union[int, float]]]:
    """
    Load the table produced by the boundaries module

    :param boundaries_filepath: Path to the boundaries table
    :return:                    Dictionary with the table content indexed by taxa
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
                        "clusters": int(line_split[1]),
                        "references": int(line_split[2]),
                        "all_kmers": int(line_split[3]),
                        "min_kmers": int(line_split[4]),
                        "max_kmers": int(line_split[5]),
                        "min_score": float(line_split[6]),
                        "max_score": float(line_split[7]),
                    }

    return boundaries


def load_input_table(filepath: str, input_extension: str = "fna.gz") -> Dict[str, str]:
    """
    Load the input table with the list of paths to the genome files and eventually their taxonomic labels

    :param filepath:            Path to the input file
    :param input_extension:     Input genome files extension
    :return:                    A list with paths to the input genomes in case of MAGs or a dictionary with 
                                genome paths and their taxonomic labels in case of reference genomes
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


def load_manifest(manifest_filepath: str) -> Dict[str, Union[str, int, float]]:
    """
    Load the manifest file

    :param manifest_filepath:   Path to the manifest file
    :return:                    Dictionary with manifest data
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
    minv: Optional[Union[int, float]] = None,
    maxv: Optional[Union[int, float]] = None,
) -> Callable:
    """
    Take full control of input numeric types by defining custom intervals
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
    genomes: List[str],
    kl: int,
    tmpdir: str,
    closely_related: bool = False,
    nproc: int = 1,
    threads: int = 1
) -> int:
    """
    Given a set of genomes, try to define the best k-mer length with kitsune

    :param genomes:         List with genome file paths (Gzip compressed or not)
    :param kl:              kitsune tests different k-mer lengths, starting from k=4 up to kl
    :param tmpdir:          Path to the temporary folder
    :param closely_related: For closesly related genomes use this flag
    :param nproc:           Max number of processes
    :param threads:         Max number of threads
    :return:                Optimal k-mer length
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

    if not genomes_paths:
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


def println(message: str, logger: Optional[Logger] = None, verbose: bool = True) -> None:
    """
    Send messages to the logger
    It will print messages on screen, send messages to the log file, or both

    :param message:     Custom message
    :param logger:      Logger object
    :param verbose:     Print messages on sceeen if True and logger is None
    """

    if logger:
        # Redirect messages to the logger
        logger.info(message)

    elif verbose:
        # In case the logger is not defined
        # Redirect messages to the screen
        print(message)


def run(
    cmdline: List[Union[str, int, float]],
    stdout: Union[int, TextIO] = sys.stdout,
    stderr: Union[int, TextIO] = sys.stderr,
    silence: bool = False,
    extended_error: bool = False,
    retries: int = 1,
) -> None:
    """
    Wrapper for the subprocess.check_call function

    :param cmdline:         Command line list
    :param stdout:          Standard output
    :param stderr:          Standard error
    :param silence:         Redirect stdout and stderr to /dev/null
    :param extended_error:  Raise errors with traceback in case of unexpected exceptions
    :param retries:         Try running the process again in case of errors
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


def validate_url(url: str) -> bool:
    """
    Validate a URL

    :param url: Input URL to be validated
    :return:    True if validated, False otherwise
    """

    regex = re.compile(
        r'^(?:http|ftp)s?://'  # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain
        r'localhost|'  # localhost
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # or ip
        r'(?::\d+)?' # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE
    )

    return re.match(regex, url) is not None
