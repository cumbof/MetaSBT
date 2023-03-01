"""
Utility functions
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Mar 1, 2023"

import argparse as ap
import errno
import logging
import math
import os
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
    filter_size: Optional[int] = None,
    nproc: int = 1,
    action: str = "bfdistance",
    mode: str = "theta",
) -> float:
    """
    bfdistance and bfoperate wrapper

    :param genomes:     List with paths to the genome or bloom filter files
    :param tmpdir:      Path to the temporary folder with bloom filters
    :param kmer_len:    Kmer length
    :param filter_size: Bloom filter size
    :param nproc:       Make it parallel
    :param action:      "bfoperate" or "bfdistance"
    :param mode:        bfoperate modes: "and", "or", "xor", "eq", and "not"
                        bfdistance modes: "hamming", "intersect", "union", and "theta"
    :return:            Dictionary with the result of bfdistance or bfoperate
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
        filter_size = estimate_bf_size(genomes, kmer_len, os.path.join(tmpdir, "genomes"), tmpdir, nproc=nproc)

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
                os.symlink(
                    genome_path,
                    os.path.join(tmpdir, os.path.basename(genome_path)),
                )

            else:
                # Uncompress the genome file
                # It can always be Gzip compressed here
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
                    "--min=2",
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
                            str(nproc),
                            "-x",
                            file_extension,
                            "--pplacer_threads",
                            str(pplacer_threads),
                            "--tab_table",
                            "-f",
                            table_path,
                            str(bins_folder),
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
    unknown_label: str = "MSBT",
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
    :param unknown_label:       Prefix label of the newly defined clusters
    :param nproc:               Make bfdistance parallel
    :return:                    Return the assignments as a dictionary <genome_path, taxonomy>
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

    # Retrieve both the kmer length, filter size, and unknown counter from the manifest file
    manifest = load_manifest(manifest_filepath)

    try:
        # --kmer-len and --filter-size must be in manifest
        kmer_len = manifest["kmer_len"]
        filter_size = manifest["filter_size"]

        # This could be the first time --unknown-counter appears in manifest
        unknown_counter_manifest = manifest["unknown_counter"] if "unknown_counter" in manifest else 0

        # Check whether the kmer length and the filter size have been successfully retrieved
        if kmer_len == 0 or filter_size == 0:
            raise Exception('Unable to retrieve "--kmer-len" and "--filter-size" in {}'.format(manifest_filepath))

    except Exception as ex:
        raise Exception("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath)).with_traceback(
            ex.__traceback__
        )

    # Retrieve the list of input genomes
    genomes = [get_file_info(genome_path)[1] for genome_path in genomes_list]

    # Extract all the levels from the taxonomic labels in boundaries
    levels_in_boundaries = set()
    for taxonomy in boundaries:
        for level in taxonomy.split("|"):
            levels_in_boundaries.add(level)

    # Initialise variable for counting the unknown clusters
    unknown_counter = unknown_counter_manifest

    # Define the list of taxonomic levels for sorting profiles
    levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    # Keep track of the already assigned genomes
    assigned_taxa: Dict[str, List[str]] = dict()
    assigned_genomes: List[str] = list()

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
                level2score = dict()
                with open(profile) as file:
                    for line in file:
                        line = line.strip()
                        if line:
                            if not line.startswith("#"):
                                line_split = line.split("\t")
                                # Key: taxonomic level
                                # Value: kmers in common with the taxonomic level
                                level2score[line_split[2]] = int(line_split[3].split("/")[0])

                # Define the assignment
                assignment: List[str] = list()

                # Keep track of the boundary for the last identified level
                # At the end of the loop below, these numbers corresponds to the boundaries of the species
                last_known_level_mink = -1

                for level in sorted(level2score.keys(), key=lambda lv: levels.index(lv[0])):
                    # Compose the whole taxonomic label up to the current level
                    taxonomy = "{}|{}".format("|".join(assignment), level) if assignment else level

                    # Retrieve the boundaries for the current taxonomy
                    if taxonomy in boundaries:
                        # Search in the boundaries dictionary
                        mink = boundaries[taxonomy]["min_kmers"]

                    else:
                        # In case the taxonomy does not appear in the boundaries dictionary
                        # Come back to the higher level and search for it
                        higher_tax = "|".join(taxonomy.split("|")[:-1])

                        while higher_tax not in levels_in_boundaries and higher_tax:
                            # Keep truncating levels if they do not appear in the boudaries dictionary
                            higher_tax = "|".join(higher_tax.split("|")[:-1])

                        if higher_tax in levels_in_boundaries:
                            # Define min boundary
                            mink = -1

                            while mink < 0 and higher_tax:
                                # Compute the average minimum common kmers among all the
                                # genomes in all the taxonomic levels that match this search criteria
                                all_mink: List[int] = list()

                                for tax in boundaries:
                                    # In case the current taxonomy in boundaries contain the specific taxonomic level
                                    tax_level = "|{}__".format(taxonomy.split("|")[-1][0])

                                    if higher_tax in tax and tax_level in tax:
                                        # Keep track of the boundaries for computing the avarage values
                                        all_mink.append(boundaries[tax]["min_kmers"])

                                if all_mink:
                                    # Compute the boundaries
                                    mink = int(math.ceil(sum(all_mink) / len(all_mink)))

                                else:
                                    # Keep truncating levels
                                    higher_tax = "|".join(higher_tax.split("|")[:-1])

                            if mink < 0:
                                # In case computing boundaries for the current taxonomy is not possible
                                # This should never happen
                                raise Exception("Unable to assign genome {}".format(genomes[i]))

                    # At this point the boudaries are defined
                    if level2score[level] >= mink:
                        # In case the score for the current level falls in the boundaries interval
                        # Set the assignment
                        assignment.append(level)

                        # Keep track of the boundary for this last identified level
                        last_known_level_mink = mink

                    else:
                        # Does not make sense to continue with the other lower taxonomic level
                        break

                # Fill the assignment with missing levels
                assigned_levels = len(assignment)

                for pos in range(assigned_levels, len(levels)):
                    # Create new clusters
                    assignment.append("{}__{}{}".format(levels[pos][0], unknown_label, unknown_counter))

                    # Increment the unknown counter
                    unknown_counter += 1

                # Compose the assigned (partial) label
                assigned_label = "|".join(assignment)

                # Assigne current genome to the defined taxonomy
                if assigned_label not in assigned_taxa:
                    assigned_taxa[assigned_label] = list()
                assigned_taxa[assigned_label].append(genomes_list[i])

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
                        assigned_taxa[assigned_label].append(genomes_list[j])

    # Update the manifest with the new unknown counter
    if unknown_counter > unknown_counter_manifest:
        if "unknown_counter" not in manifest:
            # Append the --unknown-counter info to the manifest file
            with open(manifest_filepath, "a+") as manifest_file:
                manifest_file.write("--unknown-counter {}\n".format(unknown_counter))

        else:
            # Load the manifest file
            with open(manifest_filepath) as manifest_file:
                manifest_lines = manifest_file.readlines()

            # Update the --unknown-counter info
            with open(manifest_filepath, "w+") as manifest_file:
                for line in manifest_lines:
                    line = line.strip()
                    if line:
                        line_split = line.split(" ")
                        if line_split[0] == "--unknown-counter":
                            line_split[-1] = str(unknown_counter)
                        manifest_file.write("{}\n".format(" ".join(line_split)))

    # Mapping genome - taxonomy
    assignment = dict()

    # Dumpt the new assignments to the output file
    with open(outpath, "w+") as out:
        # Add header line
        out.write("# Genome\tAssignment\n")

        for taxonomy in sorted(assigned_taxa.keys()):
            for genome_path in sorted(assigned_taxa[taxonomy]):
                # Get genome name
                _, genome, _, _ = get_file_info(genome_path)

                out.write("{}\t{}\n".format(genome, taxonomy))

                # Take track of mapping genome - taxonomy
                assignment[genome_path] = taxonomy

    return assignment


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

            if bfdistance_theta[genome1][genome2] >= similarity:
                # Take track of excluded genomes
                filtered_genomes.append(genomes[i])

                # Also take note if the excluded genomes in the filtered file
                with open(filtered_genomes_filepath, "a+") as f:
                    f.write("{}\n".format(genomes[i]))

                break

    # Redefine the list of genomes by removing the filtered ones
    genomes = list(set(genomes).difference(set(filtered_genomes)))

    return genomes


def download(url: str, folder: str) -> str:
    """
    Download a file from URL to the specified folder

    :param url:     Source file URL
    :param folder:  Target destination folder path
    :return:        Path to the downloaded file
    """

    # Check whether the destination folder path exists
    if not os.path.isdir(folder):
        os.makedirs(folder, exist_ok=True)

    # Download file from URL to the destination folder
    run(
        ["wget", "-N", url, "-P", folder],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    return os.path.join(folder, url.split(os.sep)[-1])


def estimate_bf_size(genomes: str, kmer_len: int, prefix: str, tmp_dir: str, nproc: int = 1) -> int:
    """
    Estimate the bloom filter size with ntCard

    :param genomes:     List of paths to the genome files
    :param kmer_len:    Length of the kmers
    :param prefix:      Prefix of the output histogram file
    :param tmp_dir:     Path to the temporary folder
    :param nproc:       Make it parallel
    :return:            The estimated bloom filter size
    """

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
                "--pref={}".format(prefix),
                "@{}".format(genomes_file.name),
            ],
            silence=True,
        )

    F0 = 0
    f1 = 0

    # Read the ntCard output hist file
    with open("{}_k{}.hist".format(prefix, kmer_len)) as histfile:
        F0_found = False
        f1_found = False
        for line in histfile:
            line = line.strip()
            if line:
                if line.startswith("F0"):
                    F0 = int(line.split()[-1])
                    F0_found = True
                elif line.startswith("1"):
                    f1 = int(line.split()[-1])
                    f1_found = True
                if F0_found and f1_found:
                    break

    # Estimate the bloom filter size
    return F0 - f1


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
            if taxonomy in taxonomic_boundaries:
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
    kmer_len: int = 21,
    filter_size: int = 10000,
    nproc: int = 1,
    flat_structure: bool = False,
) -> None:
    """
    Run HowDeSBT on a specific taxonomic level
    Genomes must be in the "genomes" folder under level_dir

    :param level_dir:       Path to the taxonomic level folder
    :param kmer_len:        Length of the kmers
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
            # Genomes are always Gzip compressed .fna files here
            for genome_path in Path(genomes_folder).glob("*.fna.gz"):
                # Retrieve genome file info
                _, genome_name, extension, _ = get_file_info(genome_path)

                # Define the path to the bloom filter representation of the genome
                bf_filepath = os.path.join(filters_dir, "{}.bf".format(genome_name))

                if not os.path.isfile(bf_filepath) and not os.path.isfile("{}.gz".format(bf_filepath)):
                    # Define the uncompressed genome path
                    genome_file = os.path.join(genomes_folder, "{}{}".format(genome_name, extension))

                    # Uncompress the genome file
                    with open(genome_file, "w+") as file:
                        run(["gzip", "-dc", genome_path], stdout=file, stderr=file)

                    # Build the bloom filter file from the current genome
                    run(
                        [
                            "howdesbt",
                            "makebf",
                            "--k={}".format(kmer_len),
                            "--min=2",
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

                    # Get rid of the uncompressed genome file
                    os.unlink(genome_file)

                    # Compress the bloom filter file
                    with open("{}.gz".format(bf_filepath), "w+") as file:
                        run(["gzip", "-c", bf_filepath], stdout=file, stderr=file)

                    # Increment the genomes counter
                    how_many += 1

                elif os.path.isfile("{}.gz".format(bf_filepath)):
                    # Uncompress the bloom filter file
                    with open(bf_filepath, "w+") as file:
                        run(["gzip", "-dc", "{}.gz".format(bf_filepath)], stdout=file, stderr=file)

                # Take track of the current bloom filter file path
                with open(level_list, "a+") as file:
                    file.write("{}\n".format(bf_filepath))

    else:
        # Find all the other taxonomic levels
        # Define the new list of bloom filters
        for level in os.listdir(level_dir):
            if os.path.isdir(os.path.join(level_dir, level)):
                # Defile the path to the bloom filter file
                bf_filepath = os.path.join(level_dir, level, "{}.bf".format(level))

                if os.path.isfile(bf_filepath):
                    with open(level_list, "a+") as file:
                        file.write("{}\n".format(bf_filepath))

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

    # In case of species level or flat structure
    # Remove the uncompressed version of the bloom filter files
    if os.path.basename(level_dir).startswith("s__") or flat_structure:
        bf_filepaths = [bf.strip() for bf in open(level_list).readlines() if bf.strip()]

        for bf in bf_filepaths:
            os.unlink(bf)

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
                "-kl",
                str(kl),
                "--canonical",
                "--fast",
                "-n",
                str(nproc),
                "-t",
                str(threads),
                "-o",
                outres.name,
                "--closely_related" if closely_related else ""
            ],
            silence=True
        )

        # Get kitsune output message
        out_content = open(outres.name).read().strip()

        try:
            # Try to retrieve the optimal k
            return int(outcontent.split(" ")[-1])

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
    cmdline: List[str],
    stdout: Union[int, TextIO] = sys.stdout,
    stderr: Union[int, TextIO] = sys.stderr,
    silence: bool = False,
    extended_error: bool = False,
) -> None:
    """
    Wrapper for the subprocess.check_call function

    :param cmdline: Command line list
    :param stdout:  Standard output
    :param stderr:  Standard error
    """

    # Check whether ther is something to run
    if cmdline:
        try:
            # In case of silence
            if silence:
                # Redirect the stdout and stderr to /dev/null
                stdout = subprocess.DEVNULL
                stderr = subprocess.DEVNULL

            # Run a specific command line and redirect the stdout and stderr
            # to those specified in input
            subprocess.check_call(cmdline, stdout=stdout, stderr=stderr)

        except subprocess.CalledProcessError as e:
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

    else:
        # There is nothing to run
        raise Exception("Empty command line!")
