__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jul 6, 2022"

import sys, os, io, errno, logging, math, subprocess, shutil
import numpy as np
from pathlib import Path
from logging import Logger
from logging.config import dictConfig
from typing import List, Tuple

def checkm(genomes_paths: List[str], tmp_dir: str, file_extension: str="fna.gz", nproc: int=1, pplacer_threads: int=1) -> List[str]:
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
        
        # Retrieve the input genome extension

        # Iterate over the list of paths to the genome files
        for genome_path in genomes_paths:
            # Reorganise genomes in chunks with 1000 genomes at most
            if counter % 1000 > 0:
                counter = 0
                run_id += 1
                os.makedirs(os.path.join(run_tmp_dir, "bins_{}".format(run_id)), exist_ok=True)
            
            # Symlink genome files to the bins folder of the current chunk
            os.symlink(genome_path,
                       os.path.join(run_tmp_dir, "bins_{}".format(run_id), os.path.basename(genome_path)))
        
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
                    run(["checkm", "lineage_wf", "-t", str(nproc),
                                                 "-x", file_extension,
                                                 "--pplacer_threads", str(pplacer_threads),
                                                 "--tab_table", "-f", table_path,
                                                 str(bins_folder), run_dir],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    
                    # Add the output table path to the output list
                    output_tables.append(table_path)
                except:
                    pass

    return output_tables

def cluster(kmer_matrix_filepath: str, boundaries_filepath: str, manifest_filepath: str, profiles_dir: str, outpath: str,
            unknown_label: str="MSBT") -> str:
    """
    Define new clusters with the unassigned MAGs

    :param kmer_matrix_filepath:    Path to the kmers matrix (with header line) computed with kmtricks
    :param boundaries_filepath:     Path to the file with the taxonomic boundaries defined by the boudaries module
    :param manifest_filepath:       Path to the manifest file
    :param profiles_dir:            Path to the temporary folder with the genomes profiles defined by the profile module
    :param outpath:                 Path to the output file with the new assignments
    :param unknown_label:           Prefix label of the newly defined clusters
    :return:                        Return the path to the output table with the new assignments
    """

    # Check whether the output file already exists
    if it_exists(outpath, path_type="file"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), outpath)
    # Touch the output file
    # This is required in order to avoid raising an exception while checking for the input parameters
    with open(outpath, "x") as out:
        pass

    # Check whether the input parameters exists on file system
    for param in locals().values():
        # In this case, parameter values are always file and folder paths
        if not it_exists(param, path_type="file"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), param)
    
    # Retrieve the list of input genomes
    genomes = list()
    with open(kmer_matrix_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                # Fields are separated by a space
                genomes = line.split(" ")[1:]
                break
    
    # Load boundaries
    boundaries = dict()
    with open(boundaries_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                # Skip the header lines
                if not line.startswith("#"):
                    line_split = line.split("\t")
                    if line_split[0] not in boundaries:
                        boundaries[line_split[0]] = dict()
                    # Get the minimum and maximum common kmers among all the genomes
                    # in the current taxonomic level
                    boundaries[line_split[0]]["min"] = int(line_split[1])
                    boundaries[line_split[0]]["max"] = int(line_split[2])

    # Extract all the levels from the taxonomic labels in boundaries
    levels_in_boundaries = set()
    for taxonomy in boundaries:
        for level in taxonomy.split("|"):
            levels_in_boundaries.add(level)
    
    # Retrieve the unknown counter from the manifest file
    unknown_counter_manifest = 0
    unknown_counter_found = False
    with open(manifest_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                if line.startswith("--unknown-counter"):
                    unknown_counter_manifest = int(line.split(" ")[-1])
                    unknown_counter_found = True
    # Initialise variable for counting the unknown clusters
    unknown_counter = unknown_counter_manifest

    # Load the kmers matrix
    # Skip the header line with the list of genomes
    matrix = load_matrix(kmer_matrix_filepath, skiprows=1)

    # Define the list of level IDs for sorting profiles
    levels = ["k", "p", "c", "o", "f", "g", "s"]

    # Keep track of the already assigned MAGs
    assigned_taxa = dict()
    assigned_genomes = list()

    # Iterate over genomes
    for i, row in enumerate(matrix):
        if genomes[i] not in assigned_genomes:
            # Retrieve the genome profile
            profile = os.path.join(profiles_dir, "{}__profiles.tsv".format(genomes[i]))
            # Check whether the profile exists
            if it_exists(profile, path_type="file"):
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
                assignment = list()
                # Keep track of the boundaries for the last identified level
                # At the end of the loop below, these numbers corresponds to the boundaries of the species
                last_known_level_mink = -1
                last_known_level_maxk = -1

                for level in sorted(level2score.keys(), key=lambda l: levels.index(l[0])):
                    # Compose the whole taxonomic label up to the current level
                    taxonomy = "{}|{}".format("|".join(assignment), level) if assignment else level
                    # Retrieve the boundaries for the current taxonomy
                    if taxonomy in boundaries:
                        # Search in the boundaries dictionary
                        mink = boundaries[taxonomy]["min"]
                        maxk = boundaries[taxonomy]["max"]
                    else:
                        # In case the taxonomy does not appear in the boundaries dictionary
                        # Come back to the higher level and search for it
                        higher_tax = "|".join(taxonomy.split("|")[:-1])
                        while higher_tax not in levels_in_boundaries and higher_tax:
                            # Keep truncating levels if they do not appear in the boudaries dictionary
                            higher_tax = "|".join(higher_tax.split("|")[:-1])
                        
                        if higher_tax in levels_in_boundaries:
                            # Define boundaries
                            mink = -1
                            maxk = -1
                            while mink < 0 and maxk < 0 and higher_tax:
                                # Compute the average minimum and maximum common kmers among all the
                                # genomes in all the taxonomic levels that match this search criteria
                                all_mink = list()
                                all_maxk = list()
                                for tax in boundaries:
                                    # In case the current taxonomy in boundaries contain the specific taxonomic level
                                    if higher_tax in tax and "|{}__".format(taxonomy.split("|")[-1][0]) in tax:
                                        # Keep track of the boundaries for computing the avarage values
                                        all_mink.append(boundaries[tax]["min"])
                                        all_mink.append(boundaries[tax]["max"])
                                if all_mink:
                                    # Compute the boundaries
                                    mink = int(math.ceil(sum(all_mink)/len(all_mink)))
                                    maxk = int(math.ceil(sum(all_maxk)/len(all_maxk)))
                                else:
                                    # Keep truncating levels
                                    higher_tax = "|".join(higher_tax.split("|")[:-1])
                            
                            if mink < 0 and maxk < 0:
                                # In case computing boundaries for the current taxonomy is not possible
                                # This should never happen
                                raise Exception("Unable to assign genome {}".format(genome[i]))
                    
                    # At this point the boudaries are defined
                    if level2score[level] <= maxk and level2score[level] >= mink:
                        # In case the score for the current level falls in the boundaries interval
                        # Set the assignment
                        assignment.append(level)
                        # Keep track of the boundaries for this last identified level
                        last_known_level_mink = mink
                        last_known_level_maxk = maxk
                    else:
                        # Does not make sense to continue with the other lower taxonomic level
                        break
                
                # Fill the assignment with missing levels
                assigned_levels = len(assignment)
                for i in range(assigned_levels, len(levels)):
                    # Create new clusters
                    assignment.append("{}__{}{}".format(levels[i], unknown_label, unknown_counter))
                    # Increment the unknown counter
                    unknown_counter += 1

                # Compose the assigned (partial) label
                assigned_label = "|".join(assignment)
                # Assigne current genome to the defined taxonomy
                if assigned_label not in assigned_taxa:
                    assigned_taxa[assigned_label] = list()
                assigned_taxa[assigned_label].append(genomes[i])
                # Mark current genome as assigned
                assigned_genomes.append(genomes[i])

                # Check whether other input genomes look pretty close to the current genome by computing
                # the number of kmers in common between the current genome and all the other input genomes
                for i2, row2 in enumerate(matrix):
                    if i2 > i:
                        # Count how many times a 1 appear in the same position of both the arrays
                        common = sum([1 for pos, _ in enumerate(row) if row[pos] > 0 and row2[pos] > 0])
                        # In case this number falls into the [last_known_level_mink, last_known_level_maxk] interval
                        if common <= last_known_level_maxk and common >= last_known_level_mink:
                            # Set the second genome as assigned
                            assigned_genomes.append(genomes[i2])
                            # Also assign these genomes to the same taxonomy assigned to the current genome
                            assigned_taxa[assigned_label].append(genomes[i2])

    # Update the manifest with the new unknown counter
    if unknown_counter > unknown_counter_manifest:
        if not unknown_counter_found:
            # Append the --unknown-counter info to the manifest file
            with open(manifest_filepath, "a+") as file:
                file.write("--unknown-counter {}\n".format(unknown_counter))
        else:
            # Update the --unknown-counter info
            updated_manifest_filepath = os.path.join(os.path.dirname(manifest_filepath), "manifest2.txt")
            with open(updated_manifest_filepath, "w+") as file1:
                with open(manifest_filepath) as file2:
                    for line in file2:
                        line = line.strip()
                        if line:
                            line_split = line.split(" ")
                            if line_split[0] == "--unknown-counter":
                                line_split[-1] = str(unknown_counter)
                            file1.write("{}\n".format(" ".join(line_split)))
            # Replace the old manifest file with the updated one
            os.unlink(manifest_filepath)
            os.rename(updated_manifest_filepath, manifest_filepath)
    
    # Dumpt the new assignments to the output file
    with open(outpath, "w+") as out:
        # Add header line
        out.write("# Genome\tAssignment\n")
        for taxonomy in sorted(assigned_taxa.keys()):
            for genome in sorted(assigned_taxa[taxonomy]):
                out.write("{}\t{}\n".format(genome, taxonomy))

    return outpath

def download(url: str, folder: str) -> str:
    """
    Download a file from URL to the specified folder

    :param url:     Source file URL
    :param folder:  Target destination folder path
    :return:        Path to the downloaded file
    """

    # Check whether the destination folder path exists
    if not it_exists(folder, path_type="folder"):
        os.makedirs(folder, exist_ok=True)

    # Download file from URL to the destination folder
    run(["wget", "-N", url, "-P", folder],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    return os.path.join(folder, url.split(os.sep)[-1])

def filter_checkm_tables(checkm_tables: List[str], completeness: float=0.0, contamination: float=100.0) -> List[str]:
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
        if it_exists(filepath, path_type="file"):
            with open(filepath) as table:
                line_count = 0
                for l in table:
                    l = l.strip()
                    if l:
                        # Always skip the first header line
                        if line_count > 0:
                            l_split = l.split("\t")
                            # Check whether the current genome respect both the completeness and contamination criteria
                            if float(l_split[-3]) >= completeness and float(l_split[-2]) <= contamination:
                                genomes.append(l_split[0])
                        line_count += 1
    
    return genomes

def filter_genomes(kmer_matrix_filepath: str, outpath: str, similarity: float=100.0) -> None:
    """
    Filter genomes according to their set of kmers.
    Discard a genome if there is at least one other genome with a specific percentage of kmers in common

    :param kmer_matrix_filepath:    Path to the kmers matrix file (with header line) as output of kmtricks
    :param outpath:                 Path to the output file with the list of filtered genomes
    :param similarity:              Discard a genome if it result to have at least this percentage of common kmers with another genome
    :return:                        List of excluded genomes
    """
    
    # Check whether the output file already exists
    if it_exists(outpath, path_type="file"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), outpath)

    # Retrieve the list of input genomes
    genomes = list()
    with open(kmer_matrix_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                # Fields are separated by a space
                genomes = line.split(" ")[1:]
                break

    # Load the kmers matrix
    # Skip the header line with the list of genomes
    matrix = load_matrix(kmer_matrix_filepath, skiprows=1)

    # Take track of the excluded genomes
    excluded = list()

    # Iterate over genomes
    for i1, row1 in enumerate(matrix):
        for i2, row2 in enumerate(matrix):
            if i2 > i1 and genomes[i1] not in excluded and genomes[i2] not in excluded:
                # Count the number of kmers for the first genome
                kmers = sum([1 for k in row1 if k > 0])
                # Define the filter threshold
                threshold = int(math.ceil(kmers*similarity/100.0))
                # Count how many times a 1 appear in the same position of both the arrays
                common = sum([1 for pos, _ in enumerate(row1) if row1[i1] > 0 and row2[i2] > 0])

                # Check whether these two genomes must be dereplicated
                if common >= threshold:
                    # Add one of the two genomes to the list of exluded genomes
                    excluded.append(genomes[i2])

    # Check whether all the input genomes have been excluded
    if len(genomes) == len(excluded):
        raise Exception("All the input genomes have been excluded")
    
    # Dump the list of filtered genomes to the output file
    with open(outpath, "w+") as output:
        for genome in excluded:
            output.write("{}\n".format(genome))

def get_boundaries(kmer_matrix_filepath: str) -> Tuple[int, int]:
    """
    Return kmers boundaries for current taxonomic level defined as the minimum and
    maximum number of common kmers among all the genomes in the current taxonomic level

    :param kmer_matrix_filepath:   Path to the kmers matrix file (with no header line) as output of kmtricks
    :return:                       Return a tuple with boundaries
    """

    # Load the kmers matrix
    # It does not contain any header line
    matrix = load_matrix(kmer_matrix_filepath)

    # Search for the minimum and maximum number of common kmers among all the genomes in the kmers matrix
    kmers = 0
    minv = np.Inf
    maxv = 0

    # Iterate over genomes
    for i1, row1 in enumerate(matrix):
        for i2, row2 in enumerate(matrix):
            if i2 > i1:
                # Count how many times a value >0 appear in the same position of both the arrays
                # in the count matrix produced by kmtricks
                common = sum([1 for i, _ in enumerate(row1) if row1[i] > 0 and row2[i] > 0])
                # Update the minimum and maximum number of common kmers
                if common > maxv:
                    maxv = common
                if common < minv:
                    minv = common
        # Also take track of the total number of kmers
        if kmers == 0:
            kmers = len(row1)

    return kmers, minv, maxv

def get_level_boundaries(boundaries_filepath: str, taxonomy: str) -> Tuple[int, int]:
    """
    Retrieve boundaries for a given taxonomic label

    :param boundaries_filepath:     Path to the file with boundaries produced by the boundaries module
    :param taxonomy:                Taxonomic label
    :return:                        Boundaries
    """

    minv = 0
    maxv = 0

    # Load the boundaries file
    boundaries = dict()
    with open(boundaries_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")
                    boundaries[line_split[0]] = {
                        "min": int(line_split[3]),
                        "max": int(line_split[4])
                    }

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
                min_bounds.append(taxonomic_boundaries[tax]["min"])
                max_bounds.append(taxonomic_boundaries[tax]["max"])
            
            minv = sum(min_bounds)/len(min_bounds)
            maxv = sum(max_bounds)/len(max_bounds)
        
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

def howdesbt(level_dir: str, kmer_len: int=21, filter_size: int=10000, nproc: int=1, flat_structure: bool=False) -> None:
    """
    Run HowDeSBT on a specific taxonomic level

    :param level_dir:       Path to the taxonomic level folder
    :param kmer_len:        Length of the kmers
    :param filter_size:     Size of the bloom filters
    :param nproc:           Make it parallel
    :param flat_structure:  Genomes are not taxonomically organized
    """

    # Check whether the input folder is a valid path
    if it_exists(level_dir, path_type="folder"):
        # Extract the level name from the level folder path
        level_name = os.path.basename(level_dir)

        # Define the index folder
        index_dir = os.path.join(level_dir, "index")
        if it_exists(index_dir, path_type="folder"):
            # Remove old index folder if any
            shutil.rmtree(index_dir, ignore_errors=True)
        
        # Define the path to the file with the list of genome under the current taxonomic level
        level_list = os.path.join(level_dir, "{}.txt".format(level_name))
        if it_exists(level_list, path_type="file"):
            os.unlink(level_list)
        
        # Define the path to the bloom filter representation of the current taxonomic level
        level_filter = os.path.join(level_dir, "{}.bf".format(level_name))
        if it_exists(level_filter, path_type="file"):
            os.unlink(level_filter)

        # Define the log file
        howdesbt_log_filepath = os.path.join(level_dir, "howdesbt.log")
        howdesbt_log = open(howdesbt_log_filepath, "w+")

        # Take track of how many genomes under the specific taxonomic levels
        how_many = 0

        if os.path.basename(level_dir).startswith("s__") or flat_structure:
            # Search for all the genomes under the current taxonomic level
            genomes_folder = os.path.join(level_dir, "genomes")
            if it_exists(genomes_folder, path_type="folder"):
                # Create the filters folder
                filters_dir = os.path.join(os.path.dirname(genomes_folder), "filters")
                os.makedirs(filters_dir, exist_ok=True)

                # Iterate over the genome files
                for genome_path in Path(genomes_folder).glob("*.fna.gz"):
                    # Retrieve the genome name from the file path
                    genome_name = os.path.splitext(os.path.basename(str(genome_path)))[0]
                    if str(genome_path).endswith(".gz"):
                        genome_name = os.path.splitext(genome_name)[0]
                    
                    bf_filepath = os.path.join(filters_dir, "{}.bf".format(genome_name))
                    if not it_exists(bf_filepath, path_type="file") and not it_exists("{}.gz".format(bf_filepath), path_type="file"):
                        # Uncompress the current genome
                        genome_file = os.path.join(genomes_folder, "{}.fna".format(genome_name))
                        with open(genome_file, "w+") as file:
                            run(["gzip", "-dc", str(genome_path)], stdout=file, stderr=file)

                        # Build the bloom filter file from the current genome
                        run(["howdesbt", "makebf", "--k={}".format(kmer_len),
                                                   "--min=2",
                                                   "--bits={}".format(filter_size),
                                                   "--hashes=1",
                                                   "--seed=0,0",
                                                   genome_file,
                                                   "--out={}".format(bf_filepath),
                                                   "--threads={}".format(nproc)],
                            stdout=howdesbt_log, stderr=howdesbt_log)
                        
                        # Get rid of the uncompressed genome file
                        os.unlink(genome_file)

                        # Compress the bloom filter file
                        with open("{}.gz".format(bf_filepath), "w+") as file:
                            run(["gzip", "-c", bf_filepath], stdout=file, stderr=file)

                        # Increment the genomes counter
                        how_many += 1
                    
                    elif it_exists("{}.gz".format(bf_filepath), path_type="file"):
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
                    bf_filepath = os.path.join(level_dir, level, "{}.bf".format(level))
                    if it_exists(bf_filepath, path_type="file"):
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
            run(["howdesbt", "cluster", "--list={}".format(level_list),
                                        "--bits={}".format(filter_size),
                                        "--tree={}".format(os.path.join(index_dir, "union.sbt")),
                                        "--nodename={}".format(os.path.join(index_dir, "node{number}")),
                                        "--keepallnodes"],
                stdout=howdesbt_log, stderr=howdesbt_log)

        else:
            # With only one bloom filter it does not make sense to cluster genomes
            bf_filepath = [line.strip() for line in open(level_list).readlines() if line.strip()][0]

            # There is only one line which refers to the only bloom filter file
            shutil.copy(bf_filepath, os.path.join(level_dir, "{}.bf".format(level_name)))

            # Manually define the union.sbt file with the single node
            with open(os.path.join(index_dir, "union.sbt"), "w+") as union:
                union.write("{}\n".format(bf_filepath))
            
        # Build all the bloom filter files
        run(["howdesbt", "build", "--howde",
                                  "--tree={}".format(os.path.join(index_dir, "union.sbt")),
                                  "--outtree={}".format(os.path.join(index_dir, "index.detbrief.sbt"))],
            stdout=howdesbt_log, stderr=howdesbt_log)

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
                        file1.write("{}{}\n".format("*"*stars, node_path))
        
        # Get rid of the old tree
        os.unlink(os.path.join(index_dir, "index.detbrief.sbt"))

        # Rename the new tree
        shutil.move(os.path.join(index_dir, "index.full.detbrief.sbt"), 
                    os.path.join(index_dir, "index.detbrief.sbt"))
        
        if how_many > 1:
            # Build the bloom filter representation of the current taxonomic level
            bf_filepath = os.path.join(level_dir, "{}.bf".format(level_name))

            # Merge all the leaves together by applying the OR logic operator on the bloom filter files
            # The resulting bloom filter is the representative one, which is the same as the root node of the tree
            run(["howdesbt", "bfoperate", "--list={}".format(level_list),
                                          "--or",
                                          "--out={}".format(bf_filepath)],
                stdout=howdesbt_log, stderr=howdesbt_log)
        
        # In case of species level or flat structure
        # Remove the uncompressed version of the bloom filter files
        if os.path.basename(level_dir).startswith("s__") or flat_structure:
            bf_filepaths = [bf.strip() for bf in open(level_list).readlines() if bf.strip()]
            for bf in bf_filepaths:
                os.unlink(bf)
        
        # Close the log file handler
        howdesbt_log.close()

    else:
        raise Exception("Unable to run HowDeSBT on the following folder:\n{}".format(level_dir))

def init_logger(filepath: str=None, toolid: str=None, verbose: bool=True) -> Logger:
    """
    Define a logger to print on console, on file, or both

    :param filepath:    Path to the log file
    :param verbose:     Print on screen
    :return:            Logger object or None
    """

    # Define the logger config
    logging_config = dict(
        version = 1,
        formatters = {
            "verbose": {
                "format": "[%(toolid)s][%(levelname)s][%(asctime)s] %(message)s",
                "datefmt": "%d/%b/%Y %H:%M:%S"
            }
        },
        handlers = {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "verbose",
                "stream": sys.stdout
            },
            "file": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "INFO",
                "formatter": "verbose",
                "filename": os.devnull,
                "maxBytes": 52428800,
                "backupCount": 7
            }
        },
        loggers = {
            "console": {
                "handlers": ["console"],
                "level": logging.INFO
            },
            "file": {
                "handlers": ["file"],
                "level": logging.INFO
            },
            "full": {
                "handlers": ["console", "file"],
                "level": logging.INFO
            },
        }
    )

    # In case of log file
    if filepath:
        # Check whether its folder exists
        log_dir = os.path.dirname(filepath)
        if not it_exists(log_dir, path_type="folder"):
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

def it_exists(path: str, path_type: str="file") -> bool:
    """
    Check whether a file or folder exists on the file system

    :param path:        File or folder path
    :param path_type:   Type of the input path (file, folder)
    :return:            Always return True if the input path exists
                        Otherwise, raise an exception
    """

    if isinstance(path, str):
        if path_type.strip().lower() == "file":
            if os.path.isfile(path):
                return True
        
        elif path_type.strip().lower() == "folder":
            if os.path.isdir(path):
                return True
    
    return False

def kmtricks_matrix(genomes_fof: str, run_dir: str, kmer_len: int, filter_size: int, nproc: int, output_table: str) -> None:
    """
    Run kmtricks for building the kmers matrix

    :param genomes_fof:     Path to the fof file with the list of genomes
    :param run_dir:         Path to the working directory
    :param kmer_len:        Length of the kmers
    :param filter_size:     Size of the bloom filters
    :param nproc:           Make it parallel
    :param output_table:    Path to the output kmer matrix file
    """

    # Check whether the run folder exists
    if not it_exists(run_dir, path_type="folder"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), run_dir)

    # Initialise the kmtricks log
    # Both stdout and stderr will be redirected here
    kmtricks_log_filepath = os.path.join(run_dir, "kmtricks.log")
    kmtricks_log = open(kmtricks_log_filepath, "w+")

    # Run kmtricks for building the kmers matrix
    cmd = ["kmtricks", "pipeline", "--file", genomes_fof,
                                   "--run-dir", os.path.join(run_dir, "matrix"),
                                   "--kmer-size", str(kmer_len),
                                   "--mode", "kmer:count:bin",
                                   "--hard-min", "1",
                                   "--cpr",
                                   "--threads", str(nproc)]
    if filter_size:
        cmd.extend(["--bloom-size", str(filter_size)])
    run(cmd, stdout=kmtricks_log, stderr=kmtricks_log)

    # Aggregate partitions into a single kmer matrix
    run(["kmtricks", "aggregate", "--run-dir", os.path.join(run_dir, "matrix"),
                                  "--matrix", "kmer",
                                  "--format", "text",
                                  "--cpr-in",
                                  "--sorted",
                                  "--threads", str(nproc),
                                  "--output", output_table],
        stdout=kmtricks_log, stderr=kmtricks_log)
    
    # Close the log file handler
    kmtricks_log.close()

def load_manifest(manifest_filepath: str) -> dict:
    """
    Load the manifest file

    :param manifest_filepath:   Path to the manifest file
    :return:                    Dictionary with manifest data
    """

    if not it_exists(manifest_filepath, path_type="file"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    manifest = dict()

    with open(manifest_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                line_split = line.split(" ")
                if line_split[0] == "--kmer-len":
                    manifest["kmer_len"] = int(line_split[1])
                elif line_split[0] == "--filter-size":
                    manifest["filter_size"] = int(line_split[1])
    
    return manifest

def load_matrix(kmer_matrix_filepath: str, skiprows: int=0) -> np.ndarray:
    """
    Load a kmtricks kmers matrix into a numpy ndarray with a row for each genome
    and a column for each kmer

    :param kmer_matrix_filepath:  Path to the kmers matrix file as output of kmtricks
    :param skiprows:             Define how many lines must be skipped before loading the matrix
    :return:                     Return the kmers matrix
    """

    # Check whether the input file exists
    if not it_exists(kmer_matrix_filepath, path_type="file"):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), kmer_matrix_filepath)

    # Retrieve the number of genomes as the number of columns in the matrix
    # excluding the first one with the list of kmers
    columns = 0
    with open(kmer_matrix_filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                # Fields are separated by a space
                columns = len(line.split(" "))
                break

    # Check whether the number of genomes is greater than one
    if columns <= 1:
        raise Exception("Not enough genomes")

    # Load the whole matrix with numpy
    # Do not consider the first column
    matrix = np.loadtxt(kmer_matrix_filepath, delimiter=" ", usecols=np.arange(1, columns), skiprows=skiprows)
    # Transpose the kmers matrix
    # One row for each genome
    matrix = matrix.T

    return matrix

def number(typev, minv=None, maxv=None):
    """
    Take full control of input numeric types by defining custom intervals
    """

    def type_func(value):
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
        except:
            raise ap.ArgumentTypeError("Input value must be {}".format(typev))

    return type_func

def println(message: str, logger: Logger=None, verbose: bool=True) -> None:
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

def run(cmdline: List[str], stdout: io.TextIOWrapper=sys.stdout, stderr: io.TextIOWrapper=sys.stderr, silence: bool=False, extended_error: bool=False) -> None:
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
                stdout=subprocess.DEVNULL
                stderr=subprocess.DEVNULL

            # Run a specific command line and redirect the stdout and stderr
            # to those specified in input
            subprocess.check_call(cmdline, stdout=stdout, stderr=stderr)
        
        except subprocess.CalledProcessError as e:
            # Define the error message
            error_message = "\nAn error has occurred while running the following command:\n{}\n\n".format(" ".join(cmdline))
            
            if extended_error:
                # Extend the error message
                error_message += ("If you think this is a bug and need support, please open an Issue or a new Discussion on the official GitHub repository.\n"
                                  "We would be happy to answer your questions and help you troubleshoot any kind of issue with our framework.\n")

            raise Exception(error_message)

    else:
        # There is nothing to run
        raise Exception("Empty command line!")
