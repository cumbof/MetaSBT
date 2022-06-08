#!/usr/bin/env python

__author__ = ('Fabio Cumbo (fabio.cumbo@gmail.com)')
__version__ = '0.1.0'
__date__ = 'Jun 7, 2022'

import sys, os, errno, math
import numpy as np
from inspect import getmembers, isfunction, signature
from typing import Tuple

def cluster(filepath: str, boundaries_filepath: str, manifest_filepath: str, profiles_dir: str, outpath: str) -> str:
    """
    Define new clusters with the unassigned MAGs

    :param filepath:                Path to the kmers matrix (with header line) computed with kmtricks
    :param boundaries_filepath:     Path to the file with the taxonomic boundaries defined by the boudaries module
    :param manifest_filepath:       Path to the manifest file
    :param profiles_dir:            Path to the temporary folder with the genomes profiles defined by the profile module
    :param outpath:                 Path to the output file with the new assignments
    :return:                        Return the path to the output table with the new assignments
    """

    # Check whether the output file already exists
    if os.path.exists(outpath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), outpath)
    # Touch the output file
    # This is required in order to avoid raising an exception while checking for the input parameters
    with open(outpath, "x") as out:
        pass

    # Check whether the input parameters exists on file system
    for param in locals().values():
        # In this case, parameter values are always file and folder paths
        if not os.path.exists(param):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), param)
    
    # Retrieve the list of input genomes
    genomes = list()
    with open(filepath) as file:
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
                    unknown_counter_manifest = int(line.split("=")[-1])
                    unknown_counter_found = True
    # Initialise variable for counting the unknown clusters
    unknown_counter = unknown_counter_manifest

    # Load the kmers matrix
    # Skip the header line with the list of genomes
    matrix = load_matrix(filepath, skiprows=1)

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
            if os.path.exists(profile):
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
                    assignment.append("{}__cluster{}".format(levels[i], unknown_counter))
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
                        common = sum([1 for pos, _ in enumerate(row) if row[pos]==row2[pos]==1])
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

        else:
            # Update the --unknown-counter info
            updated_manifest_filepath = os.path.join(os.path.dirname(manifest_filepath), "manifest2.txt")
            with open(updated_manifest_filepath, "w+") as file1:
                with open(manifest_filepath) as file2:
                    for line in file2:
                        line = line.strip()
                        if line:
                            line_split = line.split("=")
                            if line_split[0] == "--unknown-counter":
                                line_split[-1] = str(unknown_counter)
                            file1.write("{}\n".format("=".join(line_split)))
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

def filter_genomes(filepath: str, outpath: str) -> str:
    """
    Filter genomes according to their set of kmers.
    Discard a genome if there is at least one other genome with the same set of kmers

    :param filepath:    Path to the kmers matrix file (with header line) as output of kmtricks
    :param outpath:     Path to the output file with the list of filtered genomes
    :return:            Return the path to the output file
    """

    # Check whether the input file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)
    
    # Check whether the output file already exists
    if os.path.exists(outpath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), outpath)

    # Retrieve the list of input genomes
    genomes = list()
    with open(filepath) as file:
        for line in file:
            line = line.strip()
            if line:
                # Fields are separated by a space
                genomes = line.split(" ")[1:]
                break

    # Load the kmers matrix
    # Skip the header line with the list of genomes
    matrix = load_matrix(filepath, skiprows=1)

    # Take track of the excluded genomes
    excluded = list()

    # Iterate over genomes
    for i1, row1 in enumerate(matrix):
        for i2, row2 in enumerate(matrix):
            if i2 > i1 and genomes[i1] not in excluded and genomes[i2] not in excluded:
                # Check whether these two genomes have the same set of kmers
                if np.array_equal(row1, row2):
                    # Add one of the two genomes to the list of exluded genomes
                    excluded.append(genomes[i2])

    # Check whether all the input genomes have been excluded
    if len(genomes) == len(excluded):
        raise Exception("All the input genomes have been excluded")
    
    # Dump the list of filtered genomes to the output file
    with open(outpath, "w+") as output:
        for genome in excluded:
            output.write("{}\n".format(genome))
    
    return outpath

def get_boundaries(filepath: str) -> Tuple[int, int]:
    """
    Return kmers boundaries for current taxonomic level defined as the minimum and
    maximum number of common kmers among all the genomes in the current taxonomic level

    :param filepath:    Path to the kmers matrix file (with no header line) as output of kmtricks
    :return:            Return a tuple with boundaries
    """
    
    # Check whether the input file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

    # Load the kmers matrix
    # It does not contain any header line
    matrix = load_matrix(filepath)

    # Search for the minimum and maximum number of common kmers among all the genomes
    # in the kmers matrix
    minv = np.Inf
    maxv = 0

    # Iterate over genomes
    for i1, row1 in enumerate(matrix):
        for i2, row2 in enumerate(matrix):
            if i2 > i1:
                # Count how many times a 1 appear in the same position of both the arrays
                common = sum([1 for i, _ in enumerate(row1) if row1[i]==row2[i]==1])
                # Update the minimum and maximum common matrix
                if common > maxv:
                    maxv = common
                if common < minv:
                    minv = common

    return minv, maxv

def load_matrix(filepath: str, skiprows: int=0) -> np.ndarray:
    """
    Load a kmtricks kmers matrix into a numpy ndarray with a row for each genome
    and a column for each kmer

    :param filepath:    Path to the kmers matrix file as output of kmtricks
    :param skiprows:    Define how many lines must be skipped before loading the matrix
    :return:            Return the kmers matrix
    """

    # Retrieve the number of genomes as the number of columns in the matrix
    # excluding the first one with the list of kmers
    columns = 0
    with open(filepath) as file:
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
    matrix = np.loadtxt(filepath, delimiter=" ", usecols=np.arange(1, columns), skiprows=skiprows)
    # Transpose the kmers matrix
    # One row for each genome
    matrix = matrix.T

    return matrix

if __name__ == '__main__':
    # Define a reference to the current module
    current_module = sys.modules[__name__]
    # Define a list with all the available function in this module
    functions = getmembers(current_module, isfunction)
    # Get functions names
    functions_names = [function[0] for function in functions]

    # Check whether there are enough arguments in input
    if len(sys.argv) <= 1:
        raise Exception("Not enough arguments")
    
    # Get the target function name
    function = sys.argv[1]

    # Check whether the target function exists
    if function not in functions_names:
        raise Exception("Target function does not exist")

    # Get the target function
    target_function = functions[functions_names.index(function)][1]
    # Get the list of arguments in the signature of the target function
    target_function_parameters = signature(target_function).parameters

    # Check whether the number of input arguments matches the number of parameters
    # which defines the signature of the target function
    if len(sys.argv)-2 != len(target_function_parameters):
        raise Exception("Input arguments do not match the target function parameters")
    
    # Run the target function
    result = target_function(*sys.argv[2:])
    # Print the result to the standard output
    print(result)
