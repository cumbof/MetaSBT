#!/usr/bin/env python

__author__ = ('Fabio Cumbo (fabio.cumbo@gmail.com)')
__version__ = '0.1.0'
__date__ = 'Jun 6, 2022'

import sys, os, errno
import numpy as np
from inspect import getmembers, isfunction, signature
from typing import Tuple

def filter_genomes(filepath: str, outpath: str) -> str:
    """
    Filter genomes according to their set of kmers.
    Discard a genome if there is at least one other genome with the same set of kmers

    :param filepath:    Path to the kmers matrix file as output of kmtricks
    :param outpath:     Path to the output file with the list of filtered genomes
    :return:            Return the path to the output file
    """

    # Check whether the input file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

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

    :param filepath:    Path to the kmers matrix file as output of kmtricks
    :return:            Return a tuple with boundaries
    """
    
    # Check whether the input file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

    # Load the kmers matrix
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
