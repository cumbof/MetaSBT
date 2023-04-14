#!/usr/bin/env python3
"""
Unit tests for modules/utils.py
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Mar 5, 2023"

import errno
import hashlib
import os
import sys
import tempfile
import unittest

# Define the module name
TOOL_ID = "test_utils"

# Define the list of dependencies
DEPENDENCIES = [
    "ntcard",
    "wget"
]

# Define the test root directory
TESTS_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(TESTS_DIR), "modules")

# Define the path to the utils.py
UTILS_FILEPATH = os.path.join(MODULES_DIR, "utils.py")

# Check whether utils.py exists
if not os.path.isfile(UTILS_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), UTILS_FILEPATH)

# This is required to import the functions we need to test
sys.path.append(MODULES_DIR)

# Finally import the module
import utils


class TestUtils(unittest.TestCase):
    """
    Unit tests for utils functions
    """

    def test_download(self):
        """
        Test the utils.download() function
        Dependencies: wget
        """

        # Link to the E. coli K-12 reference genome from NCBI GenBank
        genome_url = "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"

        # SHA-256 computed on GCA_000005845.2_ASM584v2_genomic.fna.gz
        sha256_hash_precomp = "a8bf0111c936605eb3b8d73d5a5c1dbb64ef39c87553608a6752719995e9d9ba"

        # Retrieve the file name after download
        genome_filename = None

        # Take track of the SHA-256 hash computed on the downloaded file
        sha256_hash = None

        with tempfile.TemporaryDirectory() as tmpdir:
            # Download genome into the temporary folder
            genome_filepath = utils.download(url=genome_url, folder=tmpdir)

            # Check whether the genome file has been downloaded
            if not os.path.isfile(genome_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), genome_filepath)

            # Get genome file name
            _, genome_filename, _, _ = utils.get_file_info(genome_filepath, check_supported=True, check_exists=True)

            # Compute the SHA-256 hash
            with open(genome_filepath, "rb") as genome_file:
                sha256_hash = hashlib.sha256(genome_file.read()).hexdigest()

        with self.subTest():
            self.assertEqual(genome_filename, "GCA_000005845.2_ASM584v2_genomic")

        with self.subTest():
            self.assertEqual(sha256_hash, sha256_hash_precomp)

    def test_estimate_bf_size(self):
        """
        Test the utils.estimate_bf_size() function
        """

        # Link to the E. coli K-12 reference genome from NCBI GenBank
        genome_url = "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"

        # Bloom filter size
        filter_size = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            # Download genome into the temporary folder
            genome_filepath = utils.download(url=genome_url, folder=tmpdir)

            # Run ntCard for estimating the number of unique kmers in the input genome
            filter_size = utils.estimate_bf_size(
                [genome_filepath],
                kmer_len=21,
                prefix="genome",
                tmp_dir=tmpdir,
                nproc=1
            )

        self.assertTrue(filter_size > 0)

    def test_get_file_info(self):
        """
        Test the utils.get_file_info() function
        """

        # Test get_file_info() on this python file
        filepath = os.path.realpath(__file__)

        # It decomposes the input filepath by reporting
        # the absolute path to the file folder, file name,
        # file extension, and file compression
        _, filename, extension, compression = utils.get_file_info(filepath, check_supported=False, check_exists=False)

        with self.subTest():
            self.assertEqual(filename, "test_utils")

        with self.subTest():
            self.assertEqual(extension, ".py")

        with self.subTest():
            self.assertIs(compression, None)

    def test_validate_url(self):
        """
        Test the utils.validate_url() function
        """

        # URL to the NCBI GenBank Assembly Summary table
        assembly_summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

        # Define a malformed URL
        malformed_url = "http://test"

        with self.subTest():
            self.assertTrue(utils.validate_url(assembly_summary_url))

        with self.subTest():
            self.assertFalse(utils.validate_url(malformed_url))


if __name__ == "__main__":
    unittest.main()
