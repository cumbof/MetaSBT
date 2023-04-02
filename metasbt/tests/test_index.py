#!/usr/bin/env python3
"""
Unit tests for modules/index.py
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 1, 2023"

import errno
import hashlib
import os
import sys
import tempfile
import unittest

# Define the module name
TOOL_ID = "index"

# Define the list of dependencies
DEPENDENCIES = [
    "ncbitax2lin",
    "wget"
]

# Define the test root directory
TESTS_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(TESTS_DIR), "modules")

# Define the path to the index.py
INDEX_FILEPATH = os.path.join(MODULES_DIR, "index.py")

# Check whether index.py exists
if not os.path.isfile(INDEX_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), INDEX_FILEPATH)

# This is required to import the functions we need to test
sys.path.append(MODULES_DIR)

# Finally import the module
import index


class TestUtils(unittest.TestCase):
    """
    Unit tests for index functions
    """

    def test_download_assembly_summary(self):
        """
        Test the index.download_assembly_summary() function
        Dependencies: wget
        """

        # Define the URL to the NCBI GenBank Assembly Summary report table
        assembly_summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

        with tempfile.TemporaryDirectory() as tmpdir:
            # Download and load the assembly summary table
            assembly_table = index.download_assembly_summary(assembly_summary_url, tmpdir)

            with self.subTest():
                self.assertTrue(
                    os.path.isfile(
                        os.path.join(tmpdir, "{}.gz".format(os.path.basename(assembly_summary_url)))
                    )
                )

        with self.subTest():
            self.assertTrue(len(assembly_table) > 0)
        
    def test_download_taxdump(self):
        """
        Test the index.download_taxdump() function
        Dependencies: wget
        """

        # Define the URL to the NCBI taxdump
        taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

        with tempfile.TemporaryDirectory() as tmpdir:
            # Download nodes.dmp and names.dmp
            nodes_dmp, names_dmp = index.download_taxdump(taxdump_url, tmpdir)

            with self.subTest():
                self.assertTrue(os.path.isfile(nodes_dmp))

            with self.subTest():
                self.assertTrue(os.path.isfile(names_dmp))

    def test_load_taxa(self):
        """
        Test the index.load_taxa() function
        Dependencies: ncbitax2lin, wget
        """

        # Define the URL to the NCBI taxdump
        taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

        with tempfile.TemporaryDirectory() as tmpdir:
            # Download nodes.dmp and names.dmp
            nodes_dmp, names_dmp = index.download_taxdump(taxdump_url, tmpdir)

            with self.subTest():
                self.assertTrue(os.path.isfile(nodes_dmp))

            with self.subTest():
                self.assertTrue(os.path.isfile(names_dmp))

            # Build the ncbi_lineages table with ncbitax2lin
            ncbitax2lin_table = index.ncbitax2lin(nodes_dmp, names_dmp, tmpdir)

            with self.subTest():
                self.assertTrue(os.path.isfile(ncbitax2lin_table))

            tax_ids, tax_labels = index.load_taxa(
                ncbitax2lin_table,
                "Eukaryota",
                "Fungi",
                os.path.join(tmpdir, "taxamap.tsv")
            )

            with self.subTest():
                self.assertTrue(len(tax_ids) > 0)

            with self.subTest():
                self.assertTrue(len(tax_labels) > 0)

            with self.subTest():
                self.assertEqual(len(tax_labels), len(tax_labels))

            with self.subTest():
                self.assertTrue(os.path.isfile(os.path.join(tmpdir, "taxamap.tsv")))

    def test_get_sublist(self):
        """
        Test the index.get_sublist() function
        """

        # Define a list of genomes
        genomes = [
            "GCA_000005845.2",
            "GCA_000008865.2",
            "GCA_003018035.1",
            "GCA_003018455.1",
            "GCA_003697165.2",
            "GCA_024300685.1",
            "GCA_027925745.1",
            "GCA_027925765.1",
            "GCA_027925805.1",
            "GCA_027925825.1"
        ]

        # Define a subset of genomes by limiting items to a specific number
        genomes_subset_number = index.get_sublist(
            genomes,
            limit_number=4,
            flat_structure=True
        )

        # Define a subset of genomes by limiting items to a specific percentage
        genomes_subset_percentage = index.get_sublist(
            genomes,
            limit_percentage=40.0,
            flat_structure=True
        )

        with self.subTest():
            self.assertTrue(len(genomes_subset_number) > 0)

        with self.subTest():
            self.assertTrue(len(genomes_subset_percentage) > 0)

        with self.subTest():
            self.assertEqual(len(genomes_subset_number), len(genomes_subset_percentage))


if __name__ == "__main__":
    unittest.main()
