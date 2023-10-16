#!/usr/bin/env python3
"""Unit tests for modules/boundaries.py
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 1, 2023"

import errno
import os
import sys
import unittest

# Define the module name
TOOL_ID = "test_boundaries"

# Define the list of dependencies
DEPENDENCIES = list()

# Define the test root directory
TESTS_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(TESTS_DIR), "modules")

# Define the path to the boundaries.py
BOUNDARIES_FILEPATH = os.path.join(MODULES_DIR, "boundaries.py")

# Check whether boundaries.py exists
if not os.path.isfile(BOUNDARIES_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), BOUNDARIES_FILEPATH)

# This is required to import the functions we need to test
sys.path.append(MODULES_DIR)

# Finally import the module
import boundaries


class TestUtils(unittest.TestCase):
    """Unit tests for boundaries functions.
    """

    def test_get_lineage_from_path(self):
        """
        Test the boundaries.get_lineage_from_path() function
        """

        # Define a fake path to a species folder
        species_folder_path = \
            "/path/to/MetaSBT-DB/k__Bacteria/p__Bacillota/c__Bacilli/o__Bacillales/f__Bacillaceae/g__Bacillus/s__Bacillus_albus"

        self.assertEqual(
            "k__Bacteria|p__Bacillota|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_albus",
            boundaries.get_lineage_from_path(species_folder_path)
        )


if __name__ == "__main__":
    unittest.main()
