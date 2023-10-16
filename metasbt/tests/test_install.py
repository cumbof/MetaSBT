#!/usr/bin/env python3
"""Unit tests for modules/install.py
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Mar 5, 2023"

import errno
import os
import random
import sys
import tempfile
import unittest

# Define the module name
TOOL_ID = "test_install"

# Define the list of dependencies
DEPENDENCIES = list()

# Define the test root directory
TESTS_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(TESTS_DIR), "modules")

# Define the path to the install.py
INSTALL_FILEPATH = os.path.join(MODULES_DIR, "install.py")

# Check whether install.py exists
if not os.path.isfile(INSTALL_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), INSTALL_FILEPATH)

# This is required to import the functions we need to test
sys.path.append(MODULES_DIR)

# Finally import the module
import install


class TestUtils(unittest.TestCase):
    """Unit tests for install functions.
    """

    def test_convert_size(self):
        """Test the install.convert_size() function.
        """

        with self.subTest():
            self.assertEqual(install.convert_size(1024.0, "B", "KB"), 1.0)

        with self.subTest():
            self.assertEqual(install.convert_size(1.0, "GB", "MB"), 1024.0)

    def test_fix_paths(self):
        """Test the install.fix_paths() function.
        """

        # This is the folder path that must be replaced
        original_basepath = "/original/path"

        # Replace every occurrence of the original basepath with this path
        replace_with = "/replace/with"

        with tempfile.NamedTemporaryFile() as sbt:
            with open(sbt.name, "wt") as sbt_file:
                # Define a fake sbt file with 10 random nodes
                for node_count in range(10):
                    # The number of stars before the node refers to the 
                    # position of the node into the tree
                    stars = random.randrange(0, 5)

                    nodepath = "{}{}".format(
                        "*"*stars,
                        os.path.join(original_basepath, "node{}.bf".format(node_count))
                    )

                    sbt_file.write("{}\n".format(nodepath))

            # Fix filepaths in the fake sbt file by replacing
            # every occurrence of "original_basepath" with "replace_with"
            install.fix_paths(sbt.name, "path", replace_with)

            with open(sbt.name) as sbt_file:
                for line in sbt_file:
                    line = line.strip()
                    if line:
                        # Count the number of stars before the file path
                        stars = line.count("*")

                        # Trim the stars out of the line
                        filepath = line[stars:]

                        with self.subTest():
                            self.assertEqual(
                                filepath, 
                                os.path.join(replace_with, "path", os.path.basename(filepath))
                            )


if __name__ == "__main__":
    unittest.main()
