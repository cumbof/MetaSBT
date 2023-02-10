#!/usr/bin/env python3
"""
Describe your test here and locate this file under the tests folder
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Feb 8, 2023"

import unittest
from typing import List

# Define the test name
TOOL_ID = "test_template"

# Define the list of dependencies
# Add the name of an external software dependency in case it is required in your test
DEPENDENCIES: List[str] = list()


class MyUnitTest(unittest.TestCase):
    """
    Define a function for each unit test
    """

    def test_equals(self):
        """
        Just a dumb test
        """

        self.assertEqual(2 * 2, 4)


if __name__ == "__main__":
    unittest.main()
