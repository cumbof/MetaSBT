#!/usr/bin/env python3
"""A scalable framework for automatically indexing microbial genomes and accurately
characterizing metagenome-assembled genomes with Sequence Bloom Trees
"""

import argparse as ap
import errno
import importlib
import os
import platform
import pkgutil
import subprocess
import sys
from shutil import which
from typing import List

import requests

from metasbt import __date__, __version__
from metasbt.modules.utils import println, run

# Define the tool name
TOOL_ID = "MetaSBT"

# Control platform
# It works on Darwin and Linux platforms
if platform.system() not in ["Darwin", "Linux"]:
    raise Exception("{} does not work on {} platforms".format(TOOL_ID, platform.system()))

# Control current Python version
# It requires Python 3.8 or higher
if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 8):
    raise Exception(
        "{} requires Python 3.8 or higher. Your current Python version is {}.{}.{}".format(
            TOOL_ID, sys.version_info[0], sys.version_info[1], sys.version_info[2]
        )
    )

# Define the software root directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the modules folder
MODULES_DIR = os.path.join(SCRIPT_DIR, "modules")

# Define the tests folder
TESTS_DIR = os.path.join(SCRIPT_DIR, "tests")

# Define the paths to the file with Python requirements
REQUIREMENTS = os.path.join(SCRIPT_DIR, "requirements.txt")

# Define the license URL
LICENSE = "https://raw.githubusercontent.com/cumbof/{}/main/LICENSE".format(TOOL_ID)

# Define the software repository URLs
REPOSITORY_URL = "https://github.com/cumbof/{}".format(TOOL_ID)
RELEASES_API_URL = "https://api.github.com/repos/cumbof/{}/releases/latest".format(TOOL_ID)


def read_params():
    """Read and test the input arguments.

    Returns
    -------
    argparse.ArgumentParser
        The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "A scalable framework for automatically indexing microbial genomes and accurately "
            "characterizing metagenome-assembled genomes with Sequence Bloom Trees"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--check-updates",
        action="store_true",
        default=False,
        dest="check_updates",
        help="Check for software updates",
    )
    p.add_argument(
        "--citations",
        action="store_true",
        default=False,
        help=(
            "Print software citations. "
            "You are kindly asked to cite this software in your manuscript in case it results useful for your analyses"
        ),
    )
    p.add_argument(
        "--license",
        action="store_true",
        default=False,
        help="Print the software license",
    )
    p.add_argument(
        "--modules",
        action="store_true",
        default=False,
        help="List all the available modules",
    )
    p.add_argument(
        "--resolve-dependencies",
        action="store_true",
        default=False,
        dest="resolve_dependencies",
        help="Check whether all the external software dependencies are available on your system",
    )
    p.add_argument(
        "-t",
        "--test",
        type=str,
        dest="test",
        help=(
            "Run a unit test among those reported with --tests. "
            "Type \"metasbt --test all\" for running all the available unit tests at once"
        ),
    )
    p.add_argument(
        "--tests",
        action="store_true",
        default=False,
        help="List all the available unit tests",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version="{} version {} ({})".format(TOOL_ID, __version__, __date__),
        help="Print the current {} version and exit".format(TOOL_ID),
    )
    return p.parse_known_args()


def check_for_software_updates() -> None:
    """Check for new releases through the GitHub APIs.
    Print a message on the stdout in case a new release is available.

    Raises
    ------
    requests.ConnectionError
        If it is not able to connect to the GitHub APIs.
    """

    try:
        response = requests.get(RELEASES_API_URL)

        if response.status_code == 200:
            # Load the response in dictionary
            data = response.json()

            if "tag_name" in data:
                if __version__ != data["tag_name"]:
                    println("A new release is available!")
                    println("{}\n".format(REPOSITORY_URL))
    
    except requests.ConnectionError:
        pass


def get_modules(dirpath: str, test: bool=False) -> List[str]:
    """Return the list of modules under the specified directory.

    Parameters
    ----------
    dirpath : str
        Path to the folder with modules.
    test : bool, default False
        True if the folder contains unit tests

    Raises
    ------
    FileNotFoundError
        If the input path `dirpath` does not exist on the file system.

    Returns
    -------
    list
        A list with the module names in the input folder `dirpath`.
    """

    if not os.path.isdir(dirpath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), dirpath)

    modules_list = list()

    # Search for Python modules
    packages = pkgutil.walk_packages(path=[dirpath])

    for module_info in packages:
        # Exclude utils in case test if False
        if test or (module_info.name != "utils" and not test):
            # Take track of the available modules
            modules_list.append(module_info.name)

    return modules_list


def print_citations() -> None:
    """Print citations and exit."""

    println("If you are using {} for your research, please credit us in your manuscript by citing:\n".format(TOOL_ID))
    println("TBA\n")


def print_license() -> None:
    """Print the software license and exit.

    Raises
    ------
    requests.ConnectionError
        If it is not able to retrieve the license file from GitHub.
    """

    try:
        response = requests.get(LICENSE)

        if response.status_code == 200:
            # Print the license content
            println("{}\n".format(response.text))

        else:
            println("Unable to retrieve the license from the following URL:\n{}\n\nPlease try again")

    except requests.ConnectionError:
        pass


def print_modules() -> None:
    """List all the available modules and exit.

    Raises
    ------
    Exception
        If there are no modules available.
    """

    modules_list = get_modules(MODULES_DIR)

    if not modules_list:
        raise Exception("No modules available!")

    println("List of available modules:")

    for module_id in sorted(modules_list):
        println("\t{}".format(module_id))


def print_tests() -> None:
    """List all the available unit tests and exit.

    Raises
    ------
    Exception
        If there are no unit tests available.
    """

    # Unit tests are defined into python files
    # Use the get_modules() function to retrieve the list of unit tests
    tests_list = get_modules(TESTS_DIR, test=True)

    if not tests_list:
        raise Exception("No tests available!")

    println("List of available tests:")

    for test_id in sorted(tests_list):
        println("\t{}".format(test_id))


def run_test(test: str) -> None:
    """Run unit tests

    Parameters
    ----------
    test : str
        Test name among those reported with --tests.
        Use the wildcard "all" for runnig all the available tests at once.

    Raises
    ------
    ValueError
        If the input `test` is not available as a unit test.
    Exception
        If the unit test failed to run.
    """

    # Retrieve the list of the available unit tests
    tests_list = get_modules(TESTS_DIR, test=True)

    if test.lower() != "all" and test not in tests_list:
        raise ValueError("There are no unit tests for \"{}\"".format(test))

    unit_tests = tests_list if test.lower() == "all" else [test]

    for unit_test in unit_tests:
        try:
            println("Running unit tests: {}".format(unit_test))

            # Run the specified module or unit test
            run([sys.executable, os.path.join(TESTS_DIR, "{}.py".format(unit_test))])

        except Exception as ex:
            # Print the exception and run the next unit test
            println(str(ex))


def resolve_dependencies(dependencies: List[str], stop_unavailable: bool=False, verbose: bool=True) -> None:
    """Check whether all the external software dependencies and Python requirements are available.

    Parameters
    ----------
    dependencies : list
        A list of strings with the name of the external software dependencies.
    stop_unavailable : bool, default False
        Stop the execution and raise and Exception in case on of the dependencies listed in `dependencies` is not available.
    verbose : bool = True
        Print messages on the stdout if True.

    Raises
    ------
    FileNotFoundError
        If the "requirements.txt" file with the list of Python dependencies is not available.
    Exception
        - If `stop_unavailable` is True and at least one of the dependencies listed in `dependencies` is not available;
        - If HowDeSBT has not been compiled with the alternative version of the Makefile that enables advanced commands;
        - If pip failes in installing Python dependencies listed in the "requirements.txt" file.
    """

    # Sort the list of dependencies
    dependencies = sorted(list(set(dependencies)))

    println("Checking for software dependencies", verbose=verbose)

    howdesbt = False
    # Iterate over the list of external software dependencies
    for dependency in dependencies:
        available = "OK" if which(dependency) is not None else "--"
        println("\t[{}] {}".format(available, dependency), verbose=verbose)

        if dependency == "howdesbt" and available == "OK":
            howdesbt = True

        if stop_unavailable and available == "--":
            raise Exception(
                (
                    'The external software dependency "{}" is not available on this system.\n'
                    'Please run "{} --resolve-dependencies" to verify the availability of all the required dependencies'
                ).format(dependency, TOOL_ID.lower())
            )

    println(
        (
            "\n[WARNING] Please note that the versions of verified software dependencies and "
            "those required by {} could not match. This could lead to unexpected errors"
        ).format(TOOL_ID),
        verbose=verbose
    )

    if howdesbt:
        # Check whether the advanced bfoperate command is available with the current HowDeSBT installation
        try:
            run(
                ["howdesbt", "bfoperate", "--help"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        except Exception:
            println(
                (
                    "\n[WARNING] The HowDeSBT version installed on your system does not provide "
                    "advanced commands required for indexing and updating a database!\n"
                    "Please, follow the instructions on the official HowDeSBT repository on GitHub "
                    "to compile the software with its alternative version of the Makefile"
                ),
                verbose=verbose,
            )

    if not os.path.isfile(REQUIREMENTS):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), REQUIREMENTS)

    if verbose:
        println("\nChecking for Python requirements")
        println("\tThis will automatically install missing dependencies with pip")

        # Since this will automatically install the missing dependencies
        # ask for the authorization to proceed
        if input("\tDo You Want To Continue? [Y/n] ") == "Y":
            try:
                # Check whether the Python requirements are satisfied with pip
                run(
                    [
                        sys.executable,
                        "-m",
                        "pip",
                        "install",
                        "-r",
                        REQUIREMENTS,
                        "--ignore-installed",
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            except Exception as ex:
                raise Exception("An error has occurred while running pip on {}".format(REQUIREMENTS)).with_traceback(
                    ex.__traceback__
                )


def main() -> None:
    # In case no arguments are specified
    if len(sys.argv) <= 1:
        # Run the helper
        sys.argv.append("--help")

    # Both --help and --version arguments are shared among the main controller and modules
    if ("--help" in sys.argv or "--version" in sys.argv) and len(sys.argv) == 3:
        # Load the list of available modules
        modules_list = get_modules(MODULES_DIR)

        if sys.argv[1] in modules_list:
            run(
                [
                    sys.executable,
                    os.path.join(MODULES_DIR, "{}.py".format(sys.argv[1])),
                    sys.argv[2],
                ],
                extended_error=True,
            )

        else:
            raise Exception("Unrecognised module")

    else:
        # Load command line parameters
        args, unknown = read_params()

        # Always print the current version
        println("{} v{} ({})\n".format(TOOL_ID.lower(), __version__, __date__))

        if args.check_updates:
            # Check for software updates
            check_for_software_updates()

        elif args.citations:
            # Print citations and references
            print_citations()

        elif args.license:
            # Print the software license
            print_license()

        elif args.modules:
            # Print the list of available modules
            print_modules()

        elif args.resolve_dependencies:
            # Load the list of available modules
            modules_list = get_modules(MODULES_DIR)

            # Load the list of available tests
            tests_list = get_modules(TESTS_DIR, test=True)

            # Take track of all the external software dependencies
            dependencies = list()

            # Load the list of module- and test-specific dependencies
            for level, modules in zip(["modules", "tests"], [modules_list, tests_list]):
                for module_id in modules:
                    module = importlib.import_module("{}.{}.{}".format(TOOL_ID.lower(), level, module_id))
                    dependencies.extend(module.DEPENDENCIES)

            # Remove duplicates
            dependencies = list(set(dependencies))

            # Resolve external software dependencies and Python requirements
            resolve_dependencies(dependencies, stop_unavailable=False, verbose=True)

        elif args.tests:
            # Print the list of available unit tests
            print_tests()
        
        elif args.test:
            # Run a unit test among those reported with --tests
            run_test(args.test)

        else:
            # Check for software updates
            check_for_software_updates()

            # Load the list of available modules
            modules_list = get_modules(MODULES_DIR)

            # Check whether the provided command is available among modules and unit tests
            cmd_found = False

            for unknown_arg in unknown:
                if unknown_arg in modules_list:
                    # Build the command line
                    cmd_line = [
                        sys.executable,
                        os.path.join(MODULES_DIR, "{}.py".format(unknown_arg)),
                    ]

                    # Import module or unit test
                    module = importlib.import_module("{}.{}.{}".format(TOOL_ID.lower(), "modules", unknown_arg))

                    # Resolve external software dependencies
                    resolve_dependencies(module.DEPENDENCIES, stop_unavailable=True, verbose=False)

                    # Fix paths to the input files and folders
                    for pos in range(len(unknown)):
                        if unknown[pos] in module.FILES_AND_FOLDERS:
                            # Always use absolute paths
                            unknown[pos + 1] = os.path.abspath(unknown[pos + 1])

                    # Expand the command line with all the input arguments
                    cmd_line.extend(unknown)
                    cmd_line.remove(unknown_arg)

                    if unknown_arg in ["index", "update"]:
                        # Add the main MetaSBT version
                        # Supported by the index and update modules only because of the --resume arg
                        # Display a warning in case of different MetaSBT versions
                        cmd_line.extend(["--parent-version", __version__])

                    try:
                        # Run the specified module or unit test
                        run(cmd_line, extended_error=True)

                    except Exception as ex:
                        println(str(ex))
                        sys.exit(os.EX_SOFTWARE)

                    # Mark module or test as found and exit
                    cmd_found = True

                    break

            if cmd_found:
                # Print citations and credits
                println("Thanks for using {}!\n".format(TOOL_ID))
                print_citations()
                println(
                    "Remember to star the {} repository on GitHub to stay updated on its development and new features:".format(
                        TOOL_ID
                    )
                )
                println("https://github.com/cumbof/{}\n".format(TOOL_ID))
            else:
                raise Exception("Unrecognised module \"{}\"".format(unknown_arg))


if __name__ == "__main__":
    main()
