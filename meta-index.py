#!/usr/bin/env python3

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 16, 2022"

import sys

# Define the tool name
TOOL_ID = "meta-index"

# Control current Python version
# It requires Python 3 or higher
if sys.version_info[0] < 3:
    raise Exception("{} requires Python 3, your current Python version is {}.{}.{}"
                    .format(TOOL_ID, sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import os, time, errno, subprocess, requests, importlib
import argparse as ap
from pathlib import Path
from shutil import which
from modules.utils import println, run

# Define the software root directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Define the license file path
LICENSE = os.path.join(SCRIPT_DIR, "LICENSE")

# Define the modules folder
MODULES_DIR = os.path.join(SCRIPT_DIR, "modules")

# Define the paths to the file with Python requirements
REQUIREMENTS = os.path.join(SCRIPT_DIR, "requirements.txt")
# Also define the list of external software dependencies
DEPENDENCIES = [
    "checkm",                # https://github.com/Ecogenomics/CheckM
    "howdesbt",              # https://github.com/medvedevgroup/HowDeSBT
    "kmtricks",              # https://github.com/tlemane/kmtricks/
    "ncbi-genome-download",  # https://github.com/kblin/ncbi-genome-download/
    "ncbitax2lin",           # https://github.com/zyxue/ncbitax2lin
    "ntcard",                # https://github.com/bcgsc/ntCard
    "python3",               # https://www.python.org
    "wget",                  # https://www.gnu.org/software/wget/
]

# Define the software repository URLs
REPOSITORY_URL = "https://github.com/BlankenbergLab/{}".format(TOOL_ID)
RELEASES_API_URL = "https://api.github.com/repos/BlankenbergLab/{}/releases/latest".format(TOOL_ID)

def read_params():
    p = ap.ArgumentParser(prog=TOOL_ID,
                          description=("A pipeline for automatically indexing microbial genomes and accurately "
                                       "characterizing metagenome-assembled genomes with sequence bloom trees"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--check-updates",
                    action = "store_true",
                    default = False,
                    dest = "check_updates",
                    help = "Check for software updates" )
    p.add_argument( "--citations",
                    action = "store_true",
                    default = False,
                    help = ("Print software citations. "
                            "You are kindly asked to cite this software in your manuscript in case it results useful for your analyses") )
    p.add_argument( "--license",
                    action = "store_true",
                    default = False,
                    help = "Print the software license" )
    p.add_argument( "--modules",
                    action = "store_true",
                    default = False,
                    help = "List all the available modules" )
    p.add_argument( "--resolve-dependencies",
                    action = "store_true",
                    default = False,
                    dest = "resolve_dependencies",
                    help = "Check whether all the external software dependencies are available on your system" )
    p.add_argument( "-v",
                    "--version",
                    action = "version",
                    version = "{} version {} ({})".format(TOOL_ID, __version__, __date__),
                    help = "Print the current {} version and exit".format(TOOL_ID) )
    return p.parse_known_args()

def check_for_software_updates():
    """
    Check for software updates
    """

    response = requests.get(RELEASES_API_URL)
    if response.status_code == 200:
        # Load the response in dictionary
        data = response.json()
        if "tag_name" in data:
            if "v{}".format(__version__) != data["tag_name"]:
                println("A new software update is available!\n{}\n".format(REPOSITORY_URL))

def get_modules(dirpath):
    """
    Return the list of modules under the specified directory
    """

    if not os.path.exists(dirpath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), dirpath)

    modules_list = list()

    # Search for python modules
    gen = Path(dirpath).glob("*.py")
    for module_path in gen:
        module_id = os.path.splitext(os.path.basename(str(module_path)))[0]
        # Manually exclude the utilities
        if module_id != "utils":
            # Take track of the available modules
            modules_list.append(module_id)
    
    return modules_list

def print_citations():
    """
    Print citations and exit
    """
    
    println("If you are using this software for your research, please credit us in your manuscript by citing:\n")
    println("TBA\n")

def print_license():
    """
    Print the software license and exit
    """

    if not os.path.exists(LICENSE):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), LICENSE)

    with open(LICENSE, "r") as file:
        println(file.read())

def print_modules():
    """
    List all the available modules and exit
    """

    modules_list = get_modules(MODULES_DIR)

    if not modules_list:
        raise Exception("No modules available!")
    
    println("List of available modules:")
    for module_id in sorted(modules_list):
        println("\t{}".format(module_id))

def resolve_dependencies(stop_unavailable=False, verbose=True):
    """
    Check whether all the external software dependencies and Python requirements are available
    """
    
    println("Checking for software dependencies", verbose=verbose)
    # Iterate over the list of external software dependencies
    howdesbt = False
    for dependency in DEPENDENCIES:
        available = "OK" if which(dependency) is not None else "--"
        println("\t[{}] {}".format(available, dependency), verbose=verbose)
        if dependency == "howdesbt":
            howdesbt = True
        if stop_unavailable and available == "--":
            raise Exception("The external software dependency \"{}\" is not available on this system.\n"
                            "Please run \"{} --resolve_dependencies\" to verify the availability of all the required dependencies".format(dependency, TOOL_ID))

    if howdesbt:
        # Check whether the advanced bfoperate command is available with the current HowDeSBT installation
        try:
            run(["howdesbt", "bfoperate", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except:
            println(("\n[WARNING] The HowDeSBT version installed on your system does not provide advanced commands required for indexing and updating a database!\n"
                     "Please, follow the instructions on the official HowDeSBT repository on GitHub to compile the software with its alternative version of the Makefile"), 
                    verbose=verbose)

    if not os.path.exists(REQUIREMENTS):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), REQUIREMENTS)

    if verbose:
        println("\nChecking for Python requirements")
        println("\tThis will automatically install missing dependencies with pip")
        # Since checking for dependencies will automatically install the missing ones
        # Ask for the authorisation to proceed
        if input("\tDo You Want To Continue? [Y/n] ") == "Y":
            try:
                # Check whether the Python requirements are satisfied with pip
                subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", REQUIREMENTS, "--ignore-installed"])
            except:
                raise Exception("An error has occurred while running pip on {}".format(REQUIREMENTS))

def main():
    # In case of --help and --version options
    # Both the arguments are shared among the main controller and modules
    if ("--help" in sys.argv or "--version" in sys.argv) and len(sys.argv) == 3:
        # Load the list of available modules
        modules_list = get_modules(MODULES_DIR)

        if sys.argv[1] in modules_list:
            run([sys.executable, os.path.join(MODULES_DIR, "{}.py".format(sys.argv[1])), sys.argv[2]],
                extended_error=True)

        else:
            raise Exception("Unrecognised module")
    
    else:
        # Load command line parameters
        args, unknown = read_params()

        # Always print the current version
        println("{} v{} ({})\n".format(TOOL_ID, __version__, __date__))

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
            # Resolve external software dependencies and Python requirements
            resolve_dependencies(stop_unavailable=False, verbose=True)
        
        else:
            # Check for software updates
            check_for_software_updates()

            # Resolve external software dependencies
            resolve_dependencies(stop_unavailable=True, verbose=False)

            # Load the list of available modules
            modules_list = get_modules(MODULES_DIR)

            module_found = False
            for unknown_arg in unknown:
                if unknown_arg in modules_list:
                    # Build the command line
                    cmd_line = [sys.executable, os.path.join(MODULES_DIR, "{}.py".format(unknown_arg))]
                    
                    # Fix paths to the input files and folders
                    module = importlib.import_module("modules.{}".format(unknown_arg))
                    for pos in range(unknown):
                        if unknown[pos] in module.FILES_AND_FOLDERS:
                            unknown[pos+1] = os.path.abspath(unknown[pos+1])

                    # Expand the command line with all the input arguments
                    cmd_line.extend(unknown)
                    cmd_line.remove(unknown_arg)

                    try:
                        # Run the specified module
                        run(cmd_line, extended_error=True)
                    except Exception as e:
                        println(str(e))
                        sys.exit(os.EX_SOFTWARE)

                    # Mark module as found and exit
                    module_found = True

                    break

            if module_found:
                # Print citations and credits
                println("Thanks for using meta-index!\n")
                print_citations()
                println("Remember to star the meta-index repository on GitHub to stay updated on its development and new features:")
                println("https://github.com/BlankenbergLab/meta-index\n")
            else:
                raise Exception("Unrecognised module")

if __name__ == "__main__":
    main()