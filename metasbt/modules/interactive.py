#!/usr/bin/env python3
"""
MetaSBT Interactive CLI
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Nov 29, 2022"

import argparse as ap
import errno
import os
from typing import List, Optional

# Define the framework name
FRAMEWORK_ID = "MetaSBT"

# Define the module name
TOOL_ID = "interactive"

# Define the list of dependencies
DEPENDENCIES: List[str] = list()

# Define the list of input files and folders
FILES_AND_FOLDERS = [
    "--db-dir"  # Database folder path
]

# Define the list of commands
# TODO define new commands
COMMANDS = {
    "help": "Print the list of available commands",
    "quit": "Quit the interactive interface"
}


def read_params():
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="MetaSBT Interactive CLI",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required=True,
        dest="db_dir",
        help="This is the database directory with the taxonomically organised sequence bloom trees",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version="\"{}\" version {} ({})".format(TOOL_ID, __version__, __date__),
        help="Print the current \"{}\" version and exit".format(TOOL_ID),
    )
    return p.parse_args()


def interactive(db_dir: str) -> None:
    """
    Run the command-line interface

    :param db_dir:      Path to the database root folder
    """

    # Print the welcome message
    print("\n{} ({}) version {} ({})\n\nType \"help\" for a list of available commands\n".format(
        FRAMEWORK_ID, TOOL_ID, __version__, __date__
    ))

    while True:
        try:
            # Keep asking for a command 
            command = input("{}> ".format(FRAMEWORK_ID)).strip().lower()

            if command:
                # In case the user didn't provide an empty command
                print()

                if command in COMMANDS:
                    # Check whether the input is a valid command
                    
                    if command == "help":
                        # Print the list of available commands
                        print("{:<8} {}".format("Command", "Description"))
                        print("=======  ===========")
                        for cmd in sorted(COMMANDS.keys()):
                            print("{:<8} {:<15}".format(cmd, COMMANDS[cmd]))
                        
                    elif command == "quit":
                        # Quit the interactive interface
                        break
                    
                    print()

                else:
                    print("Command not available!\nType \"help\" for a list of available commands\n")
        
        except KeyboardInterrupt:
            # Disable keyboard interruption
            # Use the "quit" command to close the interactive shell
            print()
            pass


def main() -> None:
    # Load command line parameters
    args = read_params()

    # Check whether the input database folder path exists
    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    # TODO before running the interactive shell, check whether the provided --db-dir actually contains a valid MetaSBT database

    interactive(args.db_dir)


if __name__ == "__main__":
    main()
