#!/bin/bash
#title          :query
#description    :Query the sequence bloom trees at all the 7 taxonomic levels
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#============================================================================

DATE="May 4, 2022"
VERSION="0.1.0"

# Define script directory
QUERY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${QUERY_DIR}/utils.sh

# Do not expand queries by default
# If true, keep querying trees at higher taxonomic levels
EXPAND_QUERY=false

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --expand)
            # Expand the input query on all the taxonomic levels
            EXPAND_QUERY=true
            ;;
        -h|--help)
            # Print extended help
            QUERY_HELP=true
            source ${QUERY_DIR}/../HELP
            exit 0
            ;;
        --input=*)
            # Input file path with queries
            INPUT="${ARG#*=}"
            # Define helper
            if [[ "${INPUT}" =~ "?" ]]; then
                printf "query helper: --input=file\n\n"
                printf "\tThis is the input file with the set of queries.\n\n"
                exit 0
            fi
            ;;
        --license)
            # Print license
            printf "%s\n" "$(cat ${QUERY_DIR}/../LICENSE)"
            exit 0
            ;;
        --output=*)
            # Output file path with matches
            OUTPUT="${ARG#*=}"
            # Define helper
            if [[ "${OUTPUT}" =~ "?" ]]; then
                printf "query helper: --output=file\n\n"
                printf "\tThis is the output file with matches.\n\n"
                exit 0
            fi
            ;;
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies
            exit $?
            ;;
        --tree=*)
            # Tree definition file path
            TREE="${ARG#*=}"
            # Define helper
            if [[ "${TREE}" =~ "?" ]]; then
                printf "query helper: --tree=file\n\n"
                printf "\tThis is the tree definition file.\n\n"
                exit 0
            fi
            ;;
        -v|--version)
            # Print pipeline version
            printf "query version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            printf "query: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

# TODO
