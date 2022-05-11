#!/bin/bash
#title          :update
#description    :Update the index with new genomes
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="May 5, 2022"
VERSION="0.1.0"

# Define script directory
UPDATE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${UPDATE_DIR}/utils.sh

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        -h|--help)
            # Print extended help
            QUERY_HELP=true
            source ${UPDATE_DIR}/../HELP
            exit 0
            ;;
        --license)
            # Print license
            printf "%s\n" "$(cat ${UPDATE_DIR}/../LICENSE)"
            exit 0
            ;;
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies
            exit $?
            ;;
        -v|--version)
            # Print pipeline version
            printf "update version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            printf "update: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

# TODO
