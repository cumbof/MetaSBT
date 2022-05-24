#!/bin/bash
#title          :boundaries
#description    :Define taxonomy-specific boundaries based on kmers for the definition of new clusters
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#=====================================================================================================

DATE="May 18, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        -h|--help)
            # Print extended help
            UPDATE_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --license)
            # Print license
            printf "%s\n" "$(cat ${SCRIPT_DIR}/../LICENSE)"
            exit 0
            ;;
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies
            exit $?
            ;;
        -v|--version)
            # Print pipeline version
            printf "boundaries version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            printf "boundaries: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

printf "boundaries version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies
if [[ "$?" -gt "0" ]]; then
    exit 1
fi



PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0