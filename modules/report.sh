#!/bin/bash
#title          :report
#description    :Build the database report table
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="Jun 2, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --db-dir=*)
            # Database directory
            DBDIR="${ARG#*=}"
            # Define helper
            if [[ "$DBDIR" =~ "?" ]]; then
                println "report helper: --db-dir=directory\n\n"
                println "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "$DBDIR" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
            # Trim the last slash out of the path
            DBDIR="${DBDIR%/}"
            ;;
        -h|--help)
            # Print extended help
            REPORT_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --log=*)
            # Path to the log file
            LOG_FILEPATH="${ARG#*=}"
            # Define helper
            if [[ "${LOG_FILEPATH}" =~ "?" ]]; then
                println "report helper: --log=file\n\n"
                println "\tPath to the log file.\n\n"
                exit 0
            fi
            # Remove the log file if it already exists
            if [[ -f ${LOG_FILEPATH} ]]; then
                rm ${LOG_FILEPATH}
            fi
            ;;
        --output-file=*)
            # Output folder
            OUTPUTFILE="${ARG#*=}"
            # Define helper
            if [[ "$OUTPUTFILE" =~ "?" ]]; then
                println "report helper: --output-file=file\n\n"
                println "\tThis is the path to the output table.\n\n"
                exit 0
            fi
            ;;
        -v|--version)
            # Print pipeline version
            println "report version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            println "report: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

println "report version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

# Build a report table
# For each taxonomy, report the number of MAGs and reference genomes, the number of excluded MAGs and reference genomes, 
# the median completeness and contamination for both MAGs and reference genomes
# TODO

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0