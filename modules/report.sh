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

# Extract relevant information from the species folder
process_species () {
    SPECIES_FOLDER=$1   # Path to the species folder
    OUTPUT_FILE=$2      # Output table path

    # Split the species folder path
    SPECIES_FOLDER_SPLIT=($(echo ${SPECIES_FOLDER} | tr "/" " "))
    # Retrieve all the taxonomic levels from the kingdom to the species
    LEVELS=()
    COLLECT=false
    for LEVEL in "${SPECIES_FOLDER_SPLIT[@]}"; do
        # Do not collect all levels in the array
        # It contains the full path structure
        if [[ "$LEVEL" = k__* ]]; then
            # Start collecting levels from the kingdom
            COLLECT=true
        fi
        # In case it found the kingdom
        # Collect the taxonomic levels
        if $COLLECT; then
            # Take track of the taxonomic levels
            LEVELS+=($LEVEL)
        fi
    done

    # In case it was able to detect all the 7 taxonomic levels
    if [[ "${#LEVELS[@]}" -eq "7" ]]; then
        # Rebuild the lineage
        LINEAGE="$(printf "|%s" "${LEVELS[@]}")"
        # Trim the first pipe out of the lineage
        LINEAGE=${LINEAGE:1}

        # Retrieve the total number of genomes assigned to the current species
        GENOMES="$(cat "${SPECIES_FOLDER}/genomes.fof" | wc -l)"
        
        # Keep track of all the completeness and contamination percentages 
        # for all the genomes in the current species
        COMPLETENESS_ARR=()
        CONTAMINATION_ARR=()
        
        # Count the number of MAGs and reference genomes
        MAGS=0
        REFERENCES=0
        # Iterate over the genome IDs
        for GENOMENAME in $(cut -f1 "${SPECIES_FOLDER}/genomes.fof"); do
            FOUND=false
            # Search for the current genome in both references.txt and mags.txt
            if grep -q "^$GENOMENAME$" "${SPECIES_FOLDER}/references.txt"; then
                # Increment the reference genomes counter
                REFERENCES=$(($REFERENCES + 1))
                FOUND=true
            elif grep -q "^$GENOMENAME$" "${SPECIES_FOLDER}/mags.txt"; then
                # Increment the MAGs counter
                MAGS=$(($MAGS + 1))
                FOUND=true
            fi
            # In case the current genome exists in references.txt or in mags.txt
            if $FOUND && [[ -f "${SPECIES_FOLDER}/checkm.tsv" ]]; then
                # Search for current genome ID into the CheckM output table
                GENOME_DATA="$(grep "$GENOMENAME"$'\t' ${SPECIES_FOLDER}/checkm.tsv)"
                if [[ ! -z "${GENOME_DATA}" ]]; then
                    # Retrieve its completeness and contamination percentage
                    COMPLETENESS_ARR+=("$(echo ${GENOME_DATA} | cut -d$'\t' -f12)")   # Get value under column 12
                    CONTAMINATION_ARR+=("$(echo ${GENOME_DATA} | cut -d$'\t' -f13)")  # Get value under column 13
                fi
            fi
        done
        
        # Compute the mean completeness and contamination for the current species
        if [[ "${#COMPLETENESS_ARR[@]}" -gt "0" ]]; then
            # Use bc command to do math with floating points numbers
            MEAN_COMPLETENESS="$(IFS="+"; bc -l <<< (${COMPLETENESS_ARR[*]})/${#COMPLETENESS_ARR[@]})"
            MEAN_CONTAMINATION="$(IFS="+"; bc -l <<< (${CONTAMINATION_ARR[*]})/${#CONTAMINATION_ARR[@]})"
        else
            # In case no completeness and contamination statistics are available 
            # for all the genomes in the current species
            MEAN_COMPLETENESS=0
            MEAN_CONTAMINATION=0
        fi

        # Write the species report to the output file
        printf "%s\t%s\t%s\t%.2f\t%.2f\n" "$LINEAGE" "$MAGS" "$REFERENCES" "${MEAN_COMPLETENESS}" "${MEAN_CONTAMINATION}"
    fi    
}

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
# For each taxonomy, report the number of MAGs and reference genomes, and the mean completeness and contamination

# Initialise the report table
printf "# Lineage\tMAGs\tReferences\tMean Completeness\tMean Contamination\n" > $OUTPUTFILE
# Search for all the species folders
find $DBDIR -type d -name "s__*" | xargs -n 1 -I {} bash -c \
    'SPECIES={}; \
     process_species "$SPECIES" '"$OUTPUTFILE"';'

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0