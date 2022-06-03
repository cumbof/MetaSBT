#!/bin/bash
#title          :boundaries
#description    :Define taxonomy-specific boundaries based on kmers for the definition of new clusters
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#=====================================================================================================

DATE="Jun 2, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Wrapper for kmtricks
# Build a kmer matrix for establishing boundaries
# Consider clusters with a minimum number of reference genomes (100 by default)
define_boundaries () {
    LEVEL_DIR=$1    # Cluster main folder
    MIN_GENOMES=$2  # Minimum number of genomes required for running kmtricks
    TMPDIR=$3       # Temporary folder
    OUTPUT=$4       # Path to the output table with boundaries
    NPROC=$5        # Max nproc for multiprocessing
    CLEANUP=$6      # Remove temporary data at the end of the pipeline if true
    
    println "Processing %s\n" "${LEVEL_DIR}"
    # Create a temporary folder for current taxonomy
    LEVEL_NAME="$(basename ${LEVEL_DIR})"
    TMP_LEVEL_DIR=$TMPDIR/${LEVEL_NAME}
    mkdir -p ${TMP_LEVEL_DIR}
    # Search and merge reference genomes paths unders all genomes.fof files in current taxonomic level
    find ${LEVEL_DIR} -type f -iname "references.txt" -follow | xargs -n 1 -I {} bash -c \
        'REFERENCES={}; \
         SPECIES_DIR=$(dirname ${REFERENCES}); \
         while read REFNAME; do \
            REFDATA="$(grep "^$REFNAME " $(dirname $REFERENCES)/genomes.fof)"
            if [[ ! -z $REFDATA ]]; then \
                REFPATH="$(echo "$REFDATA" | rev | cut -d" " -f1 | rev)"
                echo $REFPATH >> '"${TMP_LEVEL_DIR}"'/genomes.fof;
            fi \
         done < $REFERENCES'

    if [[ -f ${TMP_LEVEL_DIR}/genomes.fof ]]; then
        # Process current level if if contains a genomes.fof file
        HOW_MANY=$(cat ${TMP_LEVEL_DIR}/genomes.fof | wc -l)
        println "\tReference genomes: %s\n" "${HOW_MANY}"
        
        # Keep processing current taxonomic level if it contains enough genomes
        if [[ "${HOW_MANY}" -ge "${MIN_GENOMES}" ]]; then
            # Run kmtricks to build the kmers matrix
            kmtricks_matrix_wrapper ${TMP_LEVEL_DIR}/genomes.fof \
                                    ${TMP_LEVEL_DIR}/matrix \
                                    $NPROC \
                                    ${TMP_LEVEL_DIR}/kmers_matrix.txt
            
            # Check whether the kmers matrix has been generated
            if [[ -f ${TMP_LEVEL_DIR}/kmers_matrix.txt ]]; then
                # Count how many columns in matrix
                HOW_MANY_COLUMNS=$(awk '{print NF; exit}' ${TMP_LEVEL_DIR}/kmers_matrix.txt)
                # Create an associative array to keep track of common kmers
                declare -A COMMON_KMERS
                # Read the kmers matrix line by line
                # Iterate over non-empty lines
                cat ${TMP_LEVEL_DIR}/kmers_matrix.txt | grep -v "^$" | while read line; do
                    KMER_LINE=($(echo $line))
                    for COLUMN1 in $(seq 2 ${HOW_MANY_COLUMNS}); do
                        for COLUMN2 in $(seq $COLUMN1 ${HOW_MANY_COLUMNS}); do
                            if [[ "$COLUMN1" -eq "$COLUMN2" ]]; then
                                # Do not compare a genome with itself
                                continue
                            fi
                            # Compare kmer absence/presence
                            # Take track of current genome couple if both the genomes have the same kmer
                            if [[ "${KMER_LINE[$COLUMN1]}" -gt "0" ]] && [[ "${KMER_LINE[$COLUMN2]}" -gt "0" ]]; then
                                # Increment counter in common kmers table
                                COMMON_KMERS[${COLUMN1}${COLUMN2}]=$((COMMON_KMERS[${COLUMN1}${COLUMN2}] + 1))
                            fi
                        done
                    done
                done

                # Check for minimum and maximum common kmers
                MIN_KMERS=-1
                MAX_KMERS=-1
                # Iterate over values in COMMON_KMERS associative array
                for count in ${COMMON_KMERS[@]}; do
                    if [[ "${MIN_KMERS}" -eq "-1" ]] && [[ "${MAX_KMERS}" -eq "-1" ]]; then
                        MIN_KMERS=$count
                        MAX_KMERS=$count
                    else
                        if [[ "$count" -lt "${MIN_KMERS}" ]]; then
                            MIN_KMERS=$count
                        elif [[ "$count" -gt "${MAX_KMERS}" ]]; then
                            MAX_KMERS=$count
                        fi
                    fi
                done

                # Report minimum and maximum common kmers to the output table
                LEVELS=($(echo ${LEVEL_DIR} | tr "\/" " "))
                LINEAGE=""
                KINGDOM_FOUND=false
                # Iterate over taxonomic levels to rebuild the full lineage
                for LEVEL in ${LEVELS[@]}; do
                    # Start keeping track of the taxonomic levels from the kingdom
                    if [[ "$LEVEL" = k__* ]]; then
                        KINGDOM_FOUND=true
                    fi
                    # Rebuild the full lineage
                    if ${KINGDOM_FOUND}; then
                        LINEAGE=$LINEAGE"|"$LEVEL
                    fi
                done
                # Remove the first pipe out of the lineage
                LINEAGE=${LINEAGE:1}
                # Dump result to the output file
                println "%s\t%s\t%s\n" "$LINEAGE" "${MIN_KMERS}" "${MAX_KMERS}" >> $OUTPUT

                # Print results
                println "\tMinimum common kmers: %s\n" "${MIN_KMERS}"
                println "\tMaximum common kmers: %s\n\n" "${MAX_KMERS}"
            else
                # kmtricks failed in building the kmers matrix
                println "\t[ERROR] An error has occurred while building the kmers matrix\n\n"
            fi

            # Cleanup temporary data
            if $CLEANUP; then
                # Removing temporary matrix folder
                rm -rf ${TMP_LEVEL_DIR}/matrix
            fi
        else
            # Stop processing current taxonomic level
            println "\t[ERROR] Not enough reference genomes found\n\n"
        fi
    else
        # Current taxonomic level cannot be processed
        println "\t[ERROR] No reference genomes found\n\n"
    fi
}

# Define default value for --nproc
NPROC=1
# Minimum number of genomes per cluster
MIN_GENOMES=100
# Remove temporary data at the end of the pipeline
CLEANUP=false

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --cleanup)
            # Remove temporary data at the end of the pipeline
            CLEANUP=true
            ;;
        --db-dir=*)
            # Database directory
            DBDIR="${ARG#*=}"
            # Define helper
            if [[ "$DBDIR" =~ "?" ]]; then
                println "boundaries helper: --db-dir=directory\n\n"
                println "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "$DBDIR" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
            # Trim the last slash out of the path
            DBDIR="${DBDIR%/}"
            # Check whether the input directory exist
            if [[ ! -d $DBDIR ]]; then
                println "Input folder does not exist!\n"
                println "--db-dir=%s\n" "$DBDIR"
                exit 1
            fi
            ;;
        -h|--help)
            # Print extended help
            BOUNDARIES_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --kingdom=*)
            # Consider genomes whose lineage belong to a specific kingdom only
            KINGDOM="${ARG#*=}"
            # Define helper
            if [[ "$KINGDOM" =~ "?" ]]; then
                println "boundaries helper: --kingdom=value\n\n"
                println "\tSelect a kingdom.\n\n"
                exit 0
            fi
            ;;
        --log=*)
            # Path to the log file
            LOG_FILEPATH="${ARG#*=}"
            # Define helper
            if [[ "${LOG_FILEPATH}" =~ "?" ]]; then
                println "boundaries helper: --log=file\n\n"
                println "\tPath to the log file.\n\n"
                exit 0
            fi
            # Remove the log file if it already exists
            if [[ -f ${LOG_FILEPATH} ]]; then
                rm ${LOG_FILEPATH}
            fi
            ;;
        --min-genomes=*)
            # Consider clusters with at least this number of genomes
            MIN_GENOMES="${ARG#*=}"
            # Define helper
            if [[ "${MIN_GENOMES}" =~ "?" ]]; then
                println "boundaries helper: --min-genomes=num\n\n"
                println "\tConsider clusters with at least this number of genomes.\n\n"
                exit 0
            fi
            # Check whether --min-genomes is an integer
            if [[ ! ${MIN_GENOMES} =~ ^[0-9]+$ ]] || [[ "${MIN_GENOMES}" -le "1" ]]; then
                println "Argument --min-genomes must be a positive integer greater than 1\n"
                exit 1
            fi
            ;;
        --nproc=*)
            # Max nproc for all parallel instructions
            NPROC="${ARG#*=}"
            # Define helper
            if [[ "$NPROC" =~ "?" ]]; then
                println "boundaries helper: --nproc=num\n\n"
                println "\tThis argument refers to the number of processors used for parallelizing the pipeline when possible.\n"
                println "\tDefault: --nproc=1\n\n"
                exit 0
            fi
            # Check whether --proc is an integer
            if [[ ! $NPROC =~ ^[0-9]+$ ]] || [[ "$NPROC" -eq "0" ]]; then
                println "Argument --nproc must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --output=*)
            # Output file
            OUTPUT="${ARG#*=}"
            # Define helper
            if [[ "$OUTPUT" =~ "?" ]]; then
                println "boundaries helper: --output=file\n\n"
                println "\tOutput file with kmer boundaries for each of the taxonomic labels in the database.\n\n"
                exit 0
            fi
            ;;
        --tmp-dir=*)
            # Temporary folder
            TMPDIR="${ARG#*=}"
            # Define helper
            if [[ "$TMPDIR" =~ "?" ]]; then
                println "boundaries helper: --tmp-dir=directory\n\n"
                println "\tPath to the folder for storing temporary data.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            TMPDIR="$( cd "$( dirname "$TMPDIR" )" &> /dev/null && pwd )"/"$( basename $TMPDIR )"
            # Trim the last slash out of the path
            TMPDIR="${TMPDIR%/}"
            ;;
        -v|--version)
            # Print pipeline version
            println "boundaries version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            println "boundaries: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

println "boundaries version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false
if [[ "$?" -gt "0" ]]; then
    println "Unsatisfied software dependencies!\n\n"
    println "Please run the following command for a list of required external software dependencies:\n\n"
    println "\t$ meta-index --resolve-dependencies\n\n"
    
    exit 1
fi

# Create temporary folder
mkdir -p $TMPDIR

# Print inputs
println "Defining boundaries:\n"
println "\t--db-dir=%s\n" "$DBDIR"
println "\t--kingdom=%s\n" "$KINGDOM"
println "\t--min-genomes=%s\n\n" "${MIN_GENOMES}"

# Check whether the output file already exists
if [[ ! -f $OUTPUT ]]; then
    # Init output table
    println "# boundaries version %s (%s)\n" "$VERSION" "$DATE" > $OUTPUT
    println "# --db-dir=%s\n" "$DBDIR" >> $OUTPUT
    println "# --kingdom=%s\n" "$KINGDOM" >> $OUTPUT
    println "# --min-genomes=%s\n" "${MIN_GENOMES}" >> $OUTPUT
    println "# Lineage\tMin kmers\tMax kmers\n" >> $OUTPUT
    # Process all taxonomic levels
    for LEVEL in "s__" "g__" "f__" "o__" "c__" "p__" "k__"; do
        find $DBDIR/k__$KINGDOM -maxdepth $DEPTH -type d -iname "${LEVEL}*" -follow -exec \
            define_boundaries {} "${MIN_GENOMES}" "$TMPDIR" "$OUTPUT" "$NPROC" "$CLEANUP" \;
    done

    # Print output table path with boundaries
    println "\nOutput table:\n"
    println "\t%s\n" "$OUTPUT"

    # Cleanup temporary data
    if $CLEANUP; then
        # Remove tmp folder
        if [[ -d $TMPDIR ]]; then
            println "Cleaning up temporary folder:\n"
            println "\t%s\n" "$TMPDIR"
            rm -rf $TMPDIR
        fi
    fi
else
    # Output file already exists
    println "[ERROR] Output file already exists\n"
    println "\t%s\n" "$OUTPUT"

    exit 1
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0