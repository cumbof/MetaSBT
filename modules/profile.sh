#!/bin/bash
#title          :profile
#description    :Query the sequence bloom trees at all the 7 taxonomic levels
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#============================================================================

DATE="Jun 8, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Requirements base path
REQUIREMENTS_DIR="$(dirname "${SCRIPT_DIR}")/requirements"
# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Do not expand queries by default
# If true, keep querying trees at higher taxonomic levels
EXPAND_QUERY=false
# Theta threshold
THRESHOLD=0.0

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --expand)
            # Expand the input query to the higher taxonomic levels
            EXPAND_QUERY=true
            ;;
        -h|--help)
            # Print extended help
            PROFILE_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --input-file=*)
            # Input file path with queries
            INPUTFILE="${ARG#*=}"
            # Define helper
            if [[ "$INPUTFILE" =~ "?" ]]; then
                println "profile helper: --input-file=file\n\n"
                println "\tPath to the input genome.\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $INPUTFILE ]]; then
                println "Input file does not exist: --input-file=%s\n" "$INPUTFILE"
                exit 1
            fi
            ;;
        --input-id=*)
            # Input id
            INPUTID="${ARG#*=}"
            # Define helper
            if [[ "$INPUTID" =~ "?" ]]; then
                println "profile helper: --input-id=value\n\n"
                println "\tUnique identifier of the input genome.\n\n"
                exit 0
            fi
            ;;
        --log=*)
            # Path to the log file
            LOG_FILEPATH="${ARG#*=}"
            # Define helper
            if [[ "${LOG_FILEPATH}" =~ "?" ]]; then
                println "profile helper: --log=file\n\n"
                println "\tPath to the log file.\n\n"
                exit 0
            fi
            # Remove the log file if it already exists
            if [[ -f ${LOG_FILEPATH} ]]; then
                rm ${LOG_FILEPATH}
            fi
            ;;
        --output-dir=*)
            # Output folder
            OUTPUTDIR="${ARG#*=}"
            # Define helper
            if [[ "$OUTPUTDIR" =~ "?" ]]; then
                println "profile helper: --output-dir=directory\n\n"
                println "\tThis is the output folder with queries results.\n\n"
                exit 0
            fi
            ;;
        --output-prefix=*)
            # Output folder
            OUTPUTPREFIX="${ARG#*=}"
            # Define helper
            if [[ "$OUTPUTPREFIX" =~ "?" ]]; then
                println "profile helper: --output-prefix=value\n\n"
                println "\tPrefix of the output files with query matches.\n\n"
                exit 0
            fi
            ;;
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies true "${REQUIREMENTS_DIR}/profile.txt"
            exit $?
            ;;
        --stop-at=*)
            # Stop expanding queries at a specific taxonomic level
            STOPAT="${ARG#*=}"
            # Define helper
            TAXALEVELS_LIST=("phylum" "class" "order" "family" "genus")   # List of available taxa levels
            TAXALEVELS_STRING="$(printf ", %s" "${TAXALEVELS_LIST[@]}")"  # Concatenate values in list
            TAXALEVELS_STRING="${TAXALEVELS_STRING:2}"                    # Trim the first 2 characters out of the string
            if [[ "$STOPAT" =~ "?" ]]; then
                println "profile helper: --stop-at=value\n\n"
                println "\tStop expanding queries at a specific taxonomic level.\n"
                println "\tPlease note that this argument works in conjunction with --expand only.\n"
                println "\tAvailable values: %s\n\n" "${TAXALEVELS_STRING}"
                exit 0
            fi
            # Check whether the provided taxa level is in the list of available values
            if [[ ! " ${TAXALEVELS_LIST[*]} " =~ " $STOPAT " ]]; then
                println "Error!\n"
                println "\"%s\" is not a valid taxonomic level.\n" "$STOPAT"
                println "Please use one of these values: %s\n" "${TAXALEVELS_STRING}"
                exit 1
            fi
            ;;
        --threshold=*)
            # Theta threshold
            THRESHOLD="${ARG#*=}"
            # Define helper
            if [[ "$THRESHOLD" =~ "?" ]]; then
                println "profile helper: --threshold=number\n\n"
                println "\tFraction of query kmers that must be present in a leaf to be considered a match.\n"
                println "\tThis must be between 0.0 and 1.0.\n"
                println "\tDefault: 0.0\n\n"
                exit 0
            fi
            ;;
        --tree=*)
            # Tree definition file path
            TREE="${ARG#*=}"
            # Define helper
            if [[ "$TREE" =~ "?" ]]; then
                println "profile helper: --tree=file\n\n"
                println "\tThis is the tree definition file.\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $TREE ]]; then
                println "Input file does not exist: --tree=%s\n" "$TREE"
                exit 1
            fi
            ;;
        --verbose)
            # Print messages
            VERBOSE=true
        -v|--version)
            # Print pipeline version
            println "profile version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            println "profile: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

println "profile version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false "${REQUIREMENTS_DIR}/profile.txt"
if [[ "$?" -gt "0" ]]; then
    println "Unsatisfied software dependencies!\n\n"
    println "Please run the following command for a list of required external software dependencies:\n\n"
    println "\t$ meta-index --resolve-dependencies\n\n"

    exit 1
fi

println "Input:\n"
println "\tFile path: %s\n" "$INPUTFILE"
println "\tInput ID: %s\n\n" "$INPUTID"

# Create the output directory if it does not exist
if [[ ! -d $OUTPUTDIR ]]; then
    mkdir -p $OUTPUTDIR
fi

# Take track of the best matches
# Closest taxa
LINEAGE=()
# Common kmers with the closest taxa
COMMON_KMERS=()
# Total kmers in query sequence
TOTAL_KMERS=()
# Scores of the closest taxa
SCORES=()

# Closest genome
CLOSEST_GENOME=""
# Common kmers with the closest genome
CLOSEST_GENOME_COMMON_KMERS=0
# Total kmers in query sequence
CLOSEST_GENOME_TOTAL_KMERS=0
# Score of the closest genome
CLOSEST_GENOME_SCORE=0

println "Querying trees:\n"
# Start querying the tree
while [[ -f $TREE ]]; do
    IDXDIR="$(dirname $TREE)"
    LEVELDIR="$(dirname $IDXDIR)"
    LEVEL="$(basename $LEVELDIR)"
    # In case of species level
    if [[ "$(basename $IDXDIR)" = "howde_index" ]]; then
        LEVELDIR="$(dirname $LEVELDIR)"
        LEVEL="$(basename $LEVELDIR)"
    fi

    # Query the sequence bloom tree
    println "\t%s\n" "$TREE"
    howdesbt query --tree=$TREE \
                   --threshold=$THRESHOLD \
                   --sort \
                   $INPUTFILE > $OUTPUTDIR/${OUTPUTPREFIX}__${LEVEL}__matches.txt

    # Take track of matches in associative array
    declare -A MATCHES_COMMON_KMERS
    declare -A MATCHES_TOTAL_KMERS
    # Take track of the best match and its number of common kmers
    BEST_MATCH=""
    BEST_KMERS=0
    BEST_TOTAL_KMERS=0
    # Read the output file with matches
    while read line; do
        if [[ ! "$line" =~ ^\*.* ]] && [[ ! "$line" =~ ^#.* ]]; then
            # Get the current match and its number of common kmers
            MATCH="$(echo $line | cut -d' ' -f1)"
            MATCH_COMMON_KMERS="$(echo $line | cut -d' ' -f2 | cut -d'/' -f1)"
            MATCH_TOTAL_KMERS="$(echo $line | cut -d' ' -f2 | cut -d'/' -f2)"
            # Update the associative array
            MATCHES_COMMON_KMERS["$MATCH"]=$((${MATCHES_COMMON_KMERS["$MATCH"]} + ${MATCH_COMMON_KMERS}))
            MATCHES_TOTAL_KMERS["$MATCH"]=$((${MATCHES_TOTAL_KMERS["$MATCH"]} + ${MATCH_TOTAL_KMERS}))
            # Update the best match
            if [[ "${MATCHES_COMMON_KMERS["$MATCH"]}" -gt "${BEST_KMERS}" ]]; then
                BEST_MATCH="$MATCH"
                BEST_KMERS=${MATCHES_COMMON_KMERS["$MATCH"]}
                BEST_TOTAL_KMERS=${MATCHES_TOTAL_KMERS["$MATCH"]}
            fi
        fi
    done < $OUTPUTDIR/${OUTPUTPREFIX}__${LEVEL}__matches.txt

    # Define the best score as the number of common kmers on the total number of kmers in the query for the best match
    BEST_SCORE="$(bc -l <<< "${BEST_KMERS} / ${BEST_TOTAL_KMERS}")"
    BEST_SCORE="$(printf "%.2f" "${BEST_SCORE}")"

    if [[ "$LEVEL" = s__* ]]; then
        # There is nothing after the species tree
        CLOSEST_GENOME="${BEST_MATCH}"
        CLOSEST_GENOME_COMMON_KMERS="${BEST_KMERS}"
        CLOSEST_GENOME_TOTAL_KMERS="${BEST_TOTAL_KMERS}"
        CLOSEST_GENOME_SCORE="${BEST_SCORE}"
        break
    else
        # Update the lineage
        LINEAGE+=("${BEST_MATCH}")
        COMMON_KMERS+=("${BEST_KMERS}")
        TOTAL_KMERS+=("${BEST_TOTAL_KMERS}")
        SCORES+=("${BEST_SCORE}")

        # Keep querying if --expand
        TREE=$LEVELDIR/${BEST_MATCH}/index/index.detbrief.sbt
        if [[ "${BEST_MATCH}" = s__* ]]; then
            TREE=$LEVELDIR/${BEST_MATCH}/index/howde_index/index.detbrief.sbt
        fi
    
        # Stop querying trees
        if ! ${EXPAND_QUERY}; then
            # In case --expand is not set
            break
        elif ${EXPAND_QUERY} && [[ "${LEVEL:0:1}" = "${STOPAT:0:1}" ]]; then
            # In case --expand is set and 
            # the taxonomic level specified with --stop-at is reached
            # Compare the first character of the current taxonomic level and the one specified with --stop-at 
            break
        fi
    fi
done

# Define the output file with a mapping between the input genome name and the taxonomic characterisation
if [[ ! -f $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv ]]; then
    # Print header lines
    println "# Tree: %s" "$TREE" > $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
    println "# Input ID\tLevel\tTaxonomy\tCommon kmers\tScore\n" >> $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
fi

# Reconstruct the lineage
if [ ${#LINEAGE[@]} -gt 0 ]; then
    # Print closest lineage
    println "\nClosest lineage:\n"
    CLOSEST_LINEAGE=$(printf "|%s" "${LINEAGE[@]}")
    CLOSEST_LINEAGE=${CLOSEST_LINEAGE:1}
    println "\t%s\n\n" "${CLOSEST_LINEAGE}"
    
    # Retrieve the maximum length of the level strings
    MAX_LEN=0
    for lineage in "${LINEAGE[@]}"; do
        if [[ "${MAX_LEN}" < "${#lineage}" ]]; then
            MAX_LEN=${#lineage}
        fi
    done

    # Print scores
    println "Score:\n"
    for i in "${!LINEAGE[@]}"; do
        println "\t%-${MAX_LEN}s\t%s\n" "${LINEAGE[i]}" "${SCORES[i]}"
        # Retrieve taxonomic level
        LEVELID=${LINEAGE[i]:0:1}
        # Expand the level ID to the full level name
        case "$LEVELID" in
            "k") LEVELNAME="kingdom" ;;
            "p") LEVELNAME="phylum" ;;
            "c") LEVELNAME="class" ;;
            "o") LEVELNAME="order" ;;
            "f") LEVELNAME="family" ;;
            "g") LEVELNAME="genus" ;;
            "s") LEVELNAME="species" ;;
            *) LEVELNAME="NA" ;;
        esac
        # Report characterisation to the output file
        println "%s\t%s\t%s\t%s\t%s\n" "$INPUTID" \
                                       "$LEVELNAME" \
                                       "${LINEAGE[i]}" \
                                       "${COMMON_KMERS[i]}/${TOTAL_KMERS[i]}" \
                                       "${SCORES[i]}" \
                                       >> $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
    done
fi

# Show the closest genome in case the input is a species tree or the query is expanded
if [ ! -z "${CLOSEST_GENOME}" ]; then
    # Print the closest genome and its score
    println "\nClosest genome:\n"
    println "\t%s\t%s\n" "${CLOSEST_GENOME}" "${CLOSEST_GENOME_SCORE}"
    # Report the closest genome to the output file
    println "%s\t%s\t%s\t%s\t%s\n" "$INPUTID" \
                                   "genome" \
                                   "${CLOSEST_GENOME}" \
                                   "${CLOSEST_GENOME_COMMON_KMERS}/${CLOSEST_GENOME_TOTAL_KMERS}" \
                                   "${CLOSEST_GENOME_SCORE}" \
                                   >> $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0