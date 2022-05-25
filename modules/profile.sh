#!/bin/bash
#title          :profile
#description    :Query the sequence bloom trees at all the 7 taxonomic levels
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#============================================================================

DATE="May 24, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Do not expand queries by default
# If true, keep querying trees at higher taxonomic levels
EXPAND_QUERY=false
# Theta threshold
THRESHOLD=0.7

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
        --input=*)
            # Input file path with queries
            INPUT="${ARG#*=}"
            # Define helper
            if [[ "${INPUT}" =~ "?" ]]; then
                printf "profile helper: --input=file\n\n"
                printf "\tThis is the input file with the set of queries.\n\n"
                exit 0
            fi
            ;;
        --output-dir=*)
            # Output folder
            OUTPUTDIR="${ARG#*=}"
            # Define helper
            if [[ "${OUTPUTDIR}" =~ "?" ]]; then
                printf "profile helper: --output-dir=directory\n\n"
                printf "\tThis is the output folder with queries results.\n\n"
                exit 0
            fi
            ;;
        --output-prefix=*)
            # Output folder
            OUTPUTPREFIX="${ARG#*=}"
            # Define helper
            if [[ "${OUTPUTPREFIX}" =~ "?" ]]; then
                printf "profile helper: --output-prefix=value\n\n"
                printf "\tPrefix of the output files with query matches.\n\n"
                exit 0
            fi
            ;;
        --threshold=*)
            # Theta threshold
            THRESHOLD="${ARG#*=}"
            # Define helper
            if [[ "${THRESHOLD}" =~ "?" ]]; then
                printf "profile helper: --threshold=number\n\n"
                printf "\tFraction of query kmers that must be present in a leaf to be considered a match.\n"
                printf "\tThis must be between 0 and 1.\n"
                printf "\tDefault: 0.7\n\n"
                exit 0
            fi
            ;;
        --tree=*)
            # Tree definition file path
            TREE="${ARG#*=}"
            # Define helper
            if [[ "${TREE}" =~ "?" ]]; then
                printf "profile helper: --tree=file\n\n"
                printf "\tThis is the tree definition file.\n\n"
                exit 0
            fi
            ;;
        *)
            printf "profile: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

printf "profile version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false
if [[ "$?" -gt "0" ]]; then
    exit 1
fi

printf "Input:\n"
printf "\t%s\n\n" "${INPUT}"

# Create the output directory if it does not exist
if [[ ! -d ${OUTPUTDIR} ]]; then
    mkdir -p ${OUTPUTDIR}
fi

# Take track of the best matches
# Closest taxa
LINEAGE=()
# Scores of the closest taxa
SCORES=()
# Closest genome
CLOSEST_GENOME=""
# Score of the closest genome
CLOSEST_GENOME_SCORE=0

printf "Querying trees:\n"
# Start querying the tree
while [[ -f ${TREE} ]]; do
    IDXDIR="$(dirname $TREE)"
    LEVELDIR="$(dirname $IDXDIR)"
    LEVEL="$(basename $LEVELDIR)"
    # Query the sequence bloom tree
    printf "\t%s\n" "${TREE}"
    howdesbt query --tree=${TREE} \
                   --threshold=${THRESHOLD} \
                   --adjust \
                   --sort \
                   ${INPUT} > ${OUTPUTDIR}/${OUTPUTPREFIX}__${LEVEL}__matches.txt
    # Get best match
    BEST=$(grep -E "^[[:alnum:]]" ${OUTPUTDIR}/${OUTPUTPREFIX}__${LEVEL}__matches.txt | head -n 1)
    MATCH=$(echo "$BEST" | cut -d' ' -f1)
    SCORE=$(echo "$BEST" | cut -d' ' -f5)

    if [[ "$LEVEL" = s__* ]]; then
        # There is nothing after the species tree
        CLOSEST_GENOME=$MATCH
        CLOSEST_GENOME_SCORE=$SCORE
        break
    else
        # Update the lineage
        LINEAGE+=("$MATCH")
        SCORES+=("$SCORE")

        # Keep querying if --expand
        TREE=${LEVELDIR}/${MATCH}/index/index.detbrief.sbt
        if [[ "$MATCH" = s__* ]]; then
            TREE=${LEVELDIR}/${MATCH}/index/howde_index/index.detbrief.sbt
        fi

        # Stop querying trees
        if [[ "${EXPAND_QUERY}" = false ]]; then
            break
        fi
    fi
done

# Define the output file with a mapping between the input genome name and the taxonomic characterisation
if [[ ! -f ${OUTPUTDIR}/${OUTPUTPREFIX}__profiles.txt ]]; then
    printf "# Input\tLevel\tTaxonomy\tScore\n" > ${OUTPUTDIR}/${OUTPUTPREFIX}__profiles.txt
fi

# Reconstruct the lineage
if [ ${#LINEAGE[@]} -gt 0 ]; then
    # Print closest lineage
    printf "Closest lineage:\n"
    CLOSEST_LINEAGE=$(printf "|%s" "${LINEAGE[@]}")
    CLOSEST_LINEAGE=${CLOSEST_LINEAGE:1}
    printf "\t%s\n\n" "${CLOSEST_LINEAGE}"

    # Print scores
    printf "Score:\n"
    for i in "${!LINEAGE[@]}"; do
        printf "\t%s\t%s\n" "${LINEAGE[i]}" "${SCORES[i]}"
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
        printf "%s\t%s\t%s\t%s\n" "$INPUT" "$LEVELNAME" "${LINEAGE[i]}" "${SCORES[i]}" >> ${OUTPUTDIR}/${OUTPUTPREFIX}__profiles.txt
    done
    printf "\n"
fi

# Show the closest genome in case the input is a species tree or the query is expanded
if [ ! -z "${CLOSEST_GENOME}" ]; then
    # Print the closest genome and its score
    printf "Closest genome:\n"
    printf "\t%s\t%s\n" "${CLOSEST_GENOME}" "${CLOSEST_GENOME_SCORE}"
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0