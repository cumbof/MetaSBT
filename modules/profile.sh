#!/bin/bash
#title          :profile
#description    :Query the sequence bloom trees at all the 7 taxonomic levels
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#============================================================================

DATE="Jun 2, 2022"
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
                println "Input file does not exist!\n"
                println "--input-file=%s\n" "$INPUTFILE"
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
        --threshold=*)
            # Theta threshold
            THRESHOLD="${ARG#*=}"
            # Define helper
            if [[ "$THRESHOLD" =~ "?" ]]; then
                println "profile helper: --threshold=number\n\n"
                println "\tFraction of query kmers that must be present in a leaf to be considered a match.\n"
                println "\tThis must be between 0 and 1.\n"
                println "\tDefault: 0.7\n\n"
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
                println "Input file does not exist!\n"
                println "--tree=%s\n" "$TREE"
                exit 1
            fi
            ;;
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

check_dependencies false
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
# Scores of the closest taxa
SCORES=()

# Closest genome
CLOSEST_GENOME=""
# Common kmers with the closest genome
CLOSEST_GENOME_COMMON_KMERS=0
# Score of the closest genome
CLOSEST_GENOME_SCORE=0

println "Querying trees:\n"
# Start querying the tree
while [[ -f $TREE ]]; do
    IDXDIR="$(dirname $TREE)"
    LEVELDIR="$(dirname $IDXDIR)"
    LEVEL="$(basename $LEVELDIR)"
    # Query the sequence bloom tree
    println "\t%s\n" "$TREE"
    howdesbt query --tree=$TREE \
                   --threshold=$THRESHOLD \
                   --adjust \
                   --sort \
                   $INPUTFILE > $OUTPUTDIR/${OUTPUTPREFIX}__${LEVEL}__matches.txt
    # Get best match line
    BEST=$(grep -E "^[[:alnum:]]" $OUTPUTDIR/${OUTPUTPREFIX}__${LEVEL}__matches.txt | head -n 1)
    # Best match ID
    MATCH=$(echo "$BEST" | cut -d' ' -f1)
    # Common kmers with the best match
    KMERS=$(echo "$BEST" | cut -d' ' -f4 | cut -d'/' -f1)
    # Score
    SCORE=$(echo "$BEST" | cut -d' ' -f5)

    if [[ "$LEVEL" = s__* ]]; then
        # There is nothing after the species tree
        CLOSEST_GENOME=$MATCH
        CLOSEST_GENOME_COMMON_KMERS=$KMERS
        CLOSEST_GENOME_SCORE=$SCORE
        break
    else
        # Update the lineage
        LINEAGE+=("$MATCH")
        COMMON_KMERS+=("$KMERS")
        SCORES+=("$SCORE")

        # Keep querying if --expand
        TREE=$LEVELDIR/$MATCH/index/index.detbrief.sbt
        if [[ "$MATCH" = s__* ]]; then
            TREE=$LEVELDIR/$MATCH/index/howde_index/index.detbrief.sbt
        fi

        # Stop querying trees
        if [[ "${EXPAND_QUERY}" = false ]]; then
            break
        fi
    fi
done

# Define the output file with a mapping between the input genome name and the taxonomic characterisation
if [[ ! -f $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv ]]; then
    println "# Input ID\tLevel\tTaxonomy\tCommon kmers\tScore\n" > $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
fi

# Reconstruct the lineage
if [ ${#LINEAGE[@]} -gt 0 ]; then
    # Print closest lineage
    println "\nClosest lineage:\n"
    CLOSEST_LINEAGE=$(printf "|%s" "${LINEAGE[@]}")
    CLOSEST_LINEAGE=${CLOSEST_LINEAGE:1}
    println "\t%s\n\n" "${CLOSEST_LINEAGE}"

    # Print scores
    println "Score:\n"
    for i in "${!LINEAGE[@]}"; do
        println "\t%s\t%s\n" "${LINEAGE[i]}" "${SCORES[i]}"
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
        println "%s\t%s\t%s\t%s\n" "$INPUTID" \
                                   "$LEVELNAME" \
                                   "${LINEAGE[i]}" \
                                   "${COMMON_KMERS[i]}" \
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
    println "%s\t%s\t%s\t%s\n" "$INPUTID" \
                               "genome" \
                               "${CLOSEST_GENOME}" \
                               "${CLOSEST_GENOME_COMMON_KMERS}" \
                               "${CLOSEST_GENOME_SCORE}" \
                               >> $OUTPUTDIR/${OUTPUTPREFIX}__profiles.tsv
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0