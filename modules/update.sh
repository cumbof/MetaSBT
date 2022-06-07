#!/bin/bash
#title          :update
#description    :Update the index with new genomes
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="Jun 7, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Requirements base path
REQUIREMENTS_DIR="$(dirname "${SCRIPT_DIR}")/requirements"
# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Get the min and max kmers boundaries for a given taxonomy
# Result is reported in MIN_BOUND and MAX_BOUND
get_boundaries () {
    BOUNDARIES=$1   # Path to the output table of the boundaries module
    TAXLABEL=$2     # Current taxonomic label
    # Split taxonomy
    TAXONOMY_SPLIT=($(echo $TAXLABEL | tr "|" " "))

    # Keep track of the min and max common kmers
    MIN_BOUNDS=()
    MAX_BOUNDS=()
    # Try searching for boundaries again with a redefined taxonomic label
    RETRY=false
    while [[ "${#MIN_BOUNDS[@]}" -eq "0" ]] && [[ "${#MAX_BOUNDS[@]}" -eq "0" ]]; do
        if ! $RETRY; then
            # Add a tab at the end to avoid multiple matched at first
            TAX_BOUNDARIES="$(grep "$TAXLABEL"$'\t' $BOUNDARIES)"
        else
            # Add a pipe at the end of the taxonomic label for all the other levels
            TAX_BOUNDARIES="$(grep "$TAXLABEL|" $BOUNDARIES)"
        fi

        if [[ ! -z "${TAX_BOUNDARIES}" ]]; then
            # In case the current taxonomy is in the boundaries file
            while read line; do
                # Get minimum and maximum common kmers for the closest taxonomy
                MIN_BOUNDS+=("$(echo "${TAX_BOUNDARIES}" | cut -d$'\t' -f2)")
                MAX_BOUNDS+=("$(echo "${TAX_BOUNDARIES}" | cut -d$'\t' -f3)")
            done <${TAX_BOUNDARIES}
        else
            # Count the current number of levels in taxonomy
            TLEN=${#TAXONOMY_SPLIT[@]}
            if [[ "$TLEN" -eq "1" ]]; then
                # Get out of the while loop if there are no other levels available
                break
            fi
            # Remove the last level
            TAXONOMY_SPLIT=("${TAXONOMY_SPLIT[@]:0:$TLEN-1}")
            # Rebuild the taxonomic label
            TAXLABEL="$(printf "|%s" "${TAXONOMY_SPLIT[@]}")"
            # Retry with a redefined taxonomic label
            RETRY=true
        fi
    done

    # Concatenate array values with the plus symbol
    # Append a zero at the end of the string because of the last plus symbol
    MIN_BOUNDS_SUM=$( echo "${MIN_BOUNDS[@]/%/+} 0" | bc -l)
    # Count how many elements in the array
    MIN_BOUNDS_LEN=${#MIN_BOUNDS[@]}
    # Compute the mean value
    MIN_BOUND=$((${MIN_BOUNDS_SUM} / ${MIN_BOUNDS_LEN}))

    # Concatenate array values with the plus symbol
    # Append a zero at the end of the string because of the last plus symbol
    MAX_BOUNDS_SUM=$( echo "${MAX_BOUNDS[@]/%/+} 0" | bc -l)
    # Count how many elements in the array
    MAX_BOUNDS_LEN=${#MAX_BOUNDS[@]}
    # Compute the mean value
    MAX_BOUND=$((${MAX_BOUNDS_SUM} / ${MAX_BOUNDS_LEN}))
}

# Define default value for --nproc
NPROC=1
# Define default boundary uncertainty percentage
BOUNDARY_UNCERTAINTY_PERC=0
# Initialize CheckM completeness and contamination to default values
# Run CheckM if provided completeness is >0 or contamination is <100
CHECKM_COMPLETENESS=0
CHECKM_CONTAMINATION=100
# CheckM cannot handle a set of genomes with mixed file extensions
# Extensions of the input genomes must be standardized before running this script
# Set a default input genome extension in case it is not explicitly specified with the --extension argument
EXTENSION="fna.gz"
# Dereplicate input genomes
# Use kmtricks to build the kmer matrix on the input set of genomes
# Remove genomes with identical set of kmers
DEREPLICATE=false
# Remove temporary data at the end of the pipeline
CLEANUP=false

# Initially set the bloom filter size and the kmer length to 0
# In case --filter-size and --kmer-len are not specified, try to load these values from the database manifest file
# The manifest file is stored under the kingdom level directory and it is created by the index module
FILTER_SIZE=0
KMER_LEN=0

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --boundaries=*)
            # Clusters boundaries computed through the boundaries module
            BOUNDARIES="${ARG#*=}"
            # Define helper
            if [[ "$BOUNDARIES" =~ "?" ]]; then
                println "update helper: --boundaries=file\n\n"
                println "\tPath to the output table produced by the boundaries module.\n"
                println "\tIt is required in case of MAGs as input genomes only\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $BOUNDARIES ]]; then
                println "Input file does not exist: --boundaries=%s\n" "$BOUNDARIES"
                exit 1
            fi
            ;;
        --boundary-uncertainty=*)
            # Boundary uncertainty percentage
            BOUNDARY_UNCERTAINTY_PERC="${ARG#*=}"
            # Define helper
            if [[ "${BOUNDARY_UNCERTAINTY_PERC}" =~ "?" ]]; then
                println "update helper: --boundary-uncertainty=num\n\n"
                println "\tDefine the percentage of kmers to enlarge and reduce boundaries\n\n"
                exit 0
            fi
            # Check whether --boundary-uncertainty is an integer between 0 and 100
            if [[ ! ${BOUNDARY_UNCERTAINTY_PERC} =~ ^[0-9]+$ ]] || [[ "${BOUNDARY_UNCERTAINTY_PERC}" -lt "0" ]] || [[ "${BOUNDARY_UNCERTAINTY_PERC}" -gt "100" ]]; then
                println "Argument --boundary-uncertainty must be a positive integer between 0 and 100\n"
                exit 1
            fi
            ;;
        --checkm-completeness=*)
            # CheckM completeness
            CHECKM_COMPLETENESS="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_COMPLETENESS}" =~ "?" ]]; then
                println "update helper: --checkm-completeness=num\n\n"
                println "\tInput genomes must have a minimum completeness percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-completeness is an integer
            if [[ ! ${CHECKM_COMPLETENESS} =~ ^[0-9]+$ ]] || [[ "${CHECKM_COMPLETENESS}" -gt "100" ]]; then
                println "Argument --checkm-completeness must be a positive integer greater than 0 (up to 100)\n"
                exit 1
            fi
            ;;
        --checkm-contamination=*)
            # CheckM contamination
            CHECKM_CONTAMINATION="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_CONTAMINATION}" =~ "?" ]]; then
                println "update helper: --checkm-contamination=num\n\n"
                println "\tInput genomes must have a maximum contamination percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-contamination is an integer
            if [[ ! ${CHECKM_CONTAMINATION} =~ ^[0-9]+$ ]] || [[ "${CHECKM_CONTAMINATION}" -gt "100" ]]; then
                println "Argument --checkm-completeness must be a positive integer lower than 100 (up to 0)\n"
                exit 1
            fi
            ;;
        --cleanup)
            # Remove temporary data at the end of the pipeline
            CLEANUP=true
            ;;
        --db-dir=*)
            # Database directory
            DBDIR="${ARG#*=}"
            # Define helper
            if [[ "$DBDIR" =~ "?" ]]; then
                println "update helper: --db-dir=directory\n\n"
                println "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "$DBDIR" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
            # Trim the last slash out of the path
            DBDIR="${DBDIR%/}"
            # Check whether the input directory exist
            if [[ ! -d $DBDIR ]]; then
                println "Input folder does not exist: --db-dir=%s\n" "$DBDIR"
                exit 1
            fi
            ;;
        --dereplicate)
            # Dereplicate input genomes
            DEREPLICATE=true
            ;;
        --extension=*)
            # Input genome files extension
            # It must be standardised before running this tool module
            # All the input genomes must have the same file extension
            EXTENSION="${ARG#*=}"
            # Define helper
            if [[ "$EXTENSION" =~ "?" ]]; then
                println "update helper: --extension=value\n\n"
                println "\tSpecify the input genome files extension.\n"
                println "\tAll the input genomes must have the same file extension before running this module.\n"
                println "\tPlease note that this module supports a limited set of file extensions for the input genomes.\n"
                println "\tIn order to list all the supported extensions, please run this module with the --supported-extensions argument.\n\n"
                exit 0
            fi
            # Check whether the specified extension is supported by querying the supported_inputs.txt file
            QUERYEXT="$EXTENSION"
            if [[ "$QUERYEXT" = *.gz ]]; then
                QUERYEXT="${QUERYEXT%.*}"
            fi
            # Search for the extension
            if ! grep -q "^$QUERYEXT"$ $SCRIPTSDIR/supported_inputs.txt ; then
                println "Error!\n"
                println "\"%s\" is not a valid extension.\n" "$QUERYEXT"
                println "Please use the --supported-extensions option for a list of supported file extensions for input genomes.\n"
                exit 1
            fi
            ;;
        --filter-size=*)
            # Bloom filter size
            FILTER_SIZE="${ARG#*=}"
            # Define helper
            if [[ "${FILTER_SIZE}" =~ "?" ]]; then
                println "update helper: --filter-size=num\n\n"
                println "\tThis is the size of the bloom filters.\n"
                println "\tIt must be the same size used while running the index module.\n\n"
                exit 0
            fi
            # Check whether --filter-size is an integer
            if [[ ! ${FILTER_SIZE} =~ ^[0-9]+$ ]] || [[ "${FILTER_SIZE}" -eq "0" ]]; then
                println "Argument --filter-size must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        -h|--help)
            # Print extended help
            UPDATE_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --input-list=*)
            # Input file with the list of paths to the new genomes
            INLIST="${ARG#*=}"
            # Define helper
            if [[ "$INLIST" =~ "?" ]]; then
                println "update helper: --input-list=file\n\n"
                println "\tThis file contains the list of paths to the new genomes that will be added to the database.\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $INLIST ]]; then
                println "Input file does not exist: --input-list=%s\n" "$INLIST"
                exit 1
            fi
            ;;
        --kingdom=*)
            # Consider genomes whose lineage belong to a specific kingdom only
            KINGDOM="${ARG#*=}"
            # Define helper
            if [[ "$KINGDOM" =~ "?" ]]; then
                println "update helper: --kingdom=value\n\n"
                println "\tSelect the kingdom that will be updated.\n\n"
                exit 0
            fi
            ;;
        --kmer-len=*)
            # Length of the kmers
            KMER_LEN="${ARG#*=}"
            # Define helper
            if [[ "${KMER_LEN}" =~ "?" ]]; then
                println "update helper: --kmer-len=num\n\n"
                println "\tThis is the length of the kmers used for building bloom filters.\n\n"
                exit 0
            fi
            # Check whether --kmer-len is an integer
            if [[ ! ${KMER_LEN} =~ ^[0-9]+$ ]] || [[ "${KMER_LEN}" -eq "0" ]]; then
                println "Argument --kmer-len must be a positive integer greater than 0\n"
                exit 1
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
        --nproc=*)
            # Max nproc for all parallel instructions
            NPROC="${ARG#*=}"
            # Define helper
            if [[ "$NPROC" =~ "?" ]]; then
                println "update helper: --nproc=num\n\n"
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
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies true "${REQUIREMENTS_DIR}/update.txt"
            exit $?
            ;;
        --supported-extensions)
            # Print the supported list of file extensions for the input sequences
            println 'List of supported file extensions for input genomes:\n'
            while IFS= read -r line || [[ -n "$line" ]]; do
                # Exclude comment lines
                if [[ ! "$line" =~ ^#.* ]]; then
                    println '\t%s\n' "$line"
                fi
            done < ${SCRIPT_DIR}/supported_inputs.txt
            exit 0
            ;;
        --taxa=*)
            # Input file with the mapping between input reference genome IDs and their taxonomic label
            TAXA="${ARG#*=}"
            # Define helper
            if [[ "$TAXA" =~ "?" ]]; then
                println "update helper: --taxa=file\n\n"
                println "\tInput file with the mapping between input genome IDs and their taxonomic label.\n"
                println "\tThis is used in case of reference genomes only \"--type=references\".\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $TAXA ]]; then
                println "Input file does not exist: --taxa=%s\n" "$TAXA"
                exit 1
            fi
            ;;
        --tmp-dir=*)
            # Temporary folder
            TMPDIR="${ARG#*=}"
            # Define helper
            if [[ "$TMPDIR" =~ "?" ]]; then
                println "update helper: --tmp-dir=directory\n\n"
                println "\tPath to the folder for storing temporary data.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            TMPDIR="$( cd "$( dirname "$TMPDIR" )" &> /dev/null && pwd )"/"$( basename $TMPDIR )"
            # Trim the last slash out of the path
            TMPDIR="${TMPDIR%/}"
            ;;
        --type=*)
            # Input genomes type: references or MAGs
            TYPE="${ARG#*=}"
            # Define helper
            if [[ "${TYPE}" =~ "?" ]]; then
                println "update helper: --type=value\n\n"
                println "\tDefine the nature of the input genomes. Only \"references\" and \"MAGs\" are allowed.\n\n"
                exit 0
            fi
            # Allowed genome types: "references" and "MAGs"
            if [[ ! "$TYPE" = "MAGs" ]] && [[ ! "$TYPE" = "references" ]]; then
                println "Input type \"%s\" is not allowed!\n" "$TYPE"
                println "Please have a look at the helper of the --type argument for a list of valid genome types.\n"
                exit 1
            fi
            ;;
        -v|--version)
            # Print pipeline version
            println "update version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        *)
            println "update: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

println "update version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false "${REQUIREMENTS_DIR}/update.txt"
if [[ "$?" -gt "0" ]]; then
    println "Unsatisfied software dependencies!\n\n"
    println "Please run the following command for a list of required external software dependencies:\n\n"
    println "\t$ meta-index --resolve-dependencies\n\n"
    
    exit 1
fi

# Retrieve the bloom filter size and the kmer length used for building the database from the manifest file
if [[ -f $DBDIR/k__$KINGDOM/menifest.txt ]]; then
    # Take track of the --filter-size and --kmer-len
    DB_FILTER_SIZE=0
    DB_KMER_LEN=0
    # Load data from the manifest file
    while read line; do
        case "$line" in
            --filter-size=*)
                # Retrieve the bloom filter size
                DB_FILTER_SIZE="${line#*=}"
                ;;
            --kmer-len=*)
                # Retrieve the kmer length
                DB_KMER_LEN="${line#*=}"
                ;;
            *)
                continue
                ;;
        esac
    done < $DBDIR/k__$KINGDOM/menifest.txt

    # Check whether --filter-size and --kmer-len were correctly retrieved from the manifest
    if [[ "${DB_FILTER_SIZE}" -eq "0" ]] || [[ "${DB_KMER_LEN}" -eq "0" ]]; then
        println "[ERROR] Data not found!\n"
        println "%s\n\n" "$DBDIR/k__$KINGDOM/menifest.txt"
        println "Unable to retrieve the bloom filter size and the kmer length used for building the database\n"
        exit 1
    fi
else
    println "[ERROR] File not found!\n"
    println "%s\n\n" "$DBDIR/k__$KINGDOM/menifest.txt"
    println "Unable to retrieve the bloom filter size and the kmer length used for building the database\n"
    exit 1
fi

# Check whether the specified --filter-size matches the one reported in the manifest file
if [[ ! "${FILTER_SIZE}" -eq "${DB_FILTER_SIZE}" ]]; then
    if [[ "${FILTER_SIZE}" -eq "0" ]]; then
        # In case the --filter-size argument is not specified
        # Use the bloom filter size reported in the manifest file of the database
        FILTER_SIZE=${DB_FILTER_SIZE}
    else
        # Otherwise, raise an error message and exit
        println "[ERROR] The specified --filter-size does not match with the size of the bloom filters in the database!\n"
        exit 1
    fi
fi

# Also check whether the specified --kmer-len matches the one reported in the manifest file
if [[ ! "${KMER_LEN}" -eq "${DB_KMER_LEN}" ]]; then
    if [[ "${KMER_LEN}" -eq "0" ]]; then
        # In case the --kmer-len argument is not specified
        # Use the kmer length reported in the manifest file of the database
        KMER_LEN=${DB_KMER_LEN}
    else
        # Otherwise, raise an error message and exit
        println "[ERROR] The specified --kmer-len does not match with the kmer length used for building the database!\n"
        exit 1
    fi
fi

# Create temporary folder
mkdir -p $TMPDIR

# Count how many genomes in input
HOW_MANY=$(cat $INLIST | wc -l)
println "Processing %s genomes\n" "${HOW_MANY}"
println "\t%s\n\n" "$INLIST"

# Input genomes must be quality-controlled before being added to the database
if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
    CHECKMTABLES="" # This is the result of the "run_checkm" function as a list of file paths separated by comma
    run_checkm $INLIST "update" $EXTENSION $TMPDIR $NPROC

    # Read all the CheckM output tables and filter genomes on their completeness and contamination scores 
    # according to the those provided in input
    touch $TMPDIR/genomes_qc.txt
    CHECKMTABLES_ARR=($(echo $CHECKMTABLES | tr "," " "))
    # In case there is at least one CheckM output table
    if [[ "${#CHECKMTABLES_ARR[@]}" -gt "0" ]]; then
        # Merge all tables
        # Retrieve header from the first file in list
        echo "$(head -n1 ${CHECKMTABLES_ARR[0]})" > $TMPDIR/checkm.tsv
        for table in ${CHECKMTABLES_ARR[@]}; do
            # Current table may not exist in case CheckM stopped working
            if [[ -f $table ]]; then
                # Get all the lines except the first one
                echo "$(tail -n +2 $table)" >> $TMPDIR/checkm.tsv
            fi
        done

        # Check whether the checkm.tsv table has been created
        if [[ -f $TAXDIR/checkm.tsv ]]; then
            # Iterate over input genomes
            for GENOMEPATH in $(cut -d" " -f1 "$INLIST"); do
                # Retrive the genome ID from the genome file path
                GENOMENAME="$(basename $GENOMEPATH)"
                GENOMEID="${GENOMENAME%".$EXTENSION"}"
                # Search for current genome ID into the CheckM output table
                GENOME_DATA="$(grep "$GENOMEID"$'\t' $TMPDIR/checkm.tsv)"
                if [[ ! -z "${GENOME_DATA}" ]]; then
                    # In case the current genome has been processed
                    # Retrieve completeness and contamination of current genome
                    COMPLETENESS="$(echo ${GENOME_DATA} | rev | cut -d' ' -f3 | rev)"  # Get completeness
                    CONTAMINATION="$(echo ${GENOME_DATA} | rev | cut -d' ' -f3 | rev)" # Get contamination
                    # Use bc command for comparing floating point numbers
                    COMPLETENESS_CHECK="$(echo "$COMPLETENESS >= ${CHECKM_COMPLETENESS}" | bc)"
                    CONTAMINATION_CHECK="$(echo "$CONTAMINATION <= ${CHECKM_CONTAMINATION}" | bc)"
                    if [[ "${COMPLETENESS_CHECK}" -eq "1" ]] && [[ "${CONTAMINATION_CHECK}" -eq "1" ]]; then
                        # Current genome passed the quality control
                        grep -w "${GENOMEID}.${EXTENSION}" $INLIST >> $TMPDIR/genomes_qc.txt
                    fi
                fi
            done
        fi
    fi

    # Create a new INLIST file with the new list of genomes
    INLIST=$TMPDIR/genomes_qc.txt
    # Count how many genomes passed the quality control
    HOW_MANY=$(cat $INLIST | wc -l)

    println "\t%s genomes passed the quality control\n" "${HOW_MANY}"
    println "\t\tMinimum completeness: %s\n" "${CHECKM_COMPLETENESS}"
    println "\t\tMaximum contamination: %s\n" "${CHECKM_CONTAMINATION}"
fi

# Use kmtricks to build a kmer matrix and compare input genomes
# Discard genomes with the same set of kmers (dereplication)
if $DEREPLICATE; then
    # Count how many input genomes survived in case they have been quality controlled
    HOW_MANY=$(cat $INLIST | wc -l)
    if [[ "${HOW_MANY}" -gt "1" ]]; then        
        GENOMES_FOF=$TMPDIR/genomes.fof
        # Build a fof file with the list of input genomes
        while read -r genomepath; do
            GENOME_NAME="$(basename $genomepath)"
            echo "${GENOME_NAME} : $genomepath" >> ${GENOMES_FOF}
        done <$INLIST

        if [[ -f ${GENOMES_FOF} ]]; then
            println "\nDereplicating %s input genomes\n" "${HOW_MANY}"
            # Run kmtricks to build the kmers matrix
            kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                    $TMPDIR/matrix \
                                    $NPROC \
                                    $TMPDIR/kmers_matrix.txt \
                                    $TMPDIR
            
            # Add header to the kmers matrix
            HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
            if [[ "${HEADER: -1}" = " " ]]; then
                # Trim the last character out of the header in xase of space
                HEADER="${HEADER:0:-1}"
            fi
            echo "#kmer $HEADER" > $TMPDIR/kmers_matrix_wh.txt
            cat $TMPDIR/kmers_matrix.txt >> $TMPDIR/kmers_matrix_wh.txt
            mv $TMPDIR/kmers_matrix_wh.txt $TMPDIR/kmers_matrix.txt

            # Dereplicate genomes
            python3 ${SCRIPT_DIR}/utils.py filter_genomes $TMPDIR/kmers_matrix.txt $TMPDIR/filtered.txt > /dev/null

            if [[ -f $TMPDIR/filtered.txt ]]; then
                # Dereplicate genomes
                FILTERED_ARRAY=()
                while read -r genome; do
                    if ! grep -q "^$genome$" $TMPDIR/filtered.txt; then
                        # Rebuild the absolute path to the genome file
                        realpath -s $(grep "^$genome " $TMPDIR/genomes.fof | cut -d' ' -f3) >> $TMPDIR/genomes_dereplicated.txt
                    fi
                done <<<"$(head -n1 $TMPDIR/kmers_matrix.txt | tr " " "\n")"

                if [[ -f $TMPDIR/genomes_dereplicated.txt ]]; then
                    # Build the new kmers matrix with the dereplicated genomes
                    FILTERED="$(printf ",%s" "${FILTERED_ARRAY[@]}")"
                    FILTERED="${FILTERED:1}"
                    csvcut --delimiter " " --not-columns $FILTERED $TMPDIR/kmers_matrix.txt > $TMPDIR/kmers_matrix_dereplicated.txt
                    # Use the new list of dereplicated genomes
                    INLIST=$TMPDIR/genomes_dereplicated.txt
                fi
            fi

            # Count how many genomes passed the dereplication
            HOW_MANY=$(cat $INLIST | wc -l)
            println "\t%s genomes passed the dereplication\n" "${HOW_MANY}"
        fi
    else
        # No input genomes survived from the quality control step
        println "No input genomes available!\n"
        PIPELINE_END_TIME="$(date +%s.%3N)"
        PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
        println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"
        
        # Print credits before stopping the pipeline execution
        credits
        
        exit 0
    fi
fi

# Count how many input genomes survived in case they have been quality controlled and dereplicated
HOW_MANY=$(cat $INLIST | wc -l)
if [[ "${HOW_MANY}" -gt "1" ]]; then
    # Create a temporary folder in case input genomes are gzip compressed
    mkdir -p $TMPDIR/genomes

    # In case one or more reference genomes are assigned to an unknown cluster
    # Keep track of the unknown taxonomy and change the lineage with 
    # the most occurring taxa among the assinged reference genomes
    declare -A TO_KNOWN_TAXA
    # Keep track of the lineages that must be rebuilt
    REBUILD=()
    # Keep track of the unassigned genomes
    # They will be assigned to new group
    UNASSIGNED=()
    # Process input genomes
    while read GENOMEPATH; do
        GENOMEEXT="${GENOMEPATH##*.}"
        GENOMENAME="$(basename ${GENOMEPATH%.*})"
        FILEPATH=$GENOMEPATH
        # Decompress input genome in case it is gzip compressed
        if [[ "$GENOMEPATH" = *.gz ]]; then
            GENOMEEXT="${GENOMENAME##*.}"
            GENOMENAME="${GENOMENAME%.*}"
            gunzip -c $GENOMEPATH > $TMPDIR/genomes/${GENOMENAME}.fna
            FILEPATH=$TMPDIR/genomes/${GENOMENAME}.fna
        fi

        println "\nProfiling %s\n" "$GENOMENAME"
        # Run the profiler to establish the closest genome and the closest group for each taxonomic level in the tree
        . ${SCRIPT_DIR}/profile.sh --input-file=$FILEPATH \
                                   --input-id=$FILEPATH \
                                   --tree=$DBDIR/k__$KINGDOM/index/index.detbrief.sbt
                                   --expand \
                                   --output-dir=$TMPDIR/profiling/ \
                                   --output-prefix=$GENOMENAME \
                                   > /dev/null 2>&1 # Silence the profiler
        # Define output file path with profiles
        PROFILE=$TMPDIR/profiling/${GENOMENAME}__profiles.tsv
        println "\t%s\n" "$PROFILE"
        
        if [[ -f $PROFILE ]]; then
            # Discard the genome if the score is 1.0 compared to the closest match at the species level
            # Do not discard the input genome if it is a reference genome and the closest genome in the database
            # is a MAG with a perfect overlap of kmers (score 1.0)
            CLOSEST_GENOME_DATA="$(grep "${GENOMENAME}.${GENOMEEXT}" $PROFILE | grep -w "genome")"
            CLOSEST_GENOME="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f3)"
            CLOSEST_GENOME_COMMON_KMERS="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f4)"
            CLOSEST_GENOME_SCORE="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f5)"
            # Print closest genome
            println "\tClosest genome: %s (common kmers: %s; score: %s)\n" "${CLOSEST_GENOME}" "${CLOSEST_GENOME_COMMON_KMERS}" "${CLOSEST_GENOME_SCORE}"
            # Reconstruct the closest taxa
            LEVELS=("kingdom" "phylum" "class" "order" "family" "genus" "species")
            CLOSEST_TAXA=""
            CLOSEST_COMMON_KMERS=""
            println "\tClosest lineage:\n"
            for level in ${LEVELS[@]}; do
                CLOSEST_LEVEL_DATA="$(grep "${GENOMENAME}.${GENOMEEXT}" $PROFILE | grep -w "$level")"
                CLOSEST_LEVEL="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f3)"
                CLOSEST_LEVEL_COMMON_KMERS="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f4)"
                CLOSEST_LEVEL_SCORE="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f5)"
                CLOSEST_TAXA=${CLOSEST_TAXA}"|"${CLOSEST_LEVEL}
                # Print profiler result
                println "\t\t%s: %s (common kmers: %s; score: %s)\n" "$level" "${CLOSEST_LEVEL}" "${CLOSEST_LEVEL_COMMON_KMERS}" "${CLOSEST_LEVEL_SCORE}"
                # Keep track of the species score
                if [[ "${CLOSEST_LEVEL}" = s__* ]]; then
                    CLOSEST_COMMON_KMERS=${CLOSEST_LEVEL_COMMON_KMERS}
                fi
            done
            # Trim the first pipe out of the closest taxonomic label
            CLOSEST_TAXA=${CLOSEST_TAXA:1}
            # Define taxa path
            CLOSEST_TAXADIR=$DBDIR/$(echo "${CLOSEST_TAXA}" | sed 's/|/\//g')

            SKIP_GENOME=false
            # Check whether the input genome must be discarded
            if [[ "${CLOSEST_GENOME_SCORE}" -eq "1.0" ]]; then
                if [[ "$TYPE" = "MAGs" ]]; then
                    if [[ -f "${CLOSEST_TAXADIR}/genomes.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/genomes.txt; then
                        # If the input genome is a MAG and the closest genome is a reference
                        # Discard the input genome
                        SKIP_GENOME=true
                        # Print the reason why the input genome has been discarded
                        println "\tDiscarding genome:\n"
                        println "\t\tInput genome is a MAG and the closest genome is a reference genome\n"
                    elif [[ -f "${CLOSEST_TAXADIR}/mags.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/mags.txt; then
                        # If the input genome is a MAG and the closest genome is a MAG
                        # Discard the input genome
                        SKIP_GENOME=true
                        # Print the reason why the input genome has been discarded
                        println "\tDiscarding genome:\n"
                        println "\t\tInput genome and the closest genome are both MAGs\n"
                    fi
                elif [[ "$TYPE" = "references" ]]; then
                    if [[ -f "${CLOSEST_TAXADIR}/genomes.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/genomes.txt; then
                        # If the input genome is a reference and the closest genome is a reference
                        # Discard the input genome
                        SKIP_GENOME=true
                        # Print the reason why the input genome has been discarded
                        println "\tDiscarding genome:\n"
                        println "\t\tInput genome and closest genome are both reference genomes\n"
                    fi
                fi
            fi

            if ! ${SKIP_GENOME}; then
                # Retrieve the min and max common kmers for the closest taxa
                MIN_BOUND=""; MAX_BOUND=""
                get_boundaries $BOUNDARIES ${CLOSEST_TAXA}

                # Process MAGs or reference genomes
                if [[ "$TYPE" = "MAGs" ]]; then
                    if [[ "${CLOSEST_COMMON_KMERS}" -le "${MAX_BOUND}" ]] && [[ "${CLOSEST_COMMON_KMERS}" -ge "${MIN_BOUND}" ]]; then
                        # Assign the current genome to the closest lineage
                        TARGETGENOME=${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                        cp $FILEPATH $TARGETGENOME
                        if [[ ! "$GENOMEEXT" = "gz" ]]; then
                            # Compress input genome
                            gzip ${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                            TARGETGENOME=${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}.gz
                        fi
                        # Add the current genome to the list of genomes in current taxonomy
                        echo "$GENOMENAME : $TARGETGENOME" >> ${CLOSEST_TAXADIR}/genomes.fof
                        echo "$TARGETGENOME" >> ${CLOSEST_TAXADIR}/mags.txt
                        # Also add it CheckM statistics if available
                        if [[ -f $TMPDIR/checkm.tsv ]]; then
                            if [[ ! -f ${CLOSEST_TAXADIR}/checkm.tsv ]]; then
                                # In case the CheckM table does not exist in the closest taxa folder
                                # Create the CheckM table and add the header line
                                echo "$(head -n1 $TMPDIR/checkm.tsv)" > ${CLOSEST_TAXADIR}/checkm.tsv
                            fi
                            # Check whether the current genome is in the CheckM table
                            CHECKM_DATA="$(grep "$GENOMENAME"$'\t' $TMPDIR/checkm.tsv)"
                            if [[ ! -z "${CHECKM_DATA}" ]]; then
                                # Append the CheckM statistics for the current genome in the CheckM table of the closest table
                                echo "${CHECKM_DATA}" >> ${CLOSEST_TAXADIR}/checkm.tsv
                            fi
                        fi
                        # Do not remove the index here because there could be other input genomes with the same taxonomic label
                        REBUILD+=(${CLOSEST_TAXA})
                        # Print assignment message
                        println "\tAssignment:\n"
                        println "\t\t%s\n\n" "${CLOSEST_TAXA}"
                    else
                        # Mark the genome as unassigned
                        UNASSIGNED+=($GENOMEPATH)
                        # The CheckM statistics for this genome will be reported after the assignment
                        # Print unassignment message
                        println "\tUnassigned\n\n"
                    fi
                elif [[ "$TYPE" = "references" ]]; then
                    # Retrieve the taxonomic label from --taxa input file
                    TAXALABEL="$(grep -w "$GENOMENAME" $TAXA | cut -d$'\t' -f2)"
                    # Define taxonomy folder
                    TAXDIR=$DBDIR/$(echo "$TAXALABEL" | sed 's/|/\//g')
                    TARGETGENOME=$TAXDIR/genomes/${GENOMENAME}.${GENOMEEXT}
                    # Print input genome taxonomic label
                    println "\tTaxonomic label:\n"
                    println "\t\t%s\n" "$TAXALABEL"

                    if [[ "${CLOSEST_COMMON_KMERS}" -le "${MAX_BOUND}" ]] && [[ "${CLOSEST_COMMON_KMERS}" -ge "${MIN_BOUND}" ]]; then
                        # Check whether the closest taxonomy contains any reference genomes
                        HOW_MANY_REFERENCES=$(cat ${CLOSEST_TAXADIR}/genomes.txt | wc -l)
                        if [[ "${HOW_MANY_REFERENCES}" -eq "0" ]]; then
                            # If the closest genome belongs to a new cluster with no reference genomes
                            # Assign the current reference genome to the new cluster and rename its lineage with the taxonomic label of the reference genome
                            # Do not change lineage here because there could be more than one reference genome assigned to the same unknown cluster
                            # Assign a new taxonomy according to a majority rule applied on the reference genomes taxa
                            if [[ ! -z "${TO_KNOWN_TAXA["${CLOSEST_TAXA}"]}" ]]; then
                                TO_KNOWN_TAXA["${CLOSEST_TAXA}"]="$GENOMENAME"
                            else
                                # If already exists in the associative array
                                # Concatenate values with a comma
                                TO_KNOWN_TAXA["${CLOSEST_TAXA}"]=${TO_KNOWN_TAXA["${CLOSEST_TAXA}"]}",$GENOMENAME"
                            fi
                            # Print the assignment
                            # Taxonomy will change according to the majority voting result
                            println "\tAssignment:\n"
                            println "\t\t%s\n\n" "${CLOSEST_TAXA}"
                        elif [[ "${HOW_MANY_REFERENCES}" -ge "1" ]]; then
                            # If the closest genome belong to a cluster with at least a reference genome
                            if [[ "$TAXALABEL" != "${CLOSEST_TAXA}" ]]; then
                                # If the taxonomic labels of the current reference genome and that one of the closest genome do not match
                                # Report the inconsistency
                                println "\tInconsistency found:\n"
                                println "\t\tInput genome: %s\n" "$TAXALABEL"
                                println "\t\tClosest lineage: %s\n" "${CLOSEST_TAXA}"
                            fi
                            # Assign the current genome to the closest lineage
                            TARGETGENOME=${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                            cp $FILEPATH $TARGETGENOME
                            if [[ ! "$GENOMEEXT" = "gz" ]]; then
                                # Compress input genome
                                gzip ${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                                TARGETGENOME=${CLOSEST_TAXADIR}/genomes/${GENOMENAME}.${GENOMEEXT}.gz
                            fi
                            # Add the current genome to the list of genomes in current taxonomy
                            echo "$GENOMENAME : $TARGETGENOME" >> ${CLOSEST_TAXADIR}/genomes.fof
                            echo "$TARGETGENOME" >> ${CLOSEST_TAXADIR}/genomes.txt
                            # Do not remove the index here because there could be other input genomes with the same taxonomic label
                            REBUILD+=(${CLOSEST_TAXA})
                            # Assign the current reference genome to the closest cluster
                            println "\tAssignment:\n"
                            println "\t\t%s\n\n" "${CLOSEST_TAXA}"
                        fi
                        
                        # Also add it CheckM statistics if available
                        if [[ -f $TMPDIR/checkm.tsv ]]; then
                            if [[ ! -f ${CLOSEST_TAXADIR}/checkm.tsv ]]; then
                                # In case the CheckM table does not exist in the closest taxa folder
                                # Create the CheckM table and add the header line
                                echo "$(head -n1 $TMPDIR/checkm.tsv)" > ${CLOSEST_TAXADIR}/checkm.tsv
                            fi
                            # Check whether the current genome is in the CheckM table
                            CHECKM_DATA="$(grep "$GENOMENAME"$'\t' $TMPDIR/checkm.tsv)"
                            if [[ ! -z "${CHECKM_DATA}" ]]; then
                                # Append the CheckM statistics for the current genome in the CheckM table of the closest table
                                echo "${CHECKM_DATA}" >> ${CLOSEST_TAXADIR}/checkm.tsv
                            fi
                        fi
                    else
                        # If nothing is close enough to the current genome and its taxonomic label does not exist in the database
                        # Create a new branch with the new taxonomy for the current genome
                        mkdir -p $TAXDIR/genomes
                        cp $FILEPATH $TARGETGENOME
                        if [[ ! "$GENOMEEXT" = "gz" ]]; then
                            # Compress input genome
                            gzip ${TAXDIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                            TARGETGENOME=${TAXDIR}/genomes/${GENOMENAME}.${GENOMEEXT}.gz
                        fi
                        # Add the current genome to the list of genomes in current taxonomy
                        echo "$GENOMENAME : $TARGETGENOME" >> $TAXDIR/genomes.fof
                        echo "$TARGETGENOME" >> $TAXDIR/genomes.txt
                        # Do not remove the index here because there could be other input genomes with the same taxonomic label
                        REBUILD+=($TAXALABEL)

                        # Also add it CheckM statistics if available
                        if [[ -f $TMPDIR/checkm.tsv ]]; then
                            if [[ ! -f $TAXDIR/checkm.tsv ]]; then
                                # In case the CheckM table does not exist in the closest taxa folder
                                # Create the CheckM table and add the header line
                                echo "$(head -n1 $TMPDIR/checkm.tsv)" > $TAXDIR/checkm.tsv
                            fi
                            # Check whether the current genome is in the CheckM table
                            CHECKM_DATA="$(grep "$GENOMENAME"$'\t' $TMPDIR/checkm.tsv)"
                            if [[ ! -z "${CHECKM_DATA}" ]]; then
                                # Append the CheckM statistics for the current genome in the CheckM table of the closest table
                                echo "${CHECKM_DATA}" >> $TAXDIR/checkm.tsv
                            fi
                        fi
                    fi
                fi
            fi
        fi

        # Remove the uncompressed genome in temporary folder
        if [[ -f $TMPDIR/genomes/${GENOMENAME}.fna ]]; then 
            rm $TMPDIR/genomes/${GENOMENAME}.fna
        fi
    done <$INLIST

    # In case a reference has been assigned to an unknown cluster
    # Update the cluster taxonomy by applying a majority rule on the reference genomes taxa
    if [[ "${#TO_KNOWN_TAXA[@]}" -gt "0" ]]; then
        # Iterate over the unknown clusters in the associative array
        for UNKNOWN_TAXONOMY in "${!TO_KNOWN_TAXA[@]}"; do
            # Rebuild the folder path to the unknown taxonomy in the database
            UNKNOWN_TAXONOMY_PATH=$DBDIR/$(echo "${UNKNOWN_TAXONOMY}" | sed 's/|/\//g')
            # Take track of the reference genomes taxonomic labels
            KNOWN_TAXA=()
            # Split the list of genomes assigned to the unknown cluster
            GENOMES=($(echo ${TO_KNOWN_TAXA["${UNKNOWN_TAXONOMY}"]} | tr "," " "))
            for GENOME in ${GENOMES[@]}; do
                # Retrieve the genome taxonomy from the --taxa input file
                KNOWN_TAXA+=("$(grep -w "$GENOME" $TAXA | cut -d$'\t' -f2)")
            done
            # Get the most occurring taxonomic label
            ASSIGNED_TAXONOMY="$(printf '%s\n' "${TO_KNOWN_TAXA[@]}" | sort | uniq -c | sort -k1,1nr -k2 | awk '{print $2; exit}')"
            # Retrieve genomes file paths from the input list of genomes and update the genomes.fof file
            while read GENOMEPATH; do
                FULLPATH=$(realpath -s $GENOMEPATH)
                GENOMENAME="$(basename $FULLPATH)"
                GENOMENAME=${GENOMENAME%".$EXTENSION"}
                # If the current genome is in the list of genomes assigned to the unknown cluster
                if [[ " ${GENOMES[@]} " =~ " $GENOMENAME " ]]; then
                    cp $GENOMEPATH ${UNKNOWN_TAXONOMY_PATH}/genomes
                    GENOMEDBPATH=${UNKNOWN_TAXONOMY_PATH}/genomes/${GENOMENAME}.${EXTENSION}
                    if [[ ! "$EXTENSION" ~= "gz" ]]; then
                        gzip $GENOMEDBPATH
                        GENOMEDBPATH=${GENOMEDBPATH}.gz
                    fi
                    echo "$GENOMENAME : $GENOMEDBPATH" >> ${UNKNOWN_TAXONOMY_PATH}/genomes.fof
                    echo "$GENOMEDBPATH" >> ${UNKNOWN_TAXONOMY_PATH}/genomes.txt
                fi
            done <$INLIST

            # Split taxonomies to levels
            UNKNOWN_LEVELS=($(echo ${UNKNOWN_TAXONOMY} | tr "|" " "))
            ASSIGNED_LEVELS=($(echo ${ASSIGNED_TAXONOMY} | tr "|" " "))

            # Characterise unknown clusters
            for POS in $(seq 7 1); do
                if [[ "${UNKNOWN_LEVELS[$POS]}" = "${ASSIGNED_LEVELS[$POS]}" ]]; then
                    break
                else
                    # Get the partial levels for the new taxonomy
                    NEW_LEVELS=("${ASSIGNED_LEVELS[@]:$POS}")
                    # Build the partial levels path for the new taxonomy
                    NEW_LEVELS_SUBPATH="$(printf "/%s" "${NEW_LEVELS[@]}")"
                    # Trim the first slash out of the partial path
                    NEW_LEVELS_SUBPATH=${NEW_LEVELS_SUBPATH:1}
                    mkdir -p $DBDIR/${NEW_LEVELS_SUBPATH}

                    # Get the partial levels for the old unknown taxonomy
                    OLD_LEVELS=("${UNKNOWN_LEVELS[@]:$POS}")
                    # Build the partial levels path for the old unknown taxonomy
                    OLD_LEVELS_SUBPATH="$(printf "/%s" "${OLD_LEVELS[@]}")"
                    # Trim the first slash out of the partial path
                    OLD_LEVELS_SUBPATH=${OLD_LEVELS_SUBPATH:1}

                    # Move everything starting with folders to the new taxonomy folder
                    find $DBDIR/${OLD_LEVELS_SUBPATH} -type d -exec mv {} $DBDIR/${NEW_LEVELS_SUBPATH} \;
                    # Also move remaining files in case of species levels (pos 7)
                    # genomes.fof, genomes.txt, and mags.txt files
                    find $DBDIR/${OLD_LEVELS_SUBPATH} -type f -exec mv {} $DBDIR/${NEW_LEVELS_SUBPATH} \;

                    # Remove the old bloom filters root nodes
                    if [[ -f $DBDIR/${NEW_LEVELS_SUBPATH}/${UNKNOWN_LEVELS[$POS]}.bf ]]; then 
                        rm -f $DBDIR/${NEW_LEVELS_SUBPATH}/${UNKNOWN_LEVELS[$POS]}.bf
                    fi

                    # Finally remove the old taxonomy folder which is now empty
                    rm -rf $DBDIR/${OLD_LEVELS_SUBPATH}
                fi
            done
            
            # Add the new renamed cluster to the REBUILD array
            REBUILD+=(${ASSIGNED_TAXONOMY})
        done
    fi

    # Cluster unassigned genomes before rebuilding the updated lineages
    if [[ "${#UNASSIGNED[@]}" -gt "0" ]]; then
        # The kmers matrix already exists in case the dereplication step has been enabled
        # Otherwise, run kmtricks
        if [[ ! -f $TMPDIR/kmers_matrix.txt ]]; then
            GENOMES_FOF=$TMPDIR/genomes.fof
            # Unassigned genomes are MAGs only
            for MAGPATH in ${UNASSIGNED[@]}; do
                # Define a genomes.fof file with the unassigned genomes
                GENOME_NAME="$(basename $MAGPATH)"
                echo "${GENOME_NAME} : $MAGPATH" >> ${GENOMES_FOF}
            done

            # Run kmtricks to build the kmers matrix
            kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                    $TMPDIR/matrix \
                                    $NPROC \
                                    $TMPDIR/kmers_matrix.txt \
                                    $TMPDIR
            
            # Add header to the kmers matrix
            HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
            if [[ "${HEADER: -1}" = " " ]]; then
                # Trim the last character out of the header in xase of space
                HEADER="${HEADER:0:-1}"
            fi
            echo "#kmer $HEADER" > $TMPDIR/kmers_matrix_wh.txt
            cat $TMPDIR/kmers_matrix.txt >> $TMPDIR/kmers_matrix_wh.txt
            mv $TMPDIR/kmers_matrix_wh.txt $TMPDIR/kmers_matrix.txt
        fi

        # Cluster genomes according to the boundaries defined by the boundaries module
        # Define a cluster for each taxomomic level
        # Look at the genomes profiles and update or build new clusters 
        # Add full taxonomy to the REBUILD list
        # Update the genomes.fof and the mags.txt tables
        # Also report the CheckM statistics of the genomes in the new clusters
        # TODO
    fi

    # Check whether there is at least one lineage that must be rebuilt
    if [[ "${#REBUILD[@]}" -gt "0" ]]; then
        # Extract taxonomic levels from the list of taxa that must be rebuilt
        # Process all the species first, then all the genera, and so on up to the kingdom
        println "Updating the database\n"
        for POS in $(seq 7 1); do 
            TAXALIST=()
            for LABEL in ${REBUILD[@]}; do
                LEVELS=($(echo $LABEL | tr "|" " "))
                # Temporarily skip partial taxonomic labels if the number of levels is lower than POS
                # Waiting for the right POS
                if [[ "${#LEVELS[@]}" -lt "$POS" ]]; then
                    continue
                fi
                # Build the partial taxonomic label
                SUBLEVELS=("${LEVELS[@]:0:$POS}")
                TAXONOMY="$(printf "|%s" "${SUBLEVELS[@]}")"
                TAXALIST+=(${TAXONOMY:1})
            done
            # Remove duplicate entries
            TAXALIST=$(printf "%s\0" "${TAXALIST[@]}" | sort -uz | xargs -0 printf "%s ")
            # Rebuild the sequence bloom trees
            for TAXONOMY in ${TAXALIST[@]}; do
                println "\t%s\n" "$TAXONOMY"
                # Remove the index
                if [[ -d $TAXONOMY/index ]]; then 
                    rm -rf rm -rf $TAXONOMY/index
                fi
                # Remove the copy of the root node
                LEVELNAME="${TAXONOMY##*\|}"
                if [[ -f $TAXONOMY/${LEVELNAME}.bf ]]; then 
                    rm -f $TAXONOMY/${LEVELNAME}.bf
                fi
                # Rebuild the index
                if [[ "$POS" -eq "7" ]]; then
                    # In case of species
                    # Rebuild the index with kmtricks
                    kmtricks_index_wrapper $TAXONOMY/genomes.fof \
                                           $DBDIR \
                                           ${KMER_LEN} \
                                           ${FILTER_SIZE} \
                                           $NPROC
                else
                    # In case of all the other taxonomic levels
                    # Rebuild the index with howdesbt
                    howdesbt_wrapper $TAXONOMY ${FILTER_SIZE}
                fi
            done
        done
    else
        # In case nothing must be rebuilt
        println "No lineages have been updated!\n"
    fi
else
    # No input genomes survived from the quality control step and dereplication
    println "No input genomes available!\n"
    PIPELINE_END_TIME="$(date +%s.%3N)"
    PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
    println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"
    
    # Print credits before stopping the pipeline execution
    credits

    exit 0
fi

# Cleanup temporary data
if $CLEANUP; then
    # Remove tmp folder
    if [[ -d $TMPDIR ]]; then
        println "Cleaning up temporary folder:\n"
        println "\t%s\n" "$TMPDIR"
        rm -rf $TMPDIR
    fi
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0