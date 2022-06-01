#!/bin/bash
#title          :update
#description    :Update the index with new genomes
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="Jun 1, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

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
            TAX_BOUNDARIES="$(grep "${TAXLABEL}"$'\t' $BOUNDARIES)"
        else
            # Add a pipe at the end of the taxonomic label for all the other levels
            TAX_BOUNDARIES="$(grep "${TAXLABEL}|" $BOUNDARIES)"
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
            TAXLABEL=$(printf "|%s" "${TAXONOMY_SPLIT[@]}")
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
# Input genome extension is forced to be fna
EXTENSION="fna"
# Dereplicate input genomes
# Use kmtricks to build the kmer matrix on the input set of genomes
# Remove genomes with identical set of kmers
DEREPLICATE=false
# Remove temporary data at the end of the pipeline
CLEANUP=false

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --boundaries=*)
            # Clusters boundaries computed through the boundaries module
            BOUNDARIES="${ARG#*=}"
            # Define helper
            if [[ "${BOUNDARIES}" =~ "?" ]]; then
                printf "update helper: --boundaries=file\n\n"
                printf "\tPath to the output table produced by the boundaries module.\n"
                printf "\tIt is required in case of MAGs as input genomes only\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $BOUNDARIES ]]; then
                printf "Input file does not exist!\n"
                printf "--boundaries=%s\n" "$BOUNDARIES"
                exit 1
            fi
            ;;
        --boundary-uncertainty=*)
            # Boundary uncertainty percentage
            BOUNDARY_UNCERTAINTY_PERC="${ARG#*=}"
            # Define helper
            if [[ "${BOUNDARY_UNCERTAINTY_PERC}" =~ "?" ]]; then
                printf "update helper: --boundary-uncertainty=num\n\n"
                printf "\tDefine the percentage of kmers to enlarge and reduce boundaries\n\n"
                exit 0
            fi
            # Check whether --boundary-uncertainty is an integer between 0 and 100
            if [[ ! ${BOUNDARY_UNCERTAINTY_PERC} =~ ^[0-9]+$ ]] || [[ "${BOUNDARY_UNCERTAINTY_PERC}" -lt "0" ]] || [[ "${BOUNDARY_UNCERTAINTY_PERC}" -gt "100" ]]; then
                printf "Argument --boundary-uncertainty must be a positive integer between 0 and 100\n"
                exit 1
            fi
            ;;
        --checkm-completeness=*)
            # CheckM completeness
            CHECKM_COMPLETENESS="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_COMPLETENESS}" =~ "?" ]]; then
                printf "update helper: --checkm-completeness=num\n\n"
                printf "\tInput genomes must have a minimum completeness percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-completeness is an integer
            if [[ ! ${CHECKM_COMPLETENESS} =~ ^[0-9]+$ ]] || [[ "${CHECKM_COMPLETENESS}" -gt "100" ]]; then
                printf "Argument --checkm-completeness must be a positive integer greater than 0 (up to 100)\n"
                exit 1
            fi
            ;;
        --checkm-contamination=*)
            # CheckM contamination
            CHECKM_CONTAMINATION="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_CONTAMINATION}" =~ "?" ]]; then
                printf "update helper: --checkm-contamination=num\n\n"
                printf "\tInput genomes must have a maximum contamination percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-contamination is an integer
            if [[ ! ${CHECKM_CONTAMINATION} =~ ^[0-9]+$ ]] || [[ "${CHECKM_CONTAMINATION}" -gt "100" ]]; then
                printf "Argument --checkm-completeness must be a positive integer lower than 100 (up to 0)\n"
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
            if [[ "${DBDIR}" =~ "?" ]]; then
                printf "update helper: --db-dir=directory\n\n"
                printf "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "${DBDIR}" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
            # Trim the last slash out of the path
            DBDIR="${DBDIR%/}"
            # Check whether the input directory exist
            if [[ ! -d $DBDIR ]]; then
                printf "Input folder does not exist!\n"
                printf "--db-dir=%s\n" "$DBDIR"
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
            if [[ "${EXTENSION}" =~ "?" ]]; then
                printf "update helper: --extension=value\n\n"
                printf "\tSpecify the input genome files extension.\n"
                printf "\tAll the input genomes must have the same file extension before running this module.\n\n"
                exit 0
            fi
            # Allowed extensions: "fa", "fasta", "fna" and gzip compressed formats
            EXTENSION_LIST=("fa" "fa.gz" "fasta" "fasta.gz" "fna" "fna.gz")
            if ! echo ${EXTENSION_LIST[@]} | grep -w -q $EXTENSION; then
                printf "File extension \"%s\" is not allowed!\n" "$TYPE"
                printf "Please have a look at the helper of the --extension argument for a list of valid input file extensions.\n"
                exit 1
            fi
            ;;
        --filter-size=*)
            # Bloom filter size
            FILTER_SIZE="${ARG#*=}"
            # Define helper
            if [[ "${FILTER_SIZE}" =~ "?" ]]; then
                printf "update helper: --filter-size=num\n\n"
                printf "\tThis is the size of the bloom filters.\n"
                printf "\tIt must be the same size used while running the index module.\n\n"
                exit 0
            fi
            # Check whether --filter-size is an integer
            if [[ ! ${FILTER_SIZE} =~ ^[0-9]+$ ]] || [[ "${FILTER_SIZE}" -eq "0" ]]; then
                printf "Argument --filter-size must be a positive integer greater than 0\n"
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
            if [[ "${INLIST}" =~ "?" ]]; then
                printf "update helper: --input-list=file\n\n"
                printf "\tThis file contains the list of paths to the new genomes that will be added to the database.\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $BOUNDARIES ]]; then
                printf "Input file does not exist!\n"
                printf "--input-list=%s\n" "$INLIST"
                exit 1
            fi
            ;;
        --kingdom=*)
            # Consider genomes whose lineage belong to a specific kingdom only
            KINGDOM="${ARG#*=}"
            # Define helper
            if [[ "${KINGDOM}" =~ "?" ]]; then
                printf "update helper: --kingdom=value\n\n"
                printf "\tSelect the kingdom that will be updated.\n\n"
                exit 0
            fi
            ;;
        --kmer-len=*)
            # Length of the kmers
            KMER_LEN="${ARG#*=}"
            # Define helper
            if [[ "${KMER_LEN}" =~ "?" ]]; then
                printf "update helper: --kmer-len=num\n\n"
                printf "\tThis is the length of the kmers used for building bloom filters.\n\n"
                exit 0
            fi
            # Check whether --kmer-len is an integer
            if [[ ! ${KMER_LEN} =~ ^[0-9]+$ ]] || [[ "${KMER_LEN}" -eq "0" ]]; then
                printf "Argument --kmer-len must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --nproc=*)
            # Max nproc for all parallel instructions
            NPROC="${ARG#*=}"
            # Define helper
            if [[ "${NPROC}" =~ "?" ]]; then
                printf "update helper: --nproc=num\n\n"
                printf "\tThis argument refers to the number of processors used for parallelizing the pipeline when possible.\n"
                printf "\tDefault: --nproc=1\n\n"
                exit 0
            fi
            # Check whether --proc is an integer
            if [[ ! $NPROC =~ ^[0-9]+$ ]] || [[ "$NPROC" -eq "0" ]]; then
                printf "Argument --nproc must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --taxa=*)
            # Input file with the mapping between input reference genome IDs and their taxonomic label
            TAXA="${ARG#*=}"
            # Define helper
            if [[ "${TAXA}" =~ "?" ]]; then
                printf "update helper: --taxa=file\n\n"
                printf "\tInput file with the mapping between input genome IDs and their taxonomic label.\n"
                printf "\tThis is used in case of reference genomes only \"--type=references\".\n\n"
                exit 0
            fi
            # Check whether the input file exists
            if [[ ! -f $BOUNDARIES ]]; then
                printf "Input file does not exist!\n"
                printf "--taxa=%s\n" "$TAXA"
                exit 1
            fi
            ;;
        --tmp-dir=*)
            # Temporary folder
            TMPDIR="${ARG#*=}"
            # Define helper
            if [[ "${TMPDIR}" =~ "?" ]]; then
                printf "update helper: --tmp-dir=directory\n\n"
                printf "\tPath to the folder for storing temporary data.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            TMPDIR="$( cd "$( dirname "${TMPDIR}" )" &> /dev/null && pwd )"/"$( basename $TMPDIR )"
            # Trim the last slash out of the path
            TMPDIR="${TMPDIR%/}"
            ;;
        --type=*)
            # Input genomes type: references or MAGs
            TYPE="${ARG#*=}"
            # Define helper
            if [[ "${TYPE}" =~ "?" ]]; then
                printf "update helper: --type=value\n\n"
                printf "\tDefine the nature of the input genomes. Only \"references\" and \"MAGs\" are allowed.\n\n"
                exit 0
            fi
            # Allowed genome types: "references" and "MAGs"
            if [[ ! "$TYPE" = "MAGs" ]] && [[ ! "$TYPE" = "references" ]]; then
                printf "Input type \"%s\" is not allowed!\n" "$TYPE"
                printf "Please have a look at the helper of the --type argument for a list of valid genome types.\n"
                exit 1
            fi
            ;;
        *)
            printf "update: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

printf "update version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false
if [[ "$?" -gt "0" ]]; then
    printf "Unsatisfied software dependencies!\n\n"
    printf "Please run the following command for a list of required external software dependencies:\n\n"
    printf "\t$ meta-index --resolve-dependencies\n\n"
    
    exit 1
fi

# Create temporary folder
mkdir -p $TMPDIR

# Count how many genomes in input
HOW_MANY=$(cat ${INLIST} | wc -l)
printf "Processing %s genomes\n" "${HOW_MANY}"
printf "\t%s\n\n" "${INLIST}"

# Input genomes must be quality-controlled before being added to the database
if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
    CHECKMTABLES="" # This is the result of the "run_checkm" function as a list of file paths separated by comma
    run_checkm $INLIST $EXTENSION $TMPDIR $NPROC
    
    # Read all the CheckM output tables and filter genomes on their completeness and contamination scores 
    # according to the those provided in input
    touch ${TMPDIR}/genomes_qc.txt
    CHECKMTABLES_ARR=($(echo $CHECKMTABLES | tr "," " "))
    for table in ${CHECKMTABLES_ARR[@]}; do
        # Skip header line with sed
        # Discard genomes according to their completeness and contamination
        sed 1d $table | while read line; do
            # Retrieve completeness and contamination of current genome
            GENOMEID="$(echo $line | cut -d$'\t' -f1)"          # Get value under column 1
            COMPLETENESS="$(echo $line | cut -d$'\t' -f12)"     # Get value under column 12
            CONTAMINATION="$(echo $line | cut -d$'\t' -f13)"    # Get value under column 13
            if [[ "$COMPLETENESS" -ge "${CHECKM_COMPLETENESS}" ]] && [[ "$CONTAMINATION" -le "${CHECKM_CONTAMINATION}" ]]; then
                # Current genome passed the quality control
                grep -w "${GENOMEID}.${EXTENSION}" ${INLIST} >> ${TMPDIR}/genomes_qc.txt
            fi
        done
    done

    # Create a new INLIST file with the new list of genomes
    INLIST=${TMPDIR}/genomes_qc.txt
    # Count how many genomes passed the quality control
    HOW_MANY=$(cat ${INLIST} | wc -l)

    printf "\t%s genomes passed the quality control\n" "${HOW_MANY}"
    printf "\t\tMinimum completeness: %s\n" "${CHECKM_COMPLETENESS}"
    printf "\t\tMaximum contamination: %s\n" "${CHECKM_CONTAMINATION}"
fi

# Use kmtricks to build a kmer matrix and compare input genomes
# Discard genomes with the same set of kmers (dereplication)
if ${DEREPLICATE}; then
    # Count how many input genomes survived in case they have been quality controlled
    HOW_MANY=$(cat ${INLIST} | wc -l)
    if [[ "${HOW_MANY}" -gt "1" ]]; then        
        GENOMES_FOF=${TMPDIR}/genomes.fof
        # Build a fof file with the list of input genomes
        while read -r genomepath; do
            GENOME_NAME="$(basename $genomepath)"
            echo "${GENOME_NAME} : $genomepath" >> ${GENOMES_FOF}
        done <$INLIST

        if [[ -f ${GENOMES_FOF} ]]; then
            printf "\nDereplicating %s input genomes\n" "${HOW_MANY}"
            # Run kmtricks to build the kmers matrix
            kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                    ${TMPDIR}/matrix \
                                    ${NPROC} \
                                    ${TMPDIR}/kmers_matrix.txt
            
            # Add header to the kmers matrix
            HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
            echo "#kmer $HEADER" > ${TMPDIR}/kmers_matrix_wh.txt
            cat ${TMPDIR}/kmers_matrix.txt >> ${TMPDIR}/kmers_matrix_wh.txt
            mv ${TMPDIR}/kmers_matrix_wh.txt ${TMPDIR}/kmers_matrix.txt
            # Remove duplicate columns
            cat ${TMPDIR}/kmers_matrix.txt | transpose | awk '!seen[substr($0, index($0, " "))]++' | transpose > ${TMPDIR}/kmers_matrix_dereplicated.txt
            
            # Dereplicate genomes
            while read -r genome; do
                if head -n1 ${TMPDIR}/kmers_matrix_dereplicated.txt | grep -q -w "$genome"; then
                    # Rebuild the absolute path to the genome file
                    realpath -s $(grep "^${genome} " ${TMPDIR}/genomes.fof | cut -d' ' -f3) >> ${TMPDIR}/genomes_dereplicated.txt
                fi
            done <<<"$(head -n1 ${TMPDIR}/kmers_matrix.txt | tr " " "\n")"
            
            # Use the new list of dereplicated genomes
            if [[ -f ${TMPDIR}/genomes_dereplicated.txt ]]; then
                INLIST=${TMPDIR}/genomes_dereplicated.txt
            fi

            # Count how many genomes passed the dereplication
            HOW_MANY=$(cat ${INLIST} | wc -l)
            printf "\t%s genomes passed the dereplication\n" "${HOW_MANY}"
        fi
    else
        # No input genomes survived from the quality control step
        printf "No input genomes available!\n"
        PIPELINE_END_TIME="$(date +%s.%3N)"
        PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
        printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"
        
        # Print credits before stopping the pipeline execution
        credits
        
        exit 0
    fi
fi

# Count how many input genomes survived in case they have been quality controlled and dereplicated
HOW_MANY=$(cat ${INLIST} | wc -l)
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
    while read GENOMEPATH; then
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

        printf "\nProfiling %s\n" "$GENOMENAME"
        # Run the profiler to establish the closest genome and the closest group for each taxonomic level in the tree
        . ${SCRIPT_DIR}/profile.sh --input-file=$FILEPATH \
                                   --input-id=$FILEPATH \
                                   --tree=${DBDIR}/k__${KINGDOM}/index/index.detbrief.sbt
                                   --expand \
                                   --output-dir=$TMPDIR/profiling/ \
                                   --output-prefix=$GENOMENAME \
                                   > /dev/null 2>&1 # Silence the profiler
        # Define output file path with profiles
        PROFILE=${TMPDIR}/profiling/${GENOMENAME}__profiles.tsv
        printf "\t%s\n" "${PROFILE}"
        
        if [[ -f $PROFILE ]]; then
            # Discard the genome if the score is 1.0 compared to the closest match at the species level
            # Do not discard the input genome if it is a reference genome and the closest genome in the database
            # is a MAG with a perfect overlap of kmers (score 1.0)
            CLOSEST_GENOME_DATA="$(grep "${GENOMENAME}.${GENOMEEXT}" $PROFILE | grep -w "genome")"
            CLOSEST_GENOME="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f3)"
            CLOSEST_GENOME_COMMON_KMERS="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f4)"
            CLOSEST_GENOME_SCORE="$(echo "${CLOSEST_GENOME_DATA}" | cut -d$'\t' -f5)"
            # Print closest genome
            printf "\tClosest genome: %s (common kmers: %s; score: %s)\n" "${CLOSEST_GENOME}" "${CLOSEST_GENOME_COMMON_KMERS}" "${CLOSEST_GENOME_SCORE}"
            # Reconstruct the closest taxa
            LEVELS=("kingdom" "phylum" "class" "order" "family" "genus" "species")
            CLOSEST_TAXA=""
            CLOSEST_COMMON_KMERS=""
            printf "\tClosest lineage:\n"
            for level in ${LEVELS[@]}; do
                CLOSEST_LEVEL_DATA="$(grep "${GENOMENAME}.${GENOMEEXT}" $PROFILE | grep -w "$level")"
                CLOSEST_LEVEL="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f3)"
                CLOSEST_LEVEL_COMMON_KMERS="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f4)"
                CLOSEST_LEVEL_SCORE="$(echo "${CLOSEST_LEVEL_DATA}" | cut -d$'\t' -f5)"
                CLOSEST_TAXA=${CLOSEST_TAXA}"|"${CLOSEST_LEVEL}
                # Print profiler result
                printf "\t\t%s: %s (common kmers: %s; score: %s)\n" "${level}" "${CLOSEST_LEVEL}" "${CLOSEST_LEVEL_COMMON_KMERS}" "${CLOSEST_LEVEL_SCORE}"
                # Keep track of the species score
                if [[ "${CLOSEST_LEVEL}" = s__* ]]; then
                    CLOSEST_COMMON_KMERS=${CLOSEST_LEVEL_COMMON_KMERS}
                fi
            done
            # Trim the first pipe out of the closest taxonomic label
            CLOSEST_TAXA=${CLOSEST_TAXA:1}
            # Define taxa path
            CLOSEST_TAXADIR=${DBDIR}/$(echo "${CLOSEST_TAXA}" | sed 's/|/\//g')

            SKIP_GENOME=false
            # Check whether the input genome must be discarded
            if [[ "${CLOSEST_GENOME_SCORE}" -eq "1.0" ]]; then
                if [[ "$TYPE" = "MAGs" ]]; then
                    if [[ -f "${CLOSEST_TAXADIR}/genomes.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/genomes.txt; then
                        # If the input genome is a MAG and the closest genome is a reference
                        # Discard the input genome
                        ${SKIP_GENOME}=true
                        # Print the reason why the input genome has been discarded
                        printf "\tDiscarding genome:\n"
                        printf "\t\tInput genome is a MAG and the closest genome is a reference genome\n"
                    elif [[ -f "${CLOSEST_TAXADIR}/mags.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/mags.txt; then
                        # If the input genome is a MAG and the closest genome is a MAG
                        # Discard the input genome
                        ${SKIP_GENOME}=true
                        # Print the reason why the input genome has been discarded
                        printf "\tDiscarding genome:\n"
                        printf "\t\tInput genome and the closest genome are both MAGs\n"
                    fi
                elif [[ "$TYPE" = "references" ]]; then
                    if [[ -f "${CLOSEST_TAXADIR}/genomes.txt" ]] && grep -q "" ${CLOSEST_TAXADIR}/genomes.txt; then
                        # If the input genome is a reference and the closest genome is a reference
                        # Discard the input genome
                        ${SKIP_GENOME}=true
                        # Print the reason why the input genome has been discarded
                        printf "\tDiscarding genome:\n"
                        printf "\t\tInput genome and closest genome are both reference genomes\n"
                    fi
                fi
            fi

            if ! ${SKIP_GENOME}; then
                # Retrieve the min and max common kmers for the closest taxa
                MIN_BOUND=""; MAX_BOUND=""
                get_boundaries $BOUNDARIES $CLOSEST_TAXA

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
                        # Do not remove the index here because there could be other input genomes with the same taxonomic label
                        REBUILD+=(${CLOSEST_TAXA})
                        # Print assignment message
                        printf "\tAssignment:\n"
                        printf "\t\t%s\n\n" "${CLOSEST_TAXA}"
                    else
                        # Mark the genome as unassigned
                        UNASSIGNED+=($GENOMEPATH)
                        # Print unassignment message
                        printf "\tUnassigned\n\n"
                    fi
                elif [[ "$TYPE" = "references" ]]; then
                    # Retrieve the taxonomic label from --taxa input file
                    TAXALABEL="$(grep -w "$GENOMENAME" $TAXA | cut -d$'\t' -f2)"
                    # Define taxonomy folder
                    TAXDIR=${DBDIR}/$(echo "${TAXALABEL}" | sed 's/|/\//g')
                    TARGETGENOME=${TAXDIR}/genomes/${GENOMENAME}.${GENOMEEXT}
                    # Print input genome taxonomic label
                    printf "\tTaxonomic label:\n"
                    printf "\t\t%s\n" "$TAXALABEL"

                    if [[ "${CLOSEST_COMMON_KMERS}" -le "${MAX_BOUND}" ]] && [[ "${CLOSEST_COMMON_KMERS}" -ge "${MIN_BOUND}" ]]; then
                        # Check whether the closest taxonomy contains any reference genomes
                        HOW_MANY_REFERENCES=$(cat ${CLOSEST_TAXADIR}/genomes.txt | wc -l)
                        if [[ "${HOW_MANY_REFERENCES}" -eq "0" ]]; then
                            # If the closest genome belongs to a new cluster with no reference genomes
                            # Assign the current reference genome to the new cluster and rename its lineage with the taxonomic label of the reference genome
                            # Do not change lineage here because there could be more than one reference genome assigned to the same unknown cluster
                            # Assign a new taxonomy according to a majority rule applied on the reference genomes taxa
                            if [[ ! -z "${TO_KNOWN_TAXA[${CLOSEST_TAXA}]}" ]]; then
                                TO_KNOWN_TAXA[${CLOSEST_TAXA}]="$GENOMENAME"
                            else
                                # If already exists in the associative array
                                # Concatenate values with a comma
                                TO_KNOWN_TAXA[${CLOSEST_TAXA}]=${TO_KNOWN_TAXA[${CLOSEST_TAXA}]}",$GENOMENAME"
                            fi
                            # Print the assignment
                            # Taxonomy will change according to the majority voting result
                            printf "\tAssignment:\n"
                            printf "\t\t%s\n\n" "${CLOSEST_TAXA}"
                        elif [[ "${HOW_MANY_REFERENCES}" -ge "1" ]]; then
                            # If the closest genome belong to a cluster with at least a reference genome
                            if [[ "${TAXALABEL}" != "${CLOSEST_TAXA}" ]]; then
                                # If the taxonomic labels of the current reference genome and that one of the closest genome do not match
                                # Report the inconsistency
                                printf "\tInconsistency found:\n"
                                printf "\t\tInput genome: %s\n" "${TAXALABEL}"
                                printf "\t\tClosest lineage: %s\n" "${CLOSEST_TAXA}"
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
                            printf "\tAssignment:\n"
                            printf "\t\t%s\n\n" "${CLOSEST_TAXA}"
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
            GENOMES=($(echo ${TO_KNOWN_TAXA[$UNKNOWN_TAXONOMY]} | tr "," " "))
            for GENOME in GENOMES; do
                # Retrieve the genome taxonomy from the --taxa input file
                KNOWN_TAXA+=("$(grep -w "$GENOME" $TAXA | cut -d$'\t' -f2)")
            done
            # Get the most occurring taxonomic label
            ASSIGNED_TAXONOMY="$(printf '%s\n' "${TO_KNOWN_TAXA[@]}" | sort | uniq -c | sort -k1,1nr -k2 | awk '{print $2; exit}')"
            # Retrieve genomes file paths from the input list of genomes and update the genomes.fof file
            while read GENOMEPATH; then
                FULLPATH=$(realpath -s $GENOMEPATH)
                GENOMENAME="$(basename $FULLPATH)"
                GENOMENAME=${GENOMENAME%"$EXTENSION"}
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
                    NEW_LEVELS_SUBPATH=$(printf "/%s" "${NEW_LEVELS[@]}")
                    # Trim the first slash out of the partial path
                    NEW_LEVELS_SUBPATH=${NEW_LEVELS_SUBPATH:1}
                    mkdir -p $DBDIR/${NEW_LEVELS_SUBPATH}

                    # Get the partial levels for the old unknown taxonomy
                    OLD_LEVELS=("${UNKNOWN_LEVELS[@]:$POS}")
                    # Build the partial levels path for the old unknown taxonomy
                    OLD_LEVELS_SUBPATH=$(printf "/%s" "${OLD_LEVELS[@]}")
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
        if [[ ! -f ${TMPDIR}/kmers_matrix.txt ]]; then
            GENOMES_FOF=${TMPDIR}/genomes.fof
            # Unassigned genomes are MAGs only
            for MAGPATH in ${UNASSIGNED[@]}; do
                # Define a genomes.fof file with the unassigned genomes
                GENOME_NAME="$(basename $MAGPATH)"
                echo "${GENOME_NAME} : $MAGPATH" >> ${GENOMES_FOF}
            done

            # Run kmtricks to build the kmers matrix
            kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                    ${TMPDIR}/matrix \
                                    ${NPROC} \
                                    ${TMPDIR}/kmers_matrix.txt
            
            # Add header to the kmers matrix
            HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
            echo "#kmer $HEADER" > ${TMPDIR}/kmers_matrix_wh.txt
            cat ${TMPDIR}/kmers_matrix.txt >> ${TMPDIR}/kmers_matrix_wh.txt
            mv ${TMPDIR}/kmers_matrix_wh.txt ${TMPDIR}/kmers_matrix.txt
        fi

        # Cluster genomes according to the boundaries defined by the boundaries module
        # Define a cluster for each taxomomic level
        # Look at the genomes profiles and update or build new clusters 
        # Add full taxonomy to the REBUILD list
        # Update the genomes.fof and the mags.txt tables
        # TODO
    fi

    # Check whether there is at least one lineage that must be rebuilt
    if [[ "${#REBUILD[@]}" -gt "0" ]]; then
        # Extract taxonomic levels from the list of taxa that must be rebuilt
        # Process all the species first, then all the genera, and so on up to the kingdom
        printf "Updating the database\n"
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
                TAXONOMY=$(printf "|%s" "${SUBLEVELS[@]}")
                TAXALIST+=(${TAXONOMY:1})
            done
            # Remove duplicate entries
            TAXALIST=$(printf "%s\0" "${TAXALIST[@]}" | sort -uz | xargs -0 printf "%s ")
            # Rebuild the sequence bloom trees
            for TAXONOMY in ${TAXALIST[@]}; do
                printf "\t%s\n" "$TAXONOMY"
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
                    kmtricks_index_wrapper ${TAXONOMY}/genomes.fof \
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
        printf "No lineages have been updated!\n"
    fi
else
    # No input genomes survived from the quality control step and dereplication
    printf "No input genomes available!\n"
    PIPELINE_END_TIME="$(date +%s.%3N)"
    PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
    printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"
    
    # Print credits before stopping the pipeline execution
    credits

    exit 0
fi

# Cleanup temporary data
if ${CLEANUP}; then
    # Remove tmp folder
    if [[ -d $TMPDIR ]]; then
        printf "Cleaning up temporary folder:\n"
        printf "\t%s\n" "${TMPDIR}"
        rm -rf ${TMPDIR}
    fi
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0