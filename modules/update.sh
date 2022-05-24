#!/bin/bash
#title          :update
#description    :Update the index with new genomes
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="May 24, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

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
            # Allowed extensions: "fa", "fasta", "fna" and GZip compressed formats
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
                printf "\tThis is the size of the bloom filters.\n\n"
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
        --license)
            # Print license
            printf "%s\n" "$(cat ${SCRIPT_DIR}/../LICENSE)"
            exit 0
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
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies
            exit $?
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

printf "update version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies
if [[ "$?" -gt "0" ]]; then
    exit 1
fi

# Create temporary folder
mkdir -p $TMPDIR

# Input genomes must be quality-controlled before being added to the database
if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
    CHECKMTABLES="" # This is the result of the "run_checkm" function as a list of file paths separated by comma
    run_checkm $INLIST $EXTENSION $TMPDIR $NPROC
    
    # Read all the CheckM output tables and filter genomes on their completeness and contamination scores 
    # according to the those provided in input
    touch ${TMPDIR}/genomes_qc.txt
    for table in ${CHECKMTABLES//,/ }; do
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
            printf "\tDereplicating %s input genomes" "${HOW_MANY}"
            # Build matrix
            kmtricks pipeline --file ${GENOMES_FOF} \
                              --run-dir ${TMPDIR}/matrix \
                              --mode kmer:count:bin \
                              --hard-min 1 \
                              --cpr \
                              --threads ${NPROC}
            # Aggregate
            kmtricks aggregate --run-dir ${TMPDIR}/matrix \
                               --matrix kmer \
                               --format text \
                               --cpr-in \
                               --sorted \
                               --threads ${NPROC} > ${TMPDIR}/kmers_matrix.txt
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
        fi
    else
        # No input genomes survived from the quality control step
        printf "\tNo input genomes available!"
        exit 0
    fi
fi

if [[ "$TYPE" = "MAGs" ]]; then
    # For each of the input genomes
    # Make a query to establish the closest group for each taxonomic level in the tree
    # If it is close enough to a group in a taxonomic level, remember that this genome must be assigned to the same group
    # Otherwise, mark the genome as unassigned
    #
    # Assigned all the genomes to the closest groups and rebuild the SBTs
    # Define new groups by looking at the kmer matrix of kmtricks for the unassigned genomes
elif [[ "$TYPE" = "references" ]]; then
    # If the taxonomic label of the input genome already exists in the database
    #   Add the new genome to the assigned taxonomic branch and rebuild the SBTs
    # Otherwise,
    #   Compare the genome with genomes in newly defined clusters at the level of the known lineage
    #   In case there is nothing close to the genome, create a new branch with the new taxonomy for that genome
    #   Otherwise, merge the unknown clusters with the new genomes and assign the new taxonomy
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