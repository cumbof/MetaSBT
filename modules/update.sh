#!/bin/bash
#title          :update
#description    :Update the index with new genomes
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="May 18, 2022"
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
CHECKM_BINEXT="fna"

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --checkm-completeness=*)
            # CheckM completeness
            CHECKM_COMPLETENESS="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_COMPLETENESS}" =~ "?" ]]; then
                printf "update helper: --checkm-completeness=num\n\n"
                printf "\tInput genomes must have a minimum completeness percentage before being processed and added to the database.\n"
                printf "\tPlease note that this argument is used only in case of MAGs as input genomes.\n\n"
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
                printf "\tInput genomes must have a maximum contamination percentage before being processed and added to the database.\n"
                printf "\tPlease note that this argument is used only in case of MAGs as input genomes.\n\n"
                exit 0
            fi
            # Check whether --checkm-contamination is an integer
            if [[ ! ${CHECKM_CONTAMINATION} =~ ^[0-9]+$ ]] || [[ "${CHECKM_CONTAMINATION}" -gt "100" ]]; then
                printf "Argument --checkm-completeness must be a positive integer lower than 100 (up to 0)\n"
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
                printf "Please have a look at the helper of the --type argument for a list of valid genome types\n"
                exit 1
            fi
            ;;
        -v|--version)
            # Print pipeline version
            printf "update version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        --work-dir=*)
            # Working directory
            WORKDIR="${ARG#*=}"
            # Define helper
            if [[ "${WORKDIR}" =~ "?" ]]; then
                printf "update helper: --work-dir=directory\n\n"
                printf "\tThis is the working directory that will contain genomes organised by species and their index produced by kmtricks.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            WORKDIR="$( cd "$( dirname "${WORKDIR}" )" &> /dev/null && pwd )"/"$( basename $WORKDIR )"
            # Trim the last slash out of the path
            WORKDIR="${WORKDIR%/}"
            ;;
        --xargs-nproc=*)
            # Max number of independent runs of kmtricks
            XARGS_NPROC="${ARG#*=}"
            # Define helper
            if [[ "${XARGS_NPROC}" =~ "?" ]]; then
                printf "update helper: --xargs-nproc=num\n\n"
                printf "\tThis refers to the number of xargs processes used for launching independent runs of kmtricks.\n"
                printf "\tDefault: --xargs-nproc=1\n\n"
                exit 0
            fi
            # Check whether --xargs-nproc is an integer
            if [[ ! $XARGS_NPROC =~ ^[0-9]+$ ]] || [[ "$XARGS_NPROC" -eq "0" ]]; then
                printf "Argument --xargs-nproc must be a positive integer greater than 0\n"
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

check_dependencies
if [[ "$?" -gt "0" ]]; then
    exit 1
fi

# Locate database directory
DBDIR=$WORKDIR/db

# Apply CheckM in case of MAGs only
if [[ "$TYPE" = "MAGs" ]]; then
    # Input genomes must be quality-controlled before being added to the database
    if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
        # Run CheckM
        printf 'Running CheckM\n'
        CHECKM_START_TIME="$(date +%s.%3N)"
        mkdir -p $WORKDIR/checkm/tmp
        # Split the set of bins in chunks with 1000 genomes at most
        println '\tOrganising genomes in chunks\n'
        split --numeric-suffixes=1 --lines=1000 --suffix-length=3 --additional-suffix=.txt $INLIST $WORKDIR/checkm/tmp/bins_
        CHUNKS=`ls "$WORKDIR"/checkm/tmp/bins_*.txt 2>/dev/null | wc -l`
        CHECKMTABLES=""
        for filepath in $WORKDIR/checkm/tmp/bins_*.txt; do
            # Retrieve chunk id
            filename="$(basename $filepath)"
            suffix="${filename#*_}"
            suffix="${suffix%.txt}"
            println '\tProcessing chunk %s/%s\n' "$((10#$suffix))" "$CHUNKS"
            # Create chunk folder
            mkdir -p $WORKDIR/checkm/tmp/bins_${suffix}
            for bin in `sed '/^$/d' $filepath`; do
                # Make a symbolic link to the 1000 genomes for the current chunk
                ln -s $bin $WORKDIR/checkm/tmp/bins_${suffix}
            done
            # Create the CheckM folder for the current chunk
            mkdir -p $WORKDIR/checkm/run_$suffix/
            # Run CheckM
            checkm lineage_wf -t ${NPROC} \
                            -x ${CHECKM_BINEXT} \
                            --pplacer_threads ${NPROC} \
                            --tab_table -f $WORKDIR/checkm/run_${suffix}.tsv \
                            $WORKDIR/checkm/tmp/bins_${suffix} $WORKDIR/checkm/run_$suffix
            # Take trace of the CheckM output tables
            CHECKMTABLES=$WORKDIR/checkm/run_${suffix}.tsv,$CHECKMTABLES
        done
        CHECKM_END_TIME="$(date +%s.%3N)"
        CHECKM_ELAPSED="$(bc <<< "${CHECKM_END_TIME}-${CHECKM_START_TIME}")"
        println '\tTotal elapsed time: %s\n\n' "$(displaytime ${CHECKM_ELAPSED})"
        # Trim the last comma out of the list of CheckM output table file paths
        CHECKMTABLES="${CHECKMTABLES%?}"

        # Read all the CheckM output tables and filter genomes on their completeness and contamination scores 
        # according to the those provided in input

        # Create a new INLIST file with the new list of genomes

    fi
fi

# Use kmtricks to build a kmer matrix and compare input genomes
# Discard genomes with the same kmers (dereplication)
# Create a new INLIST file with the new list of genomes

# In case of MAGs:
#   For each of the input genomes
#   Make a query to establish the closest group for each taxonomic level in the tree
#   If it is close enough to a group in a taxonomic level, remember that this genome must be assigned to the same group
#   Otherwise, mark the genome as unassigned
#
#   Assigned all the genomes to the closest groups and rebuild the SBTs
#   Define new groups by looking at the kmer matrix of kmtricks for the unassigned genomes

# In case of references:
#   If the taxonomic label of the input genome already exists in the database
#       Add the new genome to the assigned taxonomic branch and rebuild the SBTs
#   Otherwise,
#       Compare the genome with genomes in newly defined clusters at the level of the known lineage
#       In case there is nothing close to the genome, create a new branch with the new taxonomy for that genome
#       Otherwise, merge the unknown clusters with the new genomes and assign the new taxonomy

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

exit 0