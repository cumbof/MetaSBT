#!/bin/bash
#title          :utils
#description    :Some utility functions
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="Jun 8, 2022"
VERSION="0.1.0"

# Check for external software dependencies
check_dependencies () {
    VERBOSE=$1                   # Print messages if true
    DEPENDENCY_FILEPATHS=("$@")  # Path to the dependency file
    
    # Remove the first argument from the dependency list
    DEPENDENCY_FILEPATHS=(${DEPENDENCY_FILEPATHS[@]:1})

    # Print messages if verbose
    if $VERBOSE; then
        println "Checking for software dependencies\n"
    fi
    
    # Check whether the list of dependency files contains at least one path
    if [[ "${#DEPENDENCY_FILEPATHS[@]}" -gt "0" ]]; then
        SOFTWARE_DEPENDENCIES=()
        PYTHON_FOUND=false
        # For each of the dependency files
        for DEPENDENCY_FILEPATH in "${DEPENDENCY_FILEPATHS[@]}"; do
            # Load the set of software dependencies
            while read line; do
                SOFTWARE_DEPENDENCIES+=("$line")
                # Also check for python modules in case Python is a requirement
                if [[ "$line" = "python3" ]]; then
                    PYTHON_FOUND=true
                    # In case of python requirement
                    # Also add pip to the list of dependencies
                    SOFTWARE_DEPENDENCIES+=("pip")
                fi
            done < ${DEPENDENCY_FILEPATH}
        done
        # Remove duplicates
        SOFTWARE_DEPENDENCIES=($(printf "%s\n" "${SOFTWARE_DEPENDENCIES[@]}" | sort -u))

        # Count how many missing dependencies
        MISSING=0
        for dep in "${SOFTWARE_DEPENDENCIES[@]}"; do
            # Check for dependency
            if ! command -v $dep &> /dev/null ; then
                if $VERBOSE; then
                    println "\t[--] %s\n" "$dep"
                fi
                MISSING=$(($MISSING + 1))
            else
                if $VERBOSE; then
                    println "\t[OK] %s\n" "$dep"
                fi
            fi
        done

        # In case Python is a requirement
        if ${PYTHON_FOUND}; then
            # In case the requirements.txt file exists
            PYTHON_REQUIREMENTS="$(dirname ${DEPENDENCY_FILEPATHS[0]})/requirements.txt"
            if [[ -f ${PYTHON_REQUIREMENTS} ]]; then
                # Check for python modules
                if $VERBOSE; then
                    println "\nChecking for Python requirements\n"
                    println "\tThis will automatically install missing dependencies with pip\n"
                    read -p $'\tDo you want to continue? [Y/n] ' -n 1 -r REPLY
                    println "\n"
                    if [[ $REPLY =~ ^Y$ ]]; then
                        python3 -m pip install -r ${PYTHON_REQUIREMENTS} --ignore-installed
                    fi
                fi
            fi
        fi
        
        # Return 1 in case of missing dependencies
        if [ "$MISSING" -gt "0" ]; then
            if $VERBOSE; then
                println "\nPlease, install all the missing dependencies and try again.\n\n"
            fi
            return 1
        fi
        
        # Return 0 if all external software dependencies are satisfied
        if $VERBOSE; then
            println "\nAll required dependencies satisfied!\n\n"
        fi
        return 0
    fi

    # In case of no dependencies
    if $VERBOSE; then
        println "\nNo dependencies found!\n\n"
    fi
    return 1
}

# Check for new software release
check_for_software_updates () {
    LAST_SOFTWARE_VERSION="$(get_last_software_release)"
    VERSION=$1
    if [ ! -z ${LAST_SOFTWARE_VERSION} ]; then
        if [[ "${LAST_SOFTWARE_VERSION}" != $VERSION ]]; then
            println "A new version of meta-index is available!\n" "${LAST_SOFTWARE_VERSION}"
            println "https://github.com/BlankenbergLab/meta-index/releases\n\n"
        fi
    fi
}

# Print citations
citations () {
    println "TBA\n\n"
}

# Print credits
credits () {
    println "Thanks for using meta-index!\n\n"
    println "Please credit this tool in your manuscript by citing:\n\n"
    citations

    println "Remember to star the meta-index repository on GitHub to stay updated on its development and new features:\n"
    println "https://github.com/BlankenbergLab/meta-index\n\n"
}

# Format time in human-readable format
# Credits: https://unix.stackexchange.com/a/27014
displaytime () {
    ESTIMATED_DAYS=$(bc <<< "${1}/60/60/24")  # Days
    ESTIMATED_HOURS=$(bc <<< "${1}/60/60%24") # Hours
    ESTIMATED_MINUTES=$(bc <<< "${1}/60%60")  # Minutes
    ESTIMATED_SECONDS=$(bc <<< "${1}%60")     # Seconds
    # Print elapsed days
    if [[ "$(echo "${ESTIMATED_DAYS} > 0" | bc)" -eq "1" ]]; then println "%s days " "${ESTIMATED_DAYS}"; fi
    # Print elapsed hours
    if [[ "$(echo "${ESTIMATED_HOURS} > 0" | bc)" -eq "1" ]]; then println "%s hours " "${ESTIMATED_HOURS}"; fi
    # Print elapsed minutes
    if [[ "$(echo "${ESTIMATED_MINUTES} > 0" | bc)" -eq "1" ]]; then println "%s minutes " "${ESTIMATED_MINUTES}"; fi
    # Print "and" in case of multiple days, hours, or minutes have been printed before seconds
    if [[ "$(echo "${ESTIMATED_DAYS} > 0" | bc)" -eq "1" ]] || \
       [[ "$(echo "${ESTIMATED_HOURS} > 0" | bc)" -eq "1" ]] || \
       [[ "$(echo "${ESTIMATED_MINUTES} > 0" | bc)" -eq "1" ]]; then 
       println "and "
    fi
    # Finally print elapsed seconds
    println "%.0f seconds\n" "${ESTIMATED_SECONDS}"
}

# Download genomes from NCBI GenBank
esearch_txid () {
    DBDIR=$1            # Path to the database root directory
    HOW_MANY=$2         # Limit the number of genomes per species (-1 by default, unlimited)
    TAX_ID=$3           # NCBI taxonomy ID
    FULL_TAXONOMY=$4    # Full NCBI taxonomy
    GENOME_CATEGORY=$5  # Genome category ("references" or "mags")

    # Define a search criteria which is empty by default for MAGs
    SEARCH_CRITERIA=""
    if [[ "${GENOME_CATEGORY}" = "references" ]]; then
        # Retrieve genomes excluded from RefSeq in case of references
        SEARCH_CRITERIA="NOT excluded-from-refseq [PROP]"
    fi

    # Stop downloading if the number of genomes exceeds the limit
    COUNT_GENOMES=0
    # Download GCAs associated to a specific tax_id
    esearch -db assembly -query "txid${TAX_ID} ${SEARCH_CRITERIA}" < /dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r URL; do
            # Check whether the number of downloaded genomes exceeded the limit
            if [[ "${COUNT_GENOMES}" -ge "${HOW_MANY}" ]] && [[ "${HOW_MANY}" -gt "0" ]]; then
                # Stop downloading genomes
                break
            fi
            # Create a directory for the current taxonomy
            TAXDIR=$DBDIR/$(echo "${FULL_TAXONOMY}" | sed 's/|/\//g')
            OUTDIR=$TAXDIR/genomes
            mkdir -p $OUTDIR
            # Download GCA
            FNAME=$(echo $URL | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
            GCA=${FNAME%%.*}
            if [[ ! -f "$OUTDIR/${GCA}.fna.gz" ]]; then
                wget -q "$URL/$FNAME" -O $OUTDIR/${GCA}.fna.gz
                if [[ -f "$OUTDIR/${GCA}.fna.gz" ]]; then
                    if gzip -t $OUTDIR/${GCA}.fna.gz; then
                        # Define a file of file (fof) with the list of genomes for current species
                        GCAPATH=$(readlink -m $OUTDIR/${GCA}.fna.gz)
                        println "%s : %s\n" "$GCA" "$GCAPATH" >> $TAXDIR/genomes.fof
                        println "%s\n" "$GCAPATH" >> $TAXDIR/genomes.txt
                        println "%s\n" "$GCA" >> $TAXDIR/${GENOME_CATEGORY}.txt

                        # Increment the genome counter
                        COUNT_GENOMES=$((COUNT_GENOMES + 1))
                    else
                        # Delete corrupted genome
                        rm $OUTDIR/${GCA}.fna.gz
                        println "\t[ERROR][ID=%s][TAXID=%s] Corrupted genome\n" "$GCA" "${TAX_ID}"
                    fi
                else
                    # Report missing genomes
                    println "\t[ERROR][ID=%s][TAXID=%s] Unable to download genome\n" "$GCA" "${TAX_ID}"
                fi
            fi
          done
}

# Get the last software version
# Credits: https://gist.github.com/lukechilds/a83e1d7127b78fef38c2914c4ececc3c
get_last_software_release () {
    curl --silent "https://api.github.com/repos/BlankenbergLab/meta-index/releases/latest" | \
        grep '"tag_name":' | \
        sed -E 's/.*"([^"]+)".*/\1/'
}

# Wrapper for the howdesbt pipeline
# This is applied on all the taxonomic levels
howdesbt_wrapper () {
    LEVEL_DIR=$(readlink -m $1)     # Taxonomic level folder
    FILTER_SIZE=$2                  # Bloom filter size
    KMER_LEN=$3                     # Kmer length
    NPROC=$4                        # Max nproc for multiprocessing
    
    # Retrieve the name of the current taxonomic level
    LEVEL_NAME=$(basename ${LEVEL_DIR})
    # Define index folder path
    INDEX_DIR=${LEVEL_DIR}/index

    # Remove old list of bloom filter files and
    # the bloom filter representation of the current level
    rm -f ${LEVEL_DIR}/${LEVEL_NAME}.txt
    rm -f ${LEVEL_DIR}/${LEVEL_NAME}.bf
    # Remove the old log file
    rm -f ${LEVEL_DIR}/howdesbt.log
    # Also remove the old index
    rm -rf ${INDEX_DIR}

    if [[ "${LEVEL_NAME:0:1}" = "s" ]]; then
        # In case of species levels
        # Create the filters folders for taking track of the bloom filter files
        mkdir -p ${LEVEL_DIR}/filters
        # Build the bloom filter files
        find ${LEVEL_DIR}/genomes -type f -iname "*.gz" -follow | xargs -n 1 -I {} bash -c \
            'GENOME_FILEPATH={};
             GENOME_NAME=$(basename ${GENOME_FILEPATH});
             GENOME_FORMAT="${GENOME_NAME##*.}";
             GENOME_NAME="${GENOME_NAME%.*}";
             if [[ "${GENOME_FILEPATH}" = *.gz ]]; then
                GENOME_FORMAT="${GENOME_NAME##*.}";
                GENOME_NAME="${GENOME_NAME%.*}";
             fi

             if [[ ! -f '"${LEVEL_DIR}"'/filters/${GENOME_NAME}.bf ]]; then
                gzip -dc ${GENOME_FILEPATH} > '"${LEVEL_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT};
                howdesbt makebf --k='"${KMER_LEN}"' --min=2 --bits='"${FILTER_SIZE}"' --hashes=1 --seed=0,0 \
                                '"${LEVEL_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT} \
                                --out='"${LEVEL_DIR}"'/filters/${GENOME_NAME}.bf \
                                --threads='"$NPROC"' \
                                >> '"${LEVEL_DIR}"'/howdesbt.log 2>&1;
                rm '"${LEVEL_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT};
             fi

             readlink -m '"${LEVEL_DIR}"'/filters/${GENOME_NAME}.bf >> '"${LEVEL_DIR}"'/'"${LEVEL_NAME}"'.txt'
    else
        # For all the other taxonomic levels
        # Define the new list of bloom filters
        for DIRECTORY in ${LEVEL_DIR}/*/; do
            UPPER_LEVEL=$(basename $DIRECTORY)
            if [[ -f ${LEVEL_DIR}/${UPPER_LEVEL}/${UPPER_LEVEL}.bf ]]; then
                # Each file is the result of the OR operator on all representative bloom filters in the upper taxonomic level
                echo "$(readlink -m ${LEVEL_DIR}/${UPPER_LEVEL}/${UPPER_LEVEL}.bf)" >> ${LEVEL_DIR}/${LEVEL_NAME}.txt
            fi
        done
    fi

    # Create the index folder
    mkdir -p ${INDEX_DIR}
    # Move to the index folder
    # This will force howdesbt to build the compressed nodes into the index folder
    cd ${INDEX_DIR}
    
    # Count how many elements must be clustered
    HOWMANY=$(wc -l ${LEVEL_DIR}/${LEVEL_NAME}.txt | cut -d" " -f1)
    if [[ "$HOWMANY" -gt "1" ]]; then
        # Create a tree topology file with howdesbt
        howdesbt cluster --list=${LEVEL_DIR}/${LEVEL_NAME}.txt \
                         --bits=${FILTER_SIZE} \
                         --tree=${INDEX_DIR}/union.sbt \
                         --nodename=${INDEX_DIR}/node{number} \
                         --keepallnodes \
                         >> ${LEVEL_DIR}/howdesbt.log 2>&1
        # Build the bloom filter files for the tree
        howdesbt build --howde \
                       --tree=${INDEX_DIR}/union.sbt \
                       --outtree=${INDEX_DIR}/index.detbrief.sbt \
                       >> ${LEVEL_DIR}/howdesbt.log 2>&1
        # Remove the union.sbt file
        rm -f ${INDEX_DIR}/union.sbt
        
        # Fix node paths in the final index.detbrief.sbt file
        while read node; do
            # The number of stars defines the depth of the current node in the tree
            # Count how many stars occur before the node name
            STARS="$(grep -o '\*' <<<"$node" | grep -c .)"
            # Remove all the start from the current node
            NODE_NAME="${node:$STARS}"
            # Add the full path
            NODE_PATH="${INDEX_DIR}/${NODE_NAME}"
            # Add the depth again
            NODE_DEPTH="${node:0:$STARS}${NODE_PATH}"
            # Define the new index.detbrief.sbt
            echo "${NODE_DEPTH}" >> ${INDEX_DIR}/index.full.detbrief.sbt
        done < ${INDEX_DIR}/index.detbrief.sbt
        # Use the fixed index.detbrief.sbt
        mv ${INDEX_DIR}/index.full.detbrief.sbt ${INDEX_DIR}/index.detbrief.sbt

        # Merge all the leaves together by applying the OR logical operator on the bloom filter files
        # The resulting bloom filter is the representative one, which is the same as the root node of the tree
        while read BFPATH; do
            if [[ ! -f ${LEVEL_DIR}/${LEVEL_NAME}.bf ]]; then
                cp $BFPATH ${LEVEL_DIR}/${LEVEL_NAME}.bf
            else
                # Merge bloom filter files applying the OR logical operator
                howdesbt bfoperate $BFPATH ${LEVEL_DIR}/${LEVEL_NAME}.bf --or --out=${LEVEL_DIR}/merged.bf \
                    >> ${LEVEL_DIR}/howdesbt.log 2>&1
                # Rename the resulting file
                mv ${LEVEL_DIR}/merged.bf ${LEVEL_DIR}/${LEVEL_NAME}.bf
            fi
        done < ${LEVEL_DIR}/${LEVEL_NAME}.txt
    else
        # With only one bloom filter, it does not make sense to build an index with howdesbt
        BFPATH=$(head -n1 "${LEVEL_DIR}/${LEVEL_NAME}.txt")
        cp $BFPATH ${LEVEL_DIR}/${LEVEL_NAME}.bf
        # Manually define the union.sbt file with the single node
        echo "$BFPATH" > ${INDEX_DIR}/union.sbt
        # Build the RRR compressed bloom filter file for the node
        howdesbt build --howde \
                       --tree=${INDEX_DIR}/union.sbt \
                       --outtree=${INDEX_DIR}/index.detbrief.sbt \
                       > ${LEVEL_DIR}/howdesbt.log 2>&1
        # Remove the union.sbt file
        rm -f ${INDEX_DIR}/union.sbt
    fi
}
# Export howdesbt_wrapper to sub-shells
export -f howdesbt_wrapper

# Build a kmers presence/absence matrix with kmtricks
kmtricks_matrix_wrapper () {
    GENOMES_FOF=$1  # Path to the input fof file
    RUN_DIR=$2      # Path to the kmtricks working directory
    NPROC=$3        # Number of parallel processes
    OUT_TABLE=$4    # Path to the output kmers table
    WORKDIR=$5      # Path to the working directory

    # Build kmers presence/absence matrix
    kmtricks pipeline --file ${GENOMES_FOF} \
                      --run-dir ${RUN_DIR} \
                      --mode kmer:pa:bin \
                      --hard-min 1 \
                      --cpr \
                      --threads $NPROC \
                      >> $WORKDIR/kmtricks.log 2>&1
    
    # Aggregate sub-matrices
    kmtricks aggregate --run-dir ${RUN_DIR} \
                       --matrix kmer \
                       --format text \
                       --cpr-in \
                       --sorted \
                       --threads $NPROC \
                       --output ${OUT_TABLE} \
                       >> $WORKDIR/kmtricks.log 2>&1
}
# Export kmtricks_matrix_wrapper to sub-shells
export -f kmtricks_matrix_wrapper

# Print line on standard output and write in log file
LOG_FILEPATH=""
println () {
    LINE=$1
    VALUES=("${@:2}")
    if [[ ! -z ${LOG_FILEPATH} ]]; then
        printf "$LINE" "${VALUES[@]}" >> ${LOG_FILEPATH}
    fi
    printf "$LINE" "${VALUES[@]}"
}
# Export println to sub-shells
export -f println

# Run CheckM for quality control base on completeness and contamination
# Remember to create a global variable "CHECKMTABLES" before calling this function 
# for accessing the result as a list of file paths separated by comma
run_checkm () {
    INLIST=$1       # Path to the file with the list of paths to the input genomes
    RUNID=$2        # Unique ID of the current CheckM run
    EXTENSION=$3    # Extension of input genome files
    TMPDIR=$4       # Temporary folder
    NPROC=$5        # Number of parallel processes for CheckM

    # Run CheckM
    println "Running CheckM on %s\n" "$INLIST"
    CHECKM_START_TIME="$(date +%s.%3N)"
    # Define temporary run directory
    RUNDIR=$TMPDIR/checkm/$RUNID
    mkdir -p $RUNDIR/tmp
    # Split the set of bins in chunks with 1000 genomes at most
    println "\tOrganising genomes in chunks\n"
    split --numeric-suffixes=1 --lines=1000 --suffix-length=3 --additional-suffix=.txt $INLIST $RUNDIR/tmp/bins_
    CHUNKS=`ls "$TMPDIR"/checkm/"$RUNID"/tmp/bins_*.txt 2>/dev/null | wc -l`
    CHECKMTABLES=""
    for filepath in $RUNDIR/tmp/bins_*.txt; do
        # Retrieve chunk id
        filename="$(basename $filepath)"
        suffix="${filename#*_}"
        suffix="${suffix%.txt}"
        println "\tProcessing chunk %s/%s\n" "$((10#$suffix))" "$CHUNKS"
        # Create chunk folder
        mkdir -p $RUNDIR/tmp/bins_${suffix}
        for bin in `sed '/^$/d' $filepath`; do
            # Make a symbolic link to the 1000 genomes for the current chunk
            ln -s $bin $RUNDIR/tmp/bins_${suffix}
        done
        # Create the CheckM folder for the current chunk
        mkdir -p $RUNDIR/run_$suffix/
        # Run CheckM
        checkm lineage_wf -t $NPROC \
                          -x $EXTENSION \
                          --pplacer_threads $NPROC \
                          --tab_table -f $RUNDIR/run_${suffix}.tsv \
                          $RUNDIR/tmp/bins_${suffix} $RUNDIR/run_$suffix \
                          > $RUNDIR/checkm.log 2>&1
        # Take trace of the CheckM output tables
        CHECKMTABLES=$RUNDIR/run_${suffix}.tsv,$CHECKMTABLES
    done
    CHECKM_END_TIME="$(date +%s.%3N)"
    CHECKM_ELAPSED="$(bc <<< "${CHECKM_END_TIME}-${CHECKM_START_TIME}")"
    println "\tTotal elapsed time: %s\n\n" "$(displaytime ${CHECKM_ELAPSED})"
    # Trim the last comma out of the list of CheckM output table file paths
    CHECKMTABLES="${CHECKMTABLES%?}"
}

# Print a standard error message and redirect the user to the Issues and Discussions pages on GitHub
standard_error_message () {
    println "If you think this is a bug and need support, please open an Issue or a new Discussion on the official GitHub repository.\n"
    println "We would be happy to answer your questions and help you troubleshoot any kind of issues with the meta-index framework.\n\n"
    println "https://github.com/BlankenbergLab/meta-index/issues\n"
    println "https://github.com/BlankenbergLab/meta-index/discussions\n\n"
}

# Print the list of supported file extensions for input genomes
supported_extensions () {
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ ! "$line" =~ ^#.* ]]; then
            println '%s\n' "$line"
        fi
    done < $1
}