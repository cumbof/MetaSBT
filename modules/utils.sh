#!/bin/bash
#title          :utils
#description    :Some utility functions
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="Jun 5, 2022"
VERSION="0.1.0"

# Check for external software dependencies
check_dependencies () {
    VERBOSE=$1
    if $VERBOSE; then
        println "Checking for software dependencies\n"
    fi
    # Define the set of dependencies
    # esearch, esummary, and xtract come from the Entrez Direct (EDirect) utility
    DEPENDENCIES=("bc" "checkm" "esearch" "esummary" "gzip" "howdesbt" "kmtricks" "ncbitax2lin" "ntcard" "wget" "xtract")

    # Count how many missing dependencies
    MISSING=0
    for dep in "${DEPENDENCIES[@]}"; do
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
# This is applied on all the taxonomic level except the species
howdesbt_wrapper () {
    LEVEL_DIR=$(readlink -m $1)     # Taxonomic level folder
    FILTER_SIZE=$2                  # Bloom filter size
    
    # Retrieve the name of the current taxonomic level
    LEVEL_NAME=$(basename ${LEVEL_DIR})
    # Define index folder path
    INDEX_DIR=${LEVEL_DIR}/index

    # Remove old list of bloom filter files and
    # the bloom filter representation of the current level
    rm -f ${LEVEL_DIR}/${LEVEL_NAME}.txt
    rm -f ${LEVEL_DIR}/${LEVEL_NAME}.bf
    # Also remove the old index
    rm -rf ${INDEX_DIR}

    # Define the new list of bloom filters
    for DIRECTORY in ${LEVEL_DIR}/*/; do
        UPPER_LEVEL=$(basename $DIRECTORY)
        # Each file is the result of the OR operator on all representative bloom filters in the upper taxonomic level
        echo "$(readlink -m ${LEVEL_DIR}/${UPPER_LEVEL}/${UPPER_LEVEL}.bf)" >> ${LEVEL_DIR}/${LEVEL_NAME}.txt
    done

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
                         > ${LEVEL_DIR}/howdesbt.log 2>&1
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

# Wrapper for the kmtricks pipeline
# This is applied at the species level only
kmtricks_index_wrapper () {
    INPUT=$1        # Path to the genomes.fof file for a specific species
    DBDIR=$2        # Database root directory path
    KMER_LEN=$3     # Length of the kmers
    FILTER_SIZE=$4  # Bloom filter size
    NPROC=$5        # Max nproc for multiprocessing
    FOLDERPATH=$(dirname "$INPUT")
    if [[ ! -f "$FOLDERPATH/index/kmtricks.fof" ]]; then
        # Run the kmtricks pipeline
        kmtricks pipeline --file $INPUT \
                          --run-dir $FOLDERPATH/index \
                          --kmer-size ${KMER_LEN} \
                          --mode hash:bft:bin \
                          --hard-min 1 \
                          --bloom-size ${FILTER_SIZE} \
                          --bf-format howdesbt \
                          --cpr \
                          --skip-merge \
                          --threads $NPROC \
                          > $FOLDERPATH/index/kmtricks.log 2>&1
        
        FOLDERNAME="${FOLDERPATH##*/}"
        HOWMANY=$(wc -l $INPUT | cut -d" " -f1)
        if [[ "$HOWMANY" -gt "1" ]]; then
            # Index genomes by building a sequence bloom tree with howdesbt
            kmtricks index --run-dir $FOLDERPATH/index \
                           --howde \
                           --threads $NPROC
            # Merge all the leaves together by applying the OR logical operator on the bloom filter files
            # The resulting bloom filter is the representative one, which is the same as the root node of the tree
            while read line; do
                GENOME=$(echo "$line" | cut -d" " -f1)
                if [[ ! -f $FOLDERPATH/${FOLDERNAME}.bf ]]; then
                    cp $FOLDERPATH/index/filters/${GENOME}.bf $FOLDERPATH/${FOLDERNAME}.bf
                else
                    # Merge bloom filter files applying the OR logical operator
                    howdesbt bfoperate $FOLDERPATH/index/filters/${GENOME}.bf $FOLDERPATH/${FOLDERNAME}.bf --or --out=$FOLDERPATH/merged.bf \
                        >> $FOLDERPATH/index/kmtricks.log 2>&1
                    # Rename the resulting file
                    mv $FOLDERPATH/merged.bf $FOLDERPATH/${FOLDERNAME}.bf
                fi
            done < $INPUT
        else
            # With only one genome, it does not make sense to build an index with howdesbt
            mkdir -p $FOLDERPATH/index/howde_index
            cp $INPUT $FOLDERPATH/index/kmtricks.fof
            GENOME=$(head -n1 "$INPUT" | cut -d" " -f1)
            cp $FOLDERPATH/index/filters/${GENOME}.bf $FOLDERPATH/${FOLDERNAME}.bf
        fi
    fi
}
# Export kmtricks_index_wrapper to sub-shells
export -f kmtricks_index_wrapper

# Build a kmers matrix with kmtricks
kmtricks_matrix_wrapper () {
    GENOMES_FOF=$1
    RUN_DIR=$2
    NPROC=$3
    OUT_TABLE=$4

    # Build matrix
    kmtricks pipeline --file ${GENOMES_FOF} \
                      --run-dir ${RUN_DIR} \
                      --mode kmer:count:bin \
                      --hard-min 1 \
                      --cpr \
                      --threads $NPROC \
                      >> ${RUN_DIR}/kmtricks.log 2>&1
    # Aggregate
    kmtricks aggregate --run-dir ${RUN_DIR} \
                       --matrix kmer \
                       --format text \
                       --cpr-in \
                       --sorted \
                       --threads $NPROC \
                       --output ${OUT_TABLE} \
                       >> ${RUN_DIR}/kmtricks.log 2>&1
}

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

# Transpose matrix file
# Credits: https://stackoverflow.com/a/28167793
transpose () {
  awk '{for (i=1; i<=NF; i++) a[i,NR]=$i; max=(max<NF?NF:max)}
        END {for (i=1; i<=max; i++)
              {for (j=1; j<=NR; j++) 
                  printf "%s%s", a[i,j], (j<NR?OFS:ORS)
              }
        }'
}