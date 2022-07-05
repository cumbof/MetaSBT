#!/bin/bash
#title          :howdesbt_index
#description    :Index genomes under all the taxonomic lavels in a specific database
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================================================

DBDIR=$1        # Path to the database root folder
FILTER_SIZE=$2  # Bloom filter size
KMER_LEN=$3     # Length of the kmers
NPROC=$4        # Run HowDeSBT in parallel
XARGS_NPROC=$5  # Index clusters in parallel

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
        howdesbt bfoperate --list=${LEVEL_DIR}/${LEVEL_NAME}.txt --or --out=${LEVEL_DIR}/${LEVEL_NAME}.bf >> ${LEVEL_DIR}/howdesbt.log 2>&1
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
                       >> ${LEVEL_DIR}/howdesbt.log 2>&1
        # Remove the union.sbt file
        rm -f ${INDEX_DIR}/union.sbt
    fi
}
# Export howdesbt_wrapper to sub-shells
export -f howdesbt_wrapper

# Retrieve the absolute path of the database folder
DBDIR=$(readlink -m $DBDIR)

# Take track of the directory from which this module is launched
# This is required in order to come back on this folder after running HowDeSBT
EXEC_DIR="$PWD"

# Start with species at depth 7 up to the kingdom level
DEPTH=7
for LEVELNAME in "species" "genus" "family" "order" "class" "phylum" "kingdom"; do
    # Retrieve the level ID from the level name
    LEVEL="${LEVELNAME:0:1}__"

    # Process all the other taxonomic levels with howdesbt
    printf "Running howdesbt at the %s level\n" "$LEVELNAME"
    find $DBDIR -maxdepth $DEPTH -type d -iname "${LEVEL}*" -follow | xargs -n 1 -P ${XARGS_NPROC} -I {} bash -c \
        'LEVELDIR={};
         howdesbt_wrapper $LEVELDIR '"${FILTER_SIZE}"' '"${KMER_LEN}"' '"$NPROC"';'
    # Decrease the depth while moving to a higher taxonomic level
    DEPTH=$(expr $DEPTH - 1)
done

printf "Building the bloom filter root node of the database\n"
# Also build the index for the kingdom
howdesbt_wrapper $DBDIR ${FILTER_SIZE} ${KMER_LEN} $NPROC

# HowDeSBT calls automatically change the current folder to the taxonomy directory
# Come back to the folder from which this module has been launched
cd ${EXEC_DIR}