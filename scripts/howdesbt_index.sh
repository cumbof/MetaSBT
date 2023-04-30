#!/bin/bash
#title          :howdesbt_index
#description    :Index genomes with HowDeSBT
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================================================

DBDIR=$1        # Path to the database root folder
FILTER_SIZE=$2  # Bloom filter size
KMER_LEN=$3     # Length of the kmers
NPROC=$4        # Run HowDeSBT in parallel

# Wrapper for the howdesbt pipeline
# It assumes that the genomes have been already organised in the "genomes" folder under the database root directory
howdesbt_wrapper () {
    DB_DIR=$(readlink -m $1)    # Database root folder
    FILTER_SIZE=$2              # Bloom filter size
    KMER_LEN=$3                 # Kmer length
    NPROC=$4                    # Max nproc for multiprocessing
    
    # Define index folder path
    INDEX_DIR=${DB_DIR}/index

    # Remove the old log file
    rm -f ${DB_DIR}/howdesbt.log
    # Also remove the old index
    rm -rf ${INDEX_DIR}

    # Create the filters folders for taking track of the bloom filter files
    mkdir -p ${DB_DIR}/filters

    # Build the bloom filter files
    find ${DB_DIR}/genomes -type f -iname "*.gz" -follow | xargs -n 1 -I {} bash -c \
        'GENOME_FILEPATH={};
         GENOME_NAME=$(basename ${GENOME_FILEPATH});
         GENOME_FORMAT="${GENOME_NAME##*.}";
         GENOME_NAME="${GENOME_NAME%.*}";
         if [[ "${GENOME_FILEPATH}" = *.gz ]]; then
            GENOME_FORMAT="${GENOME_NAME##*.}";
            GENOME_NAME="${GENOME_NAME%.*}";
         fi
         if [[ ! -f '"${DB_DIR}"'/filters/${GENOME_NAME}.bf ]]; then
            if [[ -f '"${DB_DIR}"'/filters/${GENOME_NAME}.bf.gz ]]; then
                gunzip '"${DB_DIR}"'/filters/${GENOME_NAME}.bf.gz ;
            else
                gzip -dc ${GENOME_FILEPATH} > '"${DB_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT};
                howdesbt makebf --k='"${KMER_LEN}"' --min=2 --bits='"${FILTER_SIZE}"' --hashes=1 --seed=0,0 \
                                '"${DB_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT} \
                                --out='"${DB_DIR}"'/filters/${GENOME_NAME}.bf \
                                --threads='"$NPROC"' \
                                >> '"${DB_DIR}"'/howdesbt.log 2>&1;
                rm '"${DB_DIR}"'/genomes/${GENOME_NAME}.${GENOME_FORMAT};
            fi
         fi
         readlink -m '"${DB_DIR}"'/filters/${GENOME_NAME}.bf >> '"${DB_DIR}"'/db.txt'

    # Create the index folder
    mkdir -p ${INDEX_DIR}
    # Move to the index folder
    # This will force howdesbt to build the compressed nodes into the index folder
    cd ${INDEX_DIR}
    
    # Count how many elements must be clustered
    HOWMANY=$(wc -l ${DB_DIR}/db.txt | cut -d" " -f1)
    if [[ "$HOWMANY" -gt "1" ]]; then
        # Create a tree topology file with howdesbt
        howdesbt cluster --list=${DB_DIR}/db.txt \
                         --bits=${FILTER_SIZE} \
                         --tree=${INDEX_DIR}/union.sbt \
                         --nodename=${INDEX_DIR}/node{number} \
                         --keepallnodes \
                         >> ${DB_DIR}/howdesbt.log 2>&1
        # Build the bloom filter files for the tree
        howdesbt build --howde \
                       --tree=${INDEX_DIR}/union.sbt \
                       --outtree=${INDEX_DIR}/index.detbrief.sbt \
                       >> ${DB_DIR}/howdesbt.log 2>&1
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
        howdesbt bfoperate --list=${DB_DIR}/db.txt --or --out=${DB_DIR}/db.bf >> ${DB_DIR}/howdesbt.log 2>&1
    else
        # With only one bloom filter, it does not make sense to build an index with howdesbt
        BFPATH=$(head -n1 "${DB_DIR}/db.txt")
        cp $BFPATH ${DB_DIR}/db.bf
        # Manually define the union.sbt file with the single node
        echo "$BFPATH" > ${INDEX_DIR}/union.sbt
        # Build the RRR compressed bloom filter file for the node
        howdesbt build --howde \
                       --tree=${INDEX_DIR}/union.sbt \
                       --outtree=${INDEX_DIR}/index.detbrief.sbt \
                       >> ${DB_DIR}/howdesbt.log 2>&1
        # Remove the union.sbt file
        rm -f ${INDEX_DIR}/union.sbt
    fi

    # Compress filters
    find ${DB_DIR}/filters -type f -iname "*.bf" -follow -exec gzip {} \;
}
# Export howdesbt_wrapper to sub-shells
export -f howdesbt_wrapper

# Retrieve the absolute path of the database folder
DBDIR=$(readlink -m $DBDIR)

# Take track of the directory from which this module is launched
# This is required in order to come back on this folder after running HowDeSBT
EXEC_DIR="$PWD"

# Also build the index for the kingdom
howdesbt_wrapper $DBDIR ${FILTER_SIZE} ${KMER_LEN} $NPROC

# HowDeSBT calls automatically change the current folder to the taxonomy directory
# Come back to the folder from which this module has been launched
cd ${EXEC_DIR}