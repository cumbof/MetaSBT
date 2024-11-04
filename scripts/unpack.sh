#!/bin/bash
#title          :unpack
#description    :Unack a MetaSBT database to a specific location
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===============================================================

INPUT_TARBALL=$1  # Path to the input compressed tarball with a MetaSBT database
OUTPUT_DIR=$2     # Path to the output folder where a MetaSBT database will be installed in

NPROC=1

# Create a Gzip compressed tarball
tar -xzvf ${INPUT_TARBALL} -C ${OUTPUT_DIR}

# Fix the SBT nodes file paths in all the SBT definition files
find ${OUTPUT_DIR}/$(basename ${INPUT_TARBALL} .tar.gz) -type f '(' -iname "*.sbt" -o \
                                                                    -iname "*.txt" \
                                                                ')' -follow | xargs -n 1 -P $NPROC -I {} bash -c \
    'INPUT={}; \
     sed -i "s/\/.*\/'$(basename ${INPUT_TARBALL} .tar.gz)'/'${OUTPUT_DIR//\//\\\/}'\/'$(basename ${INPUT_TARBALL} .tar.gz)'/g" $INPUT'
