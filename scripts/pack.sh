#!/bin/bash
#title          :pack
#description    :Pack a MetaSBT database as a Gzip compressed tarball
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#====================================================================

DB_DIR=$1          # Path to the MetaSBT database root folder
OUTPUT_TARBALL=$2  # Path to the output compressed tarball

NPROC=1

if [[ -d "${DB_DIR}}" && "${DB_DIR}" == */ ]]; then
    DB_DIR="${DB_DIR%/*}"
fi

# Exclude genomes, temporary data, and logs from the tarball
find ${DB_DIR} -type d '(' -iname "genomes" -o \
                           -iname "k__*" -o \
                           -iname "p__*" -o \
                           -iname "c__*" -o \
                           -iname "o__*" -o \
                           -iname "f__*" -o \
                           -iname "g__*" -o \
                           -iname "s__*" -o \
                           -iname "strains" \
                       ')' -follow | xargs -n 1 -P $NPROC -I {} bash -c \
    'FOLDER={}; \
     if [[ $(basename $FOLDER) == "genomes" ]]; then \
        echo "*" > $(realpath $FOLDER)/.tagignore;
     elif [[ $(basename $FOLDER) == "strains" ]]; then \
        echo "tmp" > $(realpath $FOLDER)/.tagignore; \
        echo "*.log" >> $(realpath $FOLDER)/.tagignore; \
     else \
        echo "*.log" >> $(realpath $FOLDER)/.tagignore; \
     fi'

echo "*.log" > ${DB_DIR}/.tagignore
echo "*.sh" >> ${DB_DIR}/.tagignore

# Create a Gzip compressed tarball
tar -czvf ${OUTPUT_TARBALL} --exclude-ignore=.tagignore -C $(dirname ${DB_DIR}) $(basename ${DB_DIR})

# Get rid of the .tagignore files
find ${DB_DIR} -type f -iname ".tagignore" -exec rm {} \;
