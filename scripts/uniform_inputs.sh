#!/bin/bash
#title          :uniform_inputs
#description    :Uniform input genome files extension
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#====================================================

INPUTS_DIR=$1           # Folder path with the set of input genome files
CURRENT_EXTENSION=$2    # e.g. fa
NEW_EXTENSION=$3        # e.g. fna

find ${INPUTS_DIR} \
    -type f -iname "*.${CURRENT_EXTENSION}" -follow | xargs -n 1 -I {} bash -c \
        'INPUT={}; \
         mv "$INPUT" "${INPUT%.'"${CURRENT_EXTENSION}"'}.'"${NEW_EXTENSION}"'";'