#!/bin/bash
#title          :minimizers
#description    :Extract minimizers to fasta format
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#====================================================

INPUTS_DIR=$1           # Folder path with the set of input genome files
OUTPUT_DIR=$2           # Output folder with the new fasta files containing minimizers
EXTENSION=$3            # e.g. fna (assume it has been uniformed with uniform_inputs.sh)
FIND_NPROC=$4           # Process input files in parallel

MINIMIZERS_SIZE=$5      # Minimizers size
MINIMIZERS_WINDOW=$6    # Minimizers window
MINIMIZERS_NPROC=$4     # Process fasta records in parallel

find ${INPUTS_DIR} \
    -type f -iname "*.${EXTENSION}" -follow | xargs -n 1 -P ${FIND_NPROC} -I {} bash -c \
        'INPUT={}; \
         minimizers --input "$INPUT" \
                    --output '"$OUTPUT_DIR"'/$(basename "$INPUT") \
                    --output-type fasta \
                    --aggregate \
                    --size '"${MINIMIZERS_SIZE}"' \
                    --window '"${MINIMIZERS_WINDOW}"' \
                    --nproc '"${MINIMIZERS_NPROC}"';'