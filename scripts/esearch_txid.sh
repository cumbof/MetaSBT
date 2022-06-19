#!/bin/bash
#title          :esearch_txid
#description    :Retrieve GCAs from NCBI GenBank given a specific taxonomic ID 
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#=============================================================================

TXID=$1         # NCBI taxonomy ID
CATEGORY=$2     # Genome category ("references" or "mags")
FOLDERPATH=$3   # Path to the database root directory

# Create the target folder if it does not exist
mkdir -p $FOLDERPATH

# Define a search criteria which is empty by default for MAGs
SEARCH_CRITERIA=""
if [[ "$CATEGORY" = "references" ]]; then
    # Retrieve genomes excluded from RefSeq in case of references
    SEARCH_CRITERIA="NOT excluded-from-refseq [PROP]"
fi

# Download GCAs associated to a specific txid
esearch -db assembly -query "txid$TXID ${SEARCH_CRITERIA}" < /dev/null \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r URL; do
        # Create a directory for the current taxonomy
        OUTDIR=$FOLDERPATH/genomes
        mkdir -p $OUTDIR
        # Download GCA
        FNAME=$(echo $URL | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
        GCA=${FNAME%%.*}
        if [[ ! -f "$OUTDIR/${GCA}.fna.gz" ]]; then
            wget -q "$URL/$FNAME" -O $OUTDIR/${GCA}.fna.gz
            if [[ -f "$OUTDIR/${GCA}.fna.gz" ]]; then
                if gzip -t $OUTDIR/${GCA}.fna.gz; then
                    printf "[TXID=%s][ID=%s] Ok\n" "$TXID" "$GCA"
                    # Define a file of file (fof) with the list of genomes for current species
                    GCAPATH=$(readlink -m $OUTDIR/${GCA}.fna.gz)
                    println "%s\n" "$GCAPATH" >> $FOLDERPATH/genomes.txt
                    println "%s\n" "$GCA" >> $FOLDERPATH/${CATEGORY}.txt
                else
                    # Delete corrupted genome
                    rm $OUTDIR/${GCA}.fna.gz
                    printf "[TXID=%s][ID=%s] Corrupted genome\n" "$TXID" "$GCA"
                fi
            else
                # Report missing genomes
                printf "[TXID=%s][ID=%s] Unable to download genome\n" "$TXID" "$GCA"
            fi
        fi
        done
