#!/bin/bash
#title          :index
#description    :Retrieve reference genomes and MAGs from NCBI GenBank and build a Sequence Bloom Tree for each species with kmtricks 
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#====================================================================================================================================

DATE="May 18, 2022"
VERSION="0.1.0"

# Define script directory
INDEX_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${INDEX_DIR}/utils.sh

# Download genomes from NCBI GenBank
esearch_txid () {
    DB_DIR=$1
    TAX_ID=$2
    FULL_TAXONOMY=$3
    GENOME_CATEGORY=$4
    SEARCH_CRITERIA=$5
    # Download GCAs associated to a specific tax_id
    esearch -db assembly -query "txid${TAX_ID} ${SEARCH_CRITERIA}" < /dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r URL; do
            # Create a directory for the current taxonomy
            TAXDIR=${DB_DIR}/$(echo "${FULL_TAXONOMY}" | sed 's/|/\//g')
            OUTDIR=$TAXDIR/genomes
            mkdir -p $OUTDIR
            # Download GCA
            FNAME=$(echo $URL | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
            GCA=${FNAME%%.*}
            if [[ ! -f "${OUTDIR}/${GCA}.fna.gz" ]]; then
                wget -q "$URL/$FNAME" -O ${OUTDIR}/${GCA}.fna.gz
                if [[ -f "${OUTDIR}/${GCA}.fna.gz" ]]; then
                    if gzip -t ${OUTDIR}/${GCA}.fna.gz; then
                        # Define a file of file (fof) with the list of genomes for current species
                        GCAPATH=$(readlink -m ${OUTDIR}/${GCA}.fna.gz)
                        printf "%s : %s\n" "$GCA" "$GCAPATH" >> $TAXDIR/genomes.fof
                        printf "%s\n" "$GCA" >> $TAXDIR/${GENOME_CATEGORY}.txt
                    else
                        # Delete corrupted genome
                        rm ${OUTDIR}/${GCA}.fna.gz
                        printf "\t[ERROR][ID=%s][TAXID=%s] Corrupted genome\n" "$GCA" "${TAX_ID}"
                    fi
                else
                    # Report missing genomes
                    printf "\t[ERROR][ID=%s][TAXID=%s] Unable to download genome\n" "$GCA" "${TAX_ID}"
                fi
            fi
          done
}

# Wrapper for the howdesbt pipeline
# This is applied on all the taxonomic level except the species
howdesbt_wrapper () {
    LEVEL_DIR=$1    # Taxonomic level folder
    FILTER_SIZE=$2  # Bloom filter size
    LEVEL_NAME=$(basename ${LEVEL_DIR})
    for DIRECTORY in ${LEVEL_DIR}/*/; do
        UPPER_LEVEL=$(basename $DIRECTORY)
        # Define a list of bloom filters
        # Each file is the result of the OR operator on all representative bloom filters in the upper taxonomic level
        echo "$(readlink -m ${LEVEL_DIR}/${UPPER_LEVEL}/${UPPER_LEVEL}.bf)" >> ${LEVEL_DIR}/${LEVEL_NAME}.txt
    done

    mkdir -p ${LEVEL_DIR}/index
    HOWMANY=$(wc -l ${LEVEL_DIR}/${LEVEL_NAME}.txt | cut -d" " -f1)
    if [[ "${HOWMANY}" -gt "1" ]]; then
        # Create a tree topology file with howdesbt
        howdesbt cluster --list=${LEVEL_DIR}/${LEVEL_NAME}.txt \
                         --bits=${FILTER_SIZE} \
                         --tree=${LEVEL_DIR}/index/union.sbt \
                         --nodename=node{number} \
                         --keepallnodes
        # Build the bloom filter files for the tree
        hodesbt build --HowDe \
                      --tree=${LEVEL_DIR}/index/union.sbt \
                      --outtree=${LEVEL_DIR}/index/index.detbrief.sbt
        # Merge all the leaves together by applying the OR logical operator on the bloom filter files
        # The resulting bloom filter is the representative one, which is the same as the root node of the tree
        while read BFPATH; do
            if [[ ! -f ${LEVEL_DIR}/${LEVEL_NAME}.bf ]]; then
                cp ${BFPATH} ${LEVEL_DIR}/${LEVEL_NAME}.bf
            else
                howdesbt bfoperate ${BFPATH} ${LEVEL_DIR}/${LEVEL_NAME}.bf --or --out ${LEVEL_DIR}/merged.bf
                mv ${LEVEL_DIR}/merged.bf ${LEVEL_DIR}/${LEVEL_NAME}.bf
            fi
        done < ${LEVEL_DIR}/${LEVEL_NAME}.txt
    else
        # With only one bloom filter, it does not make sense to build an index with howdesbt
        BFPATH=$(head -n1 "${LEVEL_DIR}/${LEVEL_NAME}.txt")
        cp ${BFPATH} ${LEVEL_DIR}/${LEVEL_NAME}.bf
    fi
}
# Export howdesbt_wrapper to sub-shells
export -f howdesbt_wrapper

# Wrapper for the kmtricks pipeline
# This is applied at the species level only
kmtricks_wrapper () {
    INPUT=$1        # Path to the genomes.fof file for a specific species
    DBDIR=$2        # Database root directory path
    KMER_LEN=$3     # Length of the kmers
    FILTER_SIZE=$4  # Bloom filter size
    NPROC=$5        # Max nproc for multiprocessing
    FOLDERPATH=$(dirname "${INPUT}")
    if [[ ! -f "${FOLDERPATH}/index/kmtricks.fof" ]]; then
        # Take track of the processed species in log
        echo ${INPUT} >> ${DBDIR}/kmtricks.log

        # Run the kmtricks pipeline
        kmtricks pipeline --file ${INPUT} \
                          --run-dir ${FOLDERPATH}/index \
                          --kmer-size ${KMER_LEN} \
                          --mode hash:bft:bin \
                          --hard-min 1 \
                          --bloom-size ${FILTER_SIZE} \
                          --bf-format howdesbt \
                          --cpr \
                          --skip-merge \
                          -t ${NPROC}
        
        FOLDERNAME="${FOLDERPATH##*/}"
        HOWMANY=$(wc -l ${INPUT} | cut -d" " -f1)
        if [[ "${HOWMANY}" -gt "1" ]]; then
            # Index genomes by building a sequence bloom tree with howdesbt
            kmtricks index --run-dir ${FOLDERPATH}/index \
                           --howde \
                           --threads ${NPROC}
            # Merge all the leaves together by applying the OR logical operator on the bloom filter files
            # The resulting bloom filter is the representative one, which is the same as the root node of the tree
            while read line; do
                GENOME=$(echo "${line}" | cut -d" " -f1)
                if [[ ! -f ${FOLDERPATH}/${FOLDERNAME}.bf ]]; then
                    cp ${FOLDERPATH}/index/filters/${GENOME}.bf ${FOLDERPATH}/${FOLDERNAME}.bf
                else
                    howdesbt bfoperate ${FOLDERPATH}/index/filters/${GENOME}.bf ${FOLDERPATH}/${FOLDERNAME}.bf --or --out ${FOLDERPATH}/merged.bf
                    mv ${FOLDERPATH}/merged.bf ${FOLDERPATH}/${FOLDERNAME}.bf
                fi
            done < ${INPUT}
        else
            # With only one genome, it does not make sense to build an index with howdesbt
            mkdir -p ${FOLDERPATH}/index/howde_index
            cp ${INPUT} ${FOLDERPATH}/index/kmtricks.fof
            GENOME=$(head -n1 "${INPUT}" | cut -d" " -f1)
            cp ${FOLDERPATH}/index/filters/${GENOME}.bf ${FOLDERPATH}/${FOLDERNAME}.bf
        fi
    fi
}
# Export kmtricks_wrapper to sub-shells
export -f kmtricks_wrapper

# Define default value for --nproc and --xargs-nproc
NPROC=1
XARGS_NPROC=1
# Download and index reference genomes from isolate sequencing by default
# Disable with --no-references
INDEX_REFERENCES=true
# Do not download MAGs by default
# Enable with --index-mags
INDEX_MAGS=false

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --filter-size=*)
            # Bloom filter size
            FILTER_SIZE="${ARG#*=}"
            # Define helper
            if [[ "${FILTER_SIZE}" =~ "?" ]]; then
                printf "index helper: --filter-size=num\n\n"
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
            INDEX_HELP=true
            source ${INDEX_DIR}/../HELP
            exit 0
            ;;
        --index-mags)
            # Download and index MAGs
            INDEX_MAGs=true
            ;;
        --kingdom=*)
            # Consider genomes whose lineage belong to a specific kingdom only
            KINGDOM="${ARG#*=}"
            # Define helper
            if [[ "${KINGDOM}" =~ "?" ]]; then
                printf "index helper: --kingdom=value\n\n"
                printf "\tConsider genomes whose lineage belongs to a specific kingdom only.\n\n"
                exit 0
            fi
            ;;
        --kmer-len=*)
            # Length of the kmers
            KMER_LEN="${ARG#*=}"
            # Define helper
            if [[ "${KMER_LEN}" =~ "?" ]]; then
                printf "index helper: --kmer-len=num\n\n"
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
            printf "%s\n" "$(cat ${INDEX_DIR}/../LICENSE)"
            exit 0
            ;;
        --no-references)
            # Do not download and index reference genomes from isolate sequencing
            INDEX_REFERENCES=false
            ;;
        --nproc=*)
            # Max nproc for all parallel instructions
            NPROC="${ARG#*=}"
            # Define helper
            if [[ "${NPROC}" =~ "?" ]]; then
                printf "index helper: --nproc=num\n\n"
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
        -v|--version)
            # Print pipeline version
            printf "index version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        --work-dir=*)
            # Working directory
            WORKDIR="${ARG#*=}"
            # Define helper
            if [[ "${WORKDIR}" =~ "?" ]]; then
                printf "index helper: --work-dir=directory\n\n"
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
                printf "index helper: --xargs-nproc=num\n\n"
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
            printf "index: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

printf "index version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies
if [[ "$?" -gt "0" ]]; then
    exit 1
fi

# Init database directory
DBDIR=$WORKDIR/db

# Download taxonomy dump from NCBI
TAXDUMP="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
printf "Downloading taxonomy dump from NCBI\n"
printf "\t%s\n\n" "$TAXDUMP"

wget -N $TAXDUMP -P $WORKDIR
mkdir -p $WORKDIR/taxdump && tar zxf $WORKDIR/taxdump.tar.gz -C $WORKDIR/taxdump

# Export lineages
printf "Exporting lineages\n"
ncbitax2lin --nodes-file $WORKDIR/taxdump/nodes.dmp \
            --names-file $WORKDIR/taxdump/names.dmp \
            --output $WORKDIR/ncbi_lineages.csv.gz

# Build a mapping between tax id and full taxonomy
# Consider a specific kingdom only
# Remove special characters
# Fill empty taxa level with "unclassified"
printf "Building tax id to full taxonomy mapping\n"
zcat $WORKDIR/ncbi_lineages.csv.gz | \
    awk -v kingdom="$KINGDOM" 'BEGIN { FS=","; OFS="" } 
                                     { gsub(" ", "_"); gsub(/\.|\047|\"|\(|\)|\:/, "") }
                                     {
                                        if ($3=="") {
                                            phylum_pre=$2;
                                            phylum_suf="_unclassified";
                                        } else {
                                            phylum_pre=$3;
                                            phylum_suf="";
                                        }
                                     } {
                                        if ($4=="") {
                                            class_pre=phylum_pre; 
                                            class_suf="_unclassified";
                                        } else { 
                                            class_pre=$4; 
                                            class_suf="";
                                        }
                                     } {
                                        if ($5=="") { 
                                            order_pre=class_pre; 
                                            order_suf="_unclassified";
                                        } else { 
                                            order_pre=$5; 
                                            order_suf="";
                                        }
                                     } {
                                        if ($6=="") {
                                            family_pre=order_pre; 
                                            family_suf="_unclassified";
                                        } else { 
                                            family_pre=$6; 
                                            family_suf="";
                                        }
                                     } {
                                        if ($7=="") { 
                                            genus_pre=family_pre; 
                                            genus_suf="_unclassified"; 
                                        } else { 
                                            genus_pre=$7; 
                                            genus_suf="";
                                        }
                                     } {
                                        if ($8=="") { 
                                            species_pre=genus_pre; 
                                            species_suf="_unclassified";
                                        } else {
                                            species_pre=$8; 
                                            species_suf="";
                                        }
                                     } { 
                                        if (NR>1 && $2==kingdom) {
                                            print $1, "\t", "k__", $2,
                                                            "|p__", phylum_pre, phylum_suf,
                                                            "|c__", class_pre, class_suf,
                                                            "|o__", order_pre, order_suf,
                                                            "|f__", family_pre, family_suf,
                                                            "|g__", genus_pre, genus_suf,
                                                            "|s__", species_pre, species_suf;
                                        }
                                     }' > $WORKDIR/taxa.tsv

# Download all GCAs associated to the taxa IDs in taxa.tsv
# Use genomes that have not been excluded from RefSeq
# https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
printf "Downloading genomes from NCBI GenBank\n"
while read tax_id, taxonomy; do
    # Do not consider "unclassified" genomes since two genomes under 
    # the same "unclassified" taxonomic label could potentially be completely different
    # Unclassified genomes can be used with the "update" module with "--type=MAGs"
    if [[ ! $taxonomy == *"_unclassified"* ]]; then
        # Index reference genomes
        if ${INDEX_REFERENCES}; then
            SEARCH_CRITERIA="NOT excluded-from-refseq [PROP]"
            esearch_txid $DBDIR ${tax_id} $taxonomy "references" ${SEARCH_CRITERIA}
        fi
        # Index MAGs
        if ${INDEX_MAGs}; then
            # Run esearch_txid again for downloading genomes with no search criteria
            # Genomes already processed with the previous esearch_txid run are skipped
            esearch_txid $DBDIR ${tax_id} $taxonomy "mags"
        fi
    fi
done < $WORKDIR/taxa.tsv

# Use kmtricks to build a kmer matrix and compare input genomes
# Discard genomes with the same kmers (dereplication)
# TODO

# Run kmtricks and build a sequence bloom tree for each species
printf "Running kmtricks at the species level\n"
NFOFIN=`find $DBDIR -type f -iname "genomes.fof" | wc -l`   # Count how many genomes.fof files (one for each species)
NFOFOUT=`find $DBDIR -type f -iname "kmtricks.fof" | wc -l` # Count how many kmtricks.fof files (copy of genomes.fof generated by kmtricks)
while [ $NFOFOUT -lt $NFOFIN ]; do
    find $DBDIR -type f -iname "genomes.fof" -follow | xargs -n 1 -P ${XARGS_NPROC} -I {} bash -c \
        'kmtricks_wrapper {} '"${DBDIR}"' '"${KMER_LEN}"' '"${FILTER_SIZE}"' '"${NPROC}"
    # Search for missing species and process them one by one
    NFOFOUT=`find $DBDIR -type f -iname "kmtricks.fof" | wc -l`
    XARGS_NPROC=1
done

# Process lower taxonomic levels
LEVELS=(
    "g__" # Genus
    "f__" # Family
    "o__" # Order
    "c__" # Class
    "p__" # Phylum
    "k__" # Kingdom
)
# Start with genera at depth 6
DEPTH=6
for LEVEL in "${LEVELS[@]}"; do
    find $DBDIR -maxdepth $DEPTH -type d -iname "${LEVEL}*" -follow -exec howdesbt_wrapper {} ${FILTER_SIZE} \;
    DEPTH=$(expr $DEPTH - 1)
done

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

exit 0