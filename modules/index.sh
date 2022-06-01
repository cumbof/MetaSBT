#!/bin/bash
#title          :index
#description    :Retrieve genomes from isolate sequencing from the NCBI GenBank and build a Sequence Bloom Tree for each species with kmtricks 
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#=============================================================================================================================================

DATE="May 31, 2022"
VERSION="0.1.0"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Define default value for --nproc and --xargs-nproc
NPROC=1
XARGS_NPROC=1
# Initialize CheckM completeness and contamination to default values
# Run CheckM if provided completeness is >0 or contamination is <100
CHECKM_COMPLETENESS=0
CHECKM_CONTAMINATION=100
# Dereplicate input genomes
# Use kmtricks to build the kmer matrix on the input set of genomes
# Remove genomes with identical set of kmers
DEREPLICATE=false
# Estimate the bloom filter size
ESTIMATE_FILTER_SIZE=false
# Remove temporary data at the end of the pipeline
CLEANUP=false

# Parse input arguments
for ARG in "$@"; do
    case "$ARG" in
        --checkm-completeness=*)
            # CheckM completeness
            CHECKM_COMPLETENESS="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_COMPLETENESS}" =~ "?" ]]; then
                printf "index helper: --checkm-completeness=num\n\n"
                printf "\tInput genomes must have a minimum completeness percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-completeness is an integer
            if [[ ! ${CHECKM_COMPLETENESS} =~ ^[0-9]+$ ]] || [[ "${CHECKM_COMPLETENESS}" -gt "100" ]]; then
                printf "Argument --checkm-completeness must be a positive integer greater than 0 (up to 100)\n"
                exit 1
            fi
            ;;
        --checkm-contamination=*)
            # CheckM contamination
            CHECKM_CONTAMINATION="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_CONTAMINATION}" =~ "?" ]]; then
                printf "index helper: --checkm-contamination=num\n\n"
                printf "\tInput genomes must have a maximum contamination percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-contamination is an integer
            if [[ ! ${CHECKM_CONTAMINATION} =~ ^[0-9]+$ ]] || [[ "${CHECKM_CONTAMINATION}" -gt "100" ]]; then
                printf "Argument --checkm-completeness must be a positive integer lower than 100 (up to 0)\n"
                exit 1
            fi
            ;;
        --cleanup)
            # Remove temporary data at the end of the pipeline
            CLEANUP=true
            ;;
        --db-dir=*)
            # Database directory
            DBDIR="${ARG#*=}"
            # Define helper
            if [[ "${DBDIR}" =~ "?" ]]; then
                printf "update helper: --db-dir=directory\n\n"
                printf "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "${DBDIR}" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
            # Trim the last slash out of the path
            DBDIR="${DBDIR%/}"
            ;;
        --dereplicate)
            # Dereplicate input genomes
            DEREPLICATE=true
            ;;
        --estimate-filter-size)
            # Automatically estimate the best bloom filter size with ntCard
            ESTIMATE_FILTER_SIZE=true
            ;;
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
            source ${SCRIPT_DIR}/../HELP
            exit 0
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
        --tmp-dir=*)
            # Temporary folder
            TMPDIR="${ARG#*=}"
            # Define helper
            if [[ "${TMPDIR}" =~ "?" ]]; then
                printf "update helper: --tmp-dir=directory\n\n"
                printf "\tPath to the folder for storing temporary data.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            TMPDIR="$( cd "$( dirname "${TMPDIR}" )" &> /dev/null && pwd )"/"$( basename $TMPDIR )"
            # Trim the last slash out of the path
            TMPDIR="${TMPDIR%/}"
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

check_dependencies false
if [[ "$?" -gt "0" ]]; then
    printf "Unsatisfied software dependencies!\n\n"
    printf "Please run the following command for a list of required external software dependencies:\n\n"
    printf "\t$ meta-index --resolve-dependencies\n\n"
    
    exit 1
fi

# Init database directory
mkdir -p $DBDIR
# Init temporary folder
mkdir -p $TMPDIR

# Download taxonomy dump from NCBI
TAXDUMP="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
printf "Downloading taxonomy dump from NCBI\n"
printf "\t%s\n\n" "$TAXDUMP"

wget -N $TAXDUMP -P $TMPDIR
mkdir -p $TMPDIR/taxdump && tar zxf $TMPDIR/taxdump.tar.gz -C $TMPDIR/taxdump

# Export lineages
printf "Exporting lineages\n"
ncbitax2lin --nodes-file $TMPDIR/taxdump/nodes.dmp \
            --names-file $TMPDIR/taxdump/names.dmp \
            --output $TMPDIR/ncbi_lineages.csv.gz

# Build a mapping between tax id and full taxonomy
# Consider a specific kingdom only
# Remove special characters
# Fill empty taxa level with "unclassified"
printf "\nBuilding tax id to full taxonomy mapping\n\n"
zcat $TMPDIR/ncbi_lineages.csv.gz | \
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
                                     }' > $TMPDIR/taxa.tsv

# Download all GCAs associated to the taxa IDs in taxa.tsv
# Use genomes that have not been excluded from RefSeq
# https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
printf "Downloading genomes from NCBI GenBank\n"
if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
    printf "\tRunning CheckM to quality control input genomes\n"
    printf "\t\tMinimum completeness: %s\n" "${CHECKM_COMPLETENESS}"
    printf "\t\tMaximum contamination: %s\n\n" "${CHECKM_CONTAMINATION}"
    printf "\tDereplicating input genomes\n\n"
fi

while read tax_id, taxonomy; do
    # Do not consider "unclassified" genomes since two genomes under 
    # the same "unclassified" taxonomic label could potentially be completely different
    # Unclassified genomes can be used with the "update" module with "--type=MAGs"
    if [[ ! $taxonomy == *"_unclassified"* ]]; then
        printf "\tProcessing %s: " "$taxonomy"
        # Index reference genomes
        SEARCH_CRITERIA="NOT excluded-from-refseq [PROP]"
        esearch_txid $DBDIR ${tax_id} $taxonomy "references" ${SEARCH_CRITERIA}

        # Taxonomy folder, genomes directory, and fof file are generated by the "esearch_txid" function
        TAXDIR=$DBDIR/$(echo "$taxonomy" | sed 's/|/\//g')
        GENOMES_FOF=${TAXDIR}/genomes.fof
        GENOMES_DIR=${TAXDIR}/genomes

        # Input genomes must be quality-controlled before being added to the database
        if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
            CHECKMTABLES="" # This is the result of the "run_checkm" function as a list of file paths separated by comma
            run_checkm ${TAXDIR}/genomes.txt "fna.gz" $TMPDIR $NPROC

            # Process all the CheckM tables
            CHECKMTABLES_ARR=($(echo $CHECKMTABLES | tr "," " "))
            for table in ${CHECKMTABLES_ARR[@]}; do
                # Skip header line with sed
                # Discard genomes according to their completeness and contamination
                sed 1d $table | while read line; do
                    # Retrieve completeness and contamination of current genome
                    GENOMEID="$(echo $line | cut -d$'\t' -f1)"          # Get value under column 1
                    COMPLETENESS="$(echo $line | cut -d$'\t' -f12)"     # Get value under column 12
                    CONTAMINATION="$(echo $line | cut -d$'\t' -f13)"    # Get value under column 13
                    if [[ "$COMPLETENESS" -lt "${CHECKM_COMPLETENESS}" ]] || [[ "$CONTAMINATION" -gt "${CHECKM_CONTAMINATION}" ]]; then
                        # Remove genome from directory
                        rm -rf ${GENOMES_DIR}/${GENOMEID}.fna.gz
                    else
                        # Current genome passed the quality control
                        grep -w "${GENOMEID}.fna.gz" ${GENOMES_FOF} >> ${TAXID}/genomes_qc.fof
                    fi
                done
            done

            # Proceed if at least one input genome passed the quality control
            if [[ -f ${TAXID}/genomes_qc.fof ]]; then
                GENOMES_FOF=${TAXID}/genomes_qc.fof
            else
                # No input genomes survived the quality control step
                # Remove the genomes.fof file to avoind running kmtricks on current taxonomy
                rm ${GENOMES_FOF}
                # Skip the dereplication step and proceed with the next taxonomy
                continue
            fi
        fi

        # Use kmtricks to build a kmer matrix and compare input genomes
        # Discard genomes with the same set of kmers (dereplication)
        if ${DEREPLICATE}; then
            if [[ -f ${GENOMES_FOF} ]]; then
                HOW_MANY=$(cat ${GENOMES_FOF} | wc -l)
                if [[ "${HOW_MANY}" -gt "1" ]]; then
                    printf "\t[TAXID=%s] Dereplicating %s genomes\n" "${tax_id}" "${HOW_MANY}"
                    # Run kmtricks to build the kmers matrix
                    kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                            ${TAXDIR}/matrix \
                                            ${NPROC} \
                                            ${TAXDIR}/kmers_matrix.txt

                    # Add header to the kmers matrix
                    HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
                    echo "#kmer $HEADER" > ${TAXDIR}/kmers_matrix_wh.txt
                    cat ${TAXDIR}/kmers_matrix.txt >> ${TAXDIR}/kmers_matrix_wh.txt
                    mv ${TAXDIR}/kmers_matrix_wh.txt ${TAXDIR}/kmers_matrix.txt
                    # Remove duplicate columns
                    cat ${TAXDIR}/kmers_matrix.txt | transpose | awk '!seen[substr($0, index($0, " "))]++' | transpose > ${TAXDIR}/kmers_matrix_dereplicated.txt

                    # Dereplicate genomes
                    while read -r genome; do
                        if ! head -n1 ${TAXDIR}/kmers_matrix_dereplicated.txt | grep -q -w "$genome"; then
                            # Remove genome from directory
                            rm -rf ${GENOMES_DIR}/${genome}.fna.gz
                        else
                            grep "^${genome} " ${GENOMES_FOF} >> ${TAXDIR}/genomes_dereplicated.fof
                        fi
                    done <<<"$(head -n1 ${TAXDIR}/kmers_matrix.txt | tr " " "\n")"

                    # Cleanup space
                    rm -rf ${TAXDIR}/matrix
                    mv ${TAXDIR}/kmers_matrix_dereplicated.txt ${TAXDIR}/kmers_matrix.txt
                    if [[ -f ${TAXDIR}/genomes_dereplicated.fof ]]; then
                        mv ${TAXDIR}/genomes_dereplicated.fof ${GENOMES_FOF}
                    fi
                fi
            fi
        fi

        # Count how many genomes survived the quality control and dereplication steps
        HOW_MANY=$(cat ${GENOMES_FOF} | wc -l)
        printf "%s genomes\n" "${HOW_MANY}"
        if [[ "${HOW_MANY}" -eq "0" ]]; then
            # Remove the species folder if it does not contain any genome
            rm -rf $TAXDIR
        fi
    fi
done < $TMPDIR/taxa.tsv

# Automatically estimate the best bloom filter size with ntCard
if ${ESTIMATE_FILTER_SIZE}; then
    printf "\nRunning ntCard for estimating the bloom filter size\n"
    # Create a list with the file paths to all the downloaded genomes
    find $DBDIR -type f -iname "genomes.fof" -follow | xargs -n 1 -I {} bash -c \
        'INPUT={}; \
         cat $INPUT | cut -d" " -f3 >> '"$TMPDIR"'/genomes.txt'
    # Call ntCard
    # With --kmer=31 it will produce genomes_k21.hist
    ntcard --kmer=${KMER_LEN} --threads=$NPROC --pref==$TMPDIR/genomes @$TMPDIR/genomes.txt
    if [[ -f $TMPDIR/genomes_k${KMER_LEN}.hist ]]; then
        F0=$(grep "^F0" $TMPDIR/genomes_k${KMER_LEN}.hist | cut -f2)
        f1=$(grep "^1" $TMPDIR/genomes_k${KMER_LEN}.hist | head -n 1 | cut -f2)
        FILTER_SIZE=$(( ${F0}-${f1} ))
    else
        printf "\n[ERROR] An error has occurred while running ntCard\n"
        exit 1
    fi
fi

# Run kmtricks and build a sequence bloom tree for each species
printf "\nRunning kmtricks at the species level\n"
NFOFIN=`find $DBDIR -type f -iname "genomes.fof" | wc -l`   # Count how many genomes.fof files (one for each species)
NFOFOUT=`find $DBDIR -type f -iname "kmtricks.fof" | wc -l` # Count how many kmtricks.fof files (copy of genomes.fof generated by kmtricks)
while [ $NFOFOUT -lt $NFOFIN ]; do
    find $DBDIR -type f -iname "genomes.fof" -follow | xargs -n 1 -P ${XARGS_NPROC} -I {} bash -c \
        'INPUT={}; \
         if [[ ! -f "$(dirname $INPUT)/index/kmtricks.fof" ]]; then \
            kmtricks_index_wrapper "$INPUT" '"${DBDIR}"' '"${KMER_LEN}"' '"${FILTER_SIZE}"' '"${NPROC}"'; \
         fi'
    # Search for missing species and process them one by one
    NFOFOUT=`find $DBDIR -type f -iname "kmtricks.fof" | wc -l`
    XARGS_NPROC=1
done

# Run howdesbt to build a sequence bloom tree for each of the lower taxonomic levels
printf "\nRunning howdesbt on lower taxonomic levels\n"
# Process lower taxonomic levels
# Start with genera at depth 6
DEPTH=6
for LEVEL in "g__" "f__" "o__" "c__" "p__" "k__"; do
    find $DBDIR -maxdepth $DEPTH -type d -iname "${LEVEL}*" -follow -exec howdesbt_wrapper {} ${FILTER_SIZE} \;
    DEPTH=$(expr $DEPTH - 1)
done

# Cleanup temporary data
if ${CLEANUP}; then
    # Remove tmp folder
    if [[ -d $TMPDIR ]]; then
        printf "\nCleaning up temporary folder:\n"
        printf "\t%s\n" "${TMPDIR}"
        rm -rf ${TMPDIR}
    fi
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
printf "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0