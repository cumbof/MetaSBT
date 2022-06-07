#!/bin/bash
#title          :index
#description    :Retrieve genomes from isolate sequencing from the NCBI GenBank and build a Sequence Bloom Tree for each species with kmtricks 
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#=============================================================================================================================================

DATE="Jun 7, 2022"
VERSION="0.1.0"

# Take track of the directory from which this module is launched
# This is required in order to come back on this folder after running HowDeSBT
EXEC_DIR="$PWD"

# Define script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Requirements base path
REQUIREMENTS_DIR="$(dirname "${SCRIPT_DIR}")/requirements"
# Import utility functions
source ${SCRIPT_DIR}/utils.sh

# Define default value for --nproc and --xargs-nproc
NPROC=1
XARGS_NPROC=1
# Limit the number of genomes per species
HOW_MANY_GENOMES_PER_SPECIES=-1
# Increase the estimated bloom filter size by this percentage
INCREASE_FILTER_SIZE=0
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
                println "index helper: --checkm-completeness=num\n\n"
                println "\tInput genomes must have a minimum completeness percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-completeness is an integer
            if [[ ! ${CHECKM_COMPLETENESS} =~ ^[0-9]+$ ]] || [[ "${CHECKM_COMPLETENESS}" -gt "100" ]]; then
                println "Argument --checkm-completeness must be a positive integer greater than 0 (up to 100)\n"
                exit 1
            fi
            ;;
        --checkm-contamination=*)
            # CheckM contamination
            CHECKM_CONTAMINATION="${ARG#*=}"
            # Define helper
            if [[ "${CHECKM_CONTAMINATION}" =~ "?" ]]; then
                println "index helper: --checkm-contamination=num\n\n"
                println "\tInput genomes must have a maximum contamination percentage before being processed and added to the database\n\n"
                exit 0
            fi
            # Check whether --checkm-contamination is an integer
            if [[ ! ${CHECKM_CONTAMINATION} =~ ^[0-9]+$ ]] || [[ "${CHECKM_CONTAMINATION}" -gt "100" ]]; then
                println "Argument --checkm-completeness must be a positive integer lower than 100 (up to 0)\n"
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
            if [[ "$DBDIR" =~ "?" ]]; then
                println "update helper: --db-dir=directory\n\n"
                println "\tThis is the database directory with the taxonomically organised sequence bloom trees.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            DBDIR="$( cd "$( dirname "$DBDIR" )" &> /dev/null && pwd )"/"$( basename $DBDIR )"
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
                println "index helper: --filter-size=num\n\n"
                println "\tThis is the size of the bloom filters.\n\n"
                exit 0
            fi
            # Check whether --filter-size is an integer
            if [[ ! ${FILTER_SIZE} =~ ^[0-9]+$ ]] || [[ "${FILTER_SIZE}" -eq "0" ]]; then
                println "Argument --filter-size must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        -h|--help)
            # Print extended help
            INDEX_HELP=true
            source ${SCRIPT_DIR}/../HELP
            exit 0
            ;;
        --how-many=*)
            # Limit the number of genomes per species
            HOW_MANY_GENOMES_PER_SPECIES="${ARG#*=}"
            # Define helper
            if [[ "${HOW_MANY_GENOMES_PER_SPECIES}" =~ "?" ]]; then
                println "index helper: --how-many=num\n\n"
                println "\tLimit the number of genomes per species.\n"
                println "\tThe number of genomes per species is not limited by default.\n\n"
                exit 0
            fi
            # Check whether --how-many is an integer
            if [[ ! ${HOW_MANY_GENOMES_PER_SPECIES} =~ ^[0-9]+$ ]] || [[ "${HOW_MANY_GENOMES_PER_SPECIES}" -eq "0" ]]; then
                println "Argument --how-many must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --increase-filter-size=*)
            # Increase the estimated filter size by the specified percentage
            INCREASE_FILTER_SIZE="${ARG#*=}"
            # Define helper
            if [[ "${INCREASE_FILTER_SIZE}" =~ "?" ]]; then
                println "index helper: --increase-filter-size=num\n\n"
                println "\tIncrease the estimated filter size by the specified percentage.\n"
                println "\tThis is used in conjunction with the --estimate-filter-size argument only.\n"
                println "\tIt is highly recommended to increase the filter size by a good percentage in case you are planning to update the index with new genomes.\n"
                println "\tDefault: --increase-filter-size=0\n\n"
                exit 0
            fi
            # Check whether --increase-filter-size is an integer
            if [[ ! ${INCREASE_FILTER_SIZE} =~ ^[0-9]+$ ]] || [[ "${INCREASE_FILTER_SIZE}" -gt "100" ]]; then
                println "Argument --increase-filter-size must be a positive integer greater than 0 (up to 100)\n"
                exit 1
            fi
            ;;
        --kingdom=*)
            # Consider genomes whose lineage belong to a specific kingdom only
            KINGDOM="${ARG#*=}"
            # Define helper
            if [[ "$KINGDOM" =~ "?" ]]; then
                println "index helper: --kingdom=value\n\n"
                println "\tConsider genomes whose lineage belongs to a specific kingdom only.\n\n"
                exit 0
            fi
            ;;
        --kmer-len=*)
            # Length of the kmers
            KMER_LEN="${ARG#*=}"
            # Define helper
            if [[ "${KMER_LEN}" =~ "?" ]]; then
                println "index helper: --kmer-len=num\n\n"
                println "\tThis is the length of the kmers used for building bloom filters.\n\n"
                exit 0
            fi
            # Check whether --kmer-len is an integer
            if [[ ! ${KMER_LEN} =~ ^[0-9]+$ ]] || [[ "${KMER_LEN}" -eq "0" ]]; then
                println "Argument --kmer-len must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --log=*)
            # Path to the log file
            LOG_FILEPATH="${ARG#*=}"
            # Define helper
            if [[ "${LOG_FILEPATH}" =~ "?" ]]; then
                println "index helper: --log=file\n\n"
                println "\tPath to the log file.\n\n"
                exit 0
            fi
            # Remove the log file if it already exists
            if [[ -f ${LOG_FILEPATH} ]]; then
                rm ${LOG_FILEPATH}
            fi
            ;;
        --nproc=*)
            # Max nproc for all parallel instructions
            NPROC="${ARG#*=}"
            # Define helper
            if [[ "$NPROC" =~ "?" ]]; then
                println "index helper: --nproc=num\n\n"
                println "\tThis argument refers to the number of processors used for parallelizing the pipeline when possible.\n"
                println "\tDefault: --nproc=1\n\n"
                exit 0
            fi
            # Check whether --proc is an integer
            if [[ ! $NPROC =~ ^[0-9]+$ ]] || [[ "$NPROC" -eq "0" ]]; then
                println "Argument --nproc must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        --resolve-dependencies)
            # Check for external software dependencies and python modules
            check_dependencies true "${REQUIREMENTS_DIR}/index.txt"
            exit $?
            ;;
        --tmp-dir=*)
            # Temporary folder
            TMPDIR="${ARG#*=}"
            # Define helper
            if [[ "$TMPDIR" =~ "?" ]]; then
                println "update helper: --tmp-dir=directory\n\n"
                println "\tPath to the folder for storing temporary data.\n\n"
                exit 0
            fi
            # Reconstruct the full path
            TMPDIR="$( cd "$( dirname "$TMPDIR" )" &> /dev/null && pwd )"/"$( basename $TMPDIR )"
            # Trim the last slash out of the path
            TMPDIR="${TMPDIR%/}"
            ;;
        -v|--version)
            # Print pipeline version
            println "index version %s (%s)\n" "$VERSION" "$DATE"
            exit 0
            ;;
        --xargs-nproc=*)
            # Max number of independent runs of kmtricks and howdesbt
            XARGS_NPROC="${ARG#*=}"
            # Define helper
            if [[ "${XARGS_NPROC}" =~ "?" ]]; then
                println "index helper: --xargs-nproc=num\n\n"
                println "\tThis refers to the number of xargs processes used for launching independent runs of kmtricks and howdesbt.\n"
                println "\tDefault: --xargs-nproc=1\n\n"
                exit 0
            fi
            # Check whether --xargs-nproc is an integer
            if [[ ! $XARGS_NPROC =~ ^[0-9]+$ ]] || [[ "$XARGS_NPROC" -eq "0" ]]; then
                println "Argument --xargs-nproc must be a positive integer greater than 0\n"
                exit 1
            fi
            ;;
        *)
            println "index: invalid option -- %s\n" "$ARG"
            exit 1
            ;;
    esac
done

println "index version %s (%s)\n\n" "$VERSION" "$DATE"
PIPELINE_START_TIME="$(date +%s.%3N)"

check_dependencies false "${REQUIREMENTS_DIR}/index.txt"
if [[ "$?" -gt "0" ]]; then
    println "Unsatisfied software dependencies!\n\n"
    println "Please run the following command for a list of required external software dependencies:\n\n"
    println "\t$ meta-index --resolve-dependencies\n\n"
    
    exit 1
fi

# Init database directory
mkdir -p $DBDIR
# Init temporary folder
mkdir -p $TMPDIR

# Download taxonomy dump from NCBI
TAXDUMP="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
println "Downloading taxonomy dump from NCBI\n"
println "\t%s\n\n" "$TAXDUMP"

wget -N $TAXDUMP -P $TMPDIR
mkdir -p $TMPDIR/taxdump
# Extract NCBI taxdump for 
tar zxf $TMPDIR/taxdump.tar.gz -C $TMPDIR/taxdump

# Export lineages
println "Exporting lineages\n"
ncbitax2lin --nodes-file $TMPDIR/taxdump/nodes.dmp \
            --names-file $TMPDIR/taxdump/names.dmp \
            --output $TMPDIR/ncbi_lineages.csv.gz

# Build a mapping between tax id and full taxonomy
# Consider a specific kingdom only
# Remove special characters
# Fill empty taxa level with "unclassified"
println "\nBuilding tax id to full taxonomy mapping\n\n"
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

# Remove duplicate taxonomies
# Subspecies have the same taxonomic label but different NCBI taxa ID
# Keep the first occurrence of each taxonomic label
# They are all already sorted in ascending order on the NCBI taxa ID
awk '!seen[$2]++' $TMPDIR/taxa.tsv > $TMPDIR/taxa_unique.tsv

# Download all GCAs associated to the taxa IDs in taxa.tsv
# Use genomes that have not been excluded from RefSeq
# https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
println "Downloading genomes from NCBI GenBank\n"
if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
    println "\tRunning CheckM to quality control input genomes\n"
    println "\t\tMinimum completeness: %s\n" "${CHECKM_COMPLETENESS}"
    println "\t\tMaximum contamination: %s\n\n" "${CHECKM_CONTAMINATION}"
    println "\tDereplicating input genomes\n\n"
fi

while read tax_id, taxonomy; do
    # Do not consider "unclassified" genomes since two genomes under 
    # the same "unclassified" taxonomic label could potentially be completely different
    # Unclassified genomes can be used with the "update" module with "--type=MAGs"
    if [[ ! $taxonomy == *"_unclassified"* ]]; then
        println "\tProcessing %s: " "$taxonomy"
        # Index reference genomes
        esearch_txid "$DBDIR" ${HOW_MANY_GENOMES_PER_SPECIES} "${tax_id}" "$taxonomy" "references"

        # Taxonomy folder is generated by the "esearch_txid" function
        TAXDIR=$DBDIR/$(echo "$taxonomy" | sed 's/|/\//g')
        if [[ -d $TAXDIR ]]; then
            # Genomes directory and fof file are generated by the "esearch_txid" function
            GENOMES_FOF=$TAXDIR/genomes.fof
            GENOMES_DIR=$TAXDIR/genomes

            # Input genomes must be quality-controlled before being added to the database
            if [[ "${CHECKM_COMPLETENESS}" -gt "0" ]] && [[ "${CHECKM_CONTAMINATION}" -lt "100" ]]; then
                CHECKMTABLES="" # This is the result of the "run_checkm" function as a list of file paths separated by comma
                run_checkm $TAXDIR/genomes.txt "txid${tax_id}" "fna.gz" $TMPDIR $NPROC

                # Process all the CheckM tables
                CHECKMTABLES_ARR=($(echo $CHECKMTABLES | tr "," " "))
                # In case there is at least one CheckM output table
                if [[ "${#CHECKMTABLES_ARR[@]}" -gt "0" ]]; then
                    # Merge all tables
                    # Retrieve header from the first file in list
                    echo "$(head -n1 ${CHECKMTABLES_ARR[0]})" > $TAXDIR/checkm.tsv
                    for table in ${CHECKMTABLES_ARR[@]}; do
                        # Current table may not exist in case CheckM stopped working
                        if [[ -f $table ]]; then
                            # Get all the lines except the first one
                            echo "$(tail -n +2 $table)" >> $TAXDIR/checkm.tsv
                        fi
                    done

                    # Check whether the checkm.tsv table has been created
                    if [[ -f $TAXDIR/checkm.tsv ]]; then
                        # Iterate over input genome IDs
                        for GENOMEID in $(cut -d" " -f1 ${GENOMES_FOF}); do
                            # Search for current genome ID into the CheckM output table
                            GENOME_DATA="$(grep "$GENOMEID"$'\t' $TAXDIR/checkm.tsv)"
                            if [[ ! -z "${GENOME_DATA}" ]]; then
                                # In case the current genome has been processed
                                # Retrieve completeness and contamination of current genome
                                COMPLETENESS="$(echo ${GENOME_DATA} | rev | cut -d' ' -f3 | rev)"  # Get completeness
                                CONTAMINATION="$(echo ${GENOME_DATA} | rev | cut -d' ' -f3 | rev)" # Get contamination
                                # Use bc command for comparing floating point numbers
                                COMPLETENESS_CHECK="$(echo "$COMPLETENESS >= ${CHECKM_COMPLETENESS}" | bc)"
                                CONTAMINATION_CHECK="$(echo "$CONTAMINATION <= ${CHECKM_CONTAMINATION}" | bc)"
                                if [[ "${COMPLETENESS_CHECK}" -eq "1" ]] && [[ "${CONTAMINATION_CHECK}" -eq "1" ]]; then
                                    # Current genome passed the quality control
                                    grep -w "${GENOMEID}.fna.gz" ${GENOMES_FOF} >> $TAXDIR/genomes_qc.fof
                                else
                                    # Remove genome from directory
                                    rm -rf ${GENOMES_DIR}/${GENOMEID}.fna.gz
                                fi
                            else
                                # CheckM failed in processing the current genome
                                # Remove genome from directory
                                rm -rf ${GENOMES_DIR}/${GENOMEID}.fna.gz
                            fi
                        done
                    fi
                fi

                # Proceed if at least one input genome passed the quality control
                if [[ -f $TAXDIR/genomes_qc.fof ]]; then
                    mv $TAXDIR/genomes_qc.fof ${GENOMES_FOF}
                else
                    # No input genomes survived the quality control step
                    # Remove the species folder to avoid running kmtricks on current taxonomy
                    rm -rf $TAXDIR
                    # Skip the dereplication step and proceed with the next taxonomy
                    continue
                fi
            fi

            # Use kmtricks to build a kmer matrix and compare input genomes
            # Discard genomes with the same set of kmers (dereplication)
            if $DEREPLICATE; then
                if [[ -f ${GENOMES_FOF} ]]; then
                    HOW_MANY=$(cat ${GENOMES_FOF} | wc -l)
                    if [[ "${HOW_MANY}" -gt "1" ]]; then
                        println "\t[TAXID=%s] Dereplicating %s genomes\n" "${tax_id}" "${HOW_MANY}"
                        # Run kmtricks to build the kmers matrix
                        kmtricks_matrix_wrapper ${GENOMES_FOF} \
                                                $TAXDIR/matrix \
                                                $NPROC \
                                                $TAXDIR/kmers_matrix.txt \
                                                $TAXDIR

                        # Add header to the kmers matrix
                        HEADER=$(awk 'BEGIN {ORS = " "} {print $1}' ${GENOMES_FOF})
                        if [[ "${HEADER: -1}" = " " ]]; then
                            # Trim the last character out of the header in xase of space
                            HEADER="${HEADER:0:-1}"
                        fi
                        echo "#kmer $HEADER" > $TAXDIR/kmers_matrix_wh.txt
                        cat $TAXDIR/kmers_matrix.txt >> $TAXDIR/kmers_matrix_wh.txt
                        mv $TAXDIR/kmers_matrix_wh.txt $TAXDIR/kmers_matrix.txt

                        # Dereplicate genomes
                        python3 ${SCRIPT_DIR}/utils.py filter_genomes $TAXDIR/kmers_matrix.txt $TAXDIR/filtered.txt > /dev/null

                        if [[ -f $TAXDIR/filtered.txt ]]; then
                            # Dereplicate genomes
                            FILTERED_ARRAY=()
                            while read -r genome; do
                                if grep -q "^$genome$" $TAXDIR/filtered.txt; then
                                    # Remove genome from directory
                                    rm -rf ${GENOMES_DIR}/${genome}.fna.gz
                                    FILTERED_ARRAY+=("$genome")
                                else
                                    # Build the new fof file with the dereplicate genomes
                                    grep "^$genome " ${GENOMES_FOF} >> $TAXDIR/genomes_dereplicated.fof
                                fi
                            done <<<"$(head -n1 $TAXDIR/kmers_matrix.txt | tr " " "\n")"

                            if [[ -f $TAXDIR/genomes_dereplicated.fof ]]; then
                                # Build the new kmers matrix with the dereplicated genomes
                                FILTERED="$(printf ",%s" "${FILTERED_ARRAY[@]}")"
                                FILTERED="${FILTERED:1}"
                                csvcut --delimiter " " --not-columns $FILTERED $TAXDIR/kmers_matrix.txt > $TAXDIR/kmers_matrix_dereplicated.txt
                                # Point GENOMES_FOF to the new fof with the dereplicated genomes
                                mv $TAXDIR/genomes_dereplicated.fof ${GENOMES_FOF}
                            fi
                        fi

                        # Cleanup space
                        rm -rf $TAXDIR/matrix
                    fi
                fi
            fi

            # Count how many genomes survived the quality control and dereplication steps
            HOW_MANY=$(cat ${GENOMES_FOF} | wc -l)
            println "%s genomes\n" "${HOW_MANY}"
            if [[ "${HOW_MANY}" -eq "0" ]]; then
                # Remove the species folder if it does not contain any genome
                rm -rf $TAXDIR
            fi
        fi
    fi
done < $TMPDIR/taxa_unique.tsv

# Automatically estimate the best bloom filter size with ntCard
if ${ESTIMATE_FILTER_SIZE}; then
    println "\nRunning ntCard for estimating the bloom filter size\n"
    # Create a list with the file paths to all the downloaded genomes
    find $DBDIR -type f -iname "genomes.fof" -follow | xargs -n 1 -I {} bash -c \
        'INPUT={}
         cat $INPUT | cut -d" " -f3 >> '"$TMPDIR"'/genomes.txt'
    # Call ntCard
    # With --kmer=31 it will produce genomes_k21.hist
    ntcard --kmer=${KMER_LEN} --threads=$NPROC --pref=$TMPDIR/genomes @$TMPDIR/genomes.txt
    if [[ -f $TMPDIR/genomes_k${KMER_LEN}.hist ]]; then
        F0=$(grep "^F0"$'\t' $TMPDIR/genomes_k${KMER_LEN}.hist | cut -f2)
        f1=$(grep "^1"$'\t' $TMPDIR/genomes_k${KMER_LEN}.hist | cut -f2)
        FILTER_SIZE=$(( ${F0}-${f1} ))
        # Compute the increment
        # The result is an integer
        INCREMENT=$(( ${FILTER_SIZE}*${INCREASE_FILTER_SIZE}/100 ))
        # Increment the estimated bloom filter size
        FILTER_SIZE=$(( ${FILTER_SIZE}+${INCREMENT} ))
        println "\tEstimated bloom filter size: %s\n" "${FILTER_SIZE}"
    else
        println "\n[ERROR] An error has occurred while running ntCard\n\n"
        # Print the standard error message and exit
        standard_error_message
        exit 1
    fi
fi

# Take track of the kmer length and bloom filter size in the manifest
printf "%s\n" "--filter-size=${FILTER_SIZE}" \
              "--kingdom=$KINGDOM" \
              "--kmer-len=${KMER_LEN}" \
              > $DBDIR/k__$KINGDOM/manifest.txt

# Start with species at depth 7 up to the kingdom level
DEPTH=7
for LEVELNAME in "species" "genus" "family" "order" "class" "phylum" "kingdom"; do
    # Retrieve the level ID from the level name
    LEVEL="${LEVELNAME:0:1}__"

    # Iterate over the 7 taxonomic levels
    if [[ "$LEVELNAME" = "species" ]]; then
        println "\nRunning kmtricks at the species level\n"
        # Process species with kmtricks
        NFOFIN=$(find $DBDIR -type f -iname "genomes.fof" | wc -l)   # Count how many genomes.fof files (one for each species)
        NFOFOUT=$(find $DBDIR -type f -iname "kmtricks.fof" | wc -l) # Count how many kmtricks.fof files (copy of genomes.fof generated by kmtricks)
        # Set the number of parallel kmtricks runs
        MULTI_RUNS=${XARGS_NPROC}
        while [[ "$NFOFOUT" -lt "$NFOFIN" ]]; do
            find $DBDIR -type f -iname "genomes.fof" -follow | xargs -n 1 -P ${MULTI_RUNS} -I {} bash -c \
                'INPUT={}
                 if [[ ! -f "$(dirname $INPUT)/index/kmtricks.fof" ]]; then
                    kmtricks_index_wrapper "$INPUT" '"$DBDIR"' '"${KMER_LEN}"' '"${FILTER_SIZE}"' '"$NPROC"'
                 fi'
            # Search for missing species and process them one by one
            NFOFOUT=`find $DBDIR -type f -iname "kmtricks.fof" | wc -l`
            MULTI_RUNS=1
        done
    else
        println "Running kmtricks at the %s level\n" "$LEVELNAME"
        # Process all the other taxonomic levels with howdesbt
        NLEVELDIR=$(find $DBDIR -type d -iname "$LEVEL*" | wc -l)                       # Count the number of level folders
        NHOWDELOG=$(find $DBDIR -maxdepth $DEPTH -type f -iname "howdesbt.log" | wc -l) # Count how many times howdesbt has been lauched
        # Set the number of parallel howdesbt runs
        MULTI_RUNS=${XARGS_NPROC}
        while [[ "$NHOWDELOG" -lt "$NLEVELDIR" ]]; then
            find $DBDIR -maxdepth $DEPTH -type d -iname "${LEVEL}*" -follow | xargs -n 1 -P ${MULTI_RUNS} -I {} bash -c \
                'LEVELDIR={}
                 howdesbt_wrapper $LEVELDIR '"${FILTER_SIZE}"';'
            # Search for missing level folders and process them one by one
            NHOWDELOG=$(find $DBDIR -maxdepth $DEPTH -type f -iname "howdesbt.log" | wc -l)
            MULTI_RUNS=1
    fi
    # Decrease the depth while moving to a higher taxonomic level
    DEPTH=$(expr $DEPTH - 1)
done

# HowDeSBT calls automatically change the current folder to the taxonomy directory
# Come back to the folder from which this module has been launched
cd ${EXEC_DIR}

# Cleanup temporary data
if $CLEANUP; then
    # Remove tmp folder
    if [[ -d $TMPDIR ]]; then
        println "\nCleaning up temporary folder:\n"
        println "\t%s\n" "$TMPDIR"
        rm -rf $TMPDIR
    fi
fi

PIPELINE_END_TIME="$(date +%s.%3N)"
PIPELINE_ELAPSED="$(bc <<< "${PIPELINE_END_TIME}-${PIPELINE_START_TIME}")"
println "\nTotal elapsed time: %s\n\n" "$(displaytime ${PIPELINE_ELAPSED})"

# Print credits
credits

exit 0