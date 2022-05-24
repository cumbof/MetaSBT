#!/bin/bash
#title          :utils
#description    :Some utility functions
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="May 24, 2022"
VERSION="0.1.0"

# Check for external software dependencies
check_dependencies () {
    printf "Checking for software dependencies\n"
    # Define the set of dependencies
    DEPENDENCIES=("bc" "checkm" "howdesbt" "kmtricks" "ncbitax2lin" "wget")

    # Count how many missing dependencies
    MISSING=0
    for dep in "${DEPENDENCIES[@]}"; do
        # Check for dependency
        if ! command -v $dep &> /dev/null ; then
            printf "\t[--] %s\n" "$dep"
            MISSING=$(($MISSING + 1))
        else
            printf "\t[OK] %s\n" "$dep"
        fi
    done
    
    if [ "$MISSING" -gt "0" ]; then
        printf "\nPlease, install all the missing dependencies and try again.\n\n"
        return 1
    fi
    
    printf "\nAll required dependencies satisfied!\n\n"
    return 0
}

# Check for new software release
check_for_software_updates () {
    LAST_SOFTWARE_VERSION="$(get_last_software_release)"
    VERSION=$1
    if [ ! -z ${LAST_SOFTWARE_VERSION} ]; then
        if [[ "${LAST_SOFTWARE_VERSION}" != $VERSION ]]; then
            println "A new version of meta-index is available!\n" "${LAST_SOFTWARE_VERSION}"
            println "https://github.com/BlankenbergLab/meta-index/releases\n\n"
        fi
    fi
}

# Print credits
credits () {
    printf "Thanks for using meta-index!\n"
    printf "Please credit this tool in your manuscript by citing:\n\n"
    printf "\tTBA\n\n"

    printf "Remember to star the meta-index repository on GitHub to stay updated on its development and new features:\n"
    printf "https://github.com/BlankenbergLab/meta-index\n\n"
}

# Format seconds in human-readable format
# Credits: https://unix.stackexchange.com/a/27014
displaytime () {
    DAYS=$(bc <<< "${1}/60/60/24")
    HOURS=$(bc <<< "${1}/60/60%24")
    MINUTES=$(bc <<< "${1}/60%60")
    SECONDS=$(bc <<< "${1}%60")
    if [[ "$DAYS" -gt "0" ]]; then printf "%s days " "$DAYS"; fi
    if [[ "$HOURS" -gt "0" ]]; then printf "%s hours " "$HOURS"; fi
    if [[ "$MINUTES" -gt "0" ]]; then printf "%s minutes " "$MINUTES"; fi
    if [[ "$DAYS" -gt "0" ]] || [[ "$HOURS" -gt "0" ]] || [[ "$MINUTES" -gt "0" ]]; then printf "and "; fi
    printf "%s seconds\n" "$SECONDS"
}

# Get the last software version
# Credits: https://gist.github.com/lukechilds/a83e1d7127b78fef38c2914c4ececc3c
get_last_software_release () {
    curl --silent "https://api.github.com/repos/BlankenbergLab/meta-index/releases/latest" | \
        grep '"tag_name":' | \
        sed -E 's/.*"([^"]+)".*/\1/'
}

# Run CheckM for quality control base on completeness and contamination
# Remember to create a global variable "CHECKMTABLES" before calling this function 
# for accessing the result as a list of file paths separated by comma
run_checkm () {
    INLIST=$1       # Path to the file with the list of paths to the input genomes
    EXTENSION=$2    # Extension of input genome files
    TMPDIR=$3       # Temporary folder
    NPROC=$4        # Number of parallel processes for CheckM

    # Run CheckM
    printf 'Running CheckM\n'
    CHECKM_START_TIME="$(date +%s.%3N)"
    mkdir -p $TMPDIR/checkm/tmp
    # Split the set of bins in chunks with 1000 genomes at most
    println '\tOrganising genomes in chunks\n'
    split --numeric-suffixes=1 --lines=1000 --suffix-length=3 --additional-suffix=.txt $INLIST $TMPDIR/checkm/tmp/bins_
    CHUNKS=`ls "$TMPDIR"/checkm/tmp/bins_*.txt 2>/dev/null | wc -l`
    CHECKMTABLES=""
    for filepath in $TMPDIR/checkm/tmp/bins_*.txt; do
        # Retrieve chunk id
        filename="$(basename $filepath)"
        suffix="${filename#*_}"
        suffix="${suffix%.txt}"
        println '\tProcessing chunk %s/%s\n' "$((10#$suffix))" "$CHUNKS"
        # Create chunk folder
        mkdir -p $TMPDIR/checkm/tmp/bins_${suffix}
        for bin in `sed '/^$/d' $filepath`; do
            # Make a symbolic link to the 1000 genomes for the current chunk
            ln -s $bin $TMPDIR/checkm/tmp/bins_${suffix}
        done
        # Create the CheckM folder for the current chunk
        mkdir -p $TMPDIR/checkm/run_$suffix/
        # Run CheckM
        checkm lineage_wf -t ${NPROC} \
                          -x ${EXTENSION} \
                          --pplacer_threads ${NPROC} \
                          --tab_table -f $TMPDIR/checkm/run_${suffix}.tsv \
                          $TMPDIR/checkm/tmp/bins_${suffix} $TMPDIR/checkm/run_$suffix
        # Take trace of the CheckM output tables
        CHECKMTABLES=$TMPDIR/checkm/run_${suffix}.tsv,$CHECKMTABLES
    done
    CHECKM_END_TIME="$(date +%s.%3N)"
    CHECKM_ELAPSED="$(bc <<< "${CHECKM_END_TIME}-${CHECKM_START_TIME}")"
    println '\tTotal elapsed time: %s\n\n' "$(displaytime ${CHECKM_ELAPSED})"
    # Trim the last comma out of the list of CheckM output table file paths
    CHECKMTABLES="${CHECKMTABLES%?}"
}

# Transpose matrix file
# Credits: https://stackoverflow.com/a/28167793
transpose () {
  awk '{for (i=1; i<=NF; i++) a[i,NR]=$i; max=(max<NF?NF:max)}
        END {for (i=1; i<=max; i++)
              {for (j=1; j<=NR; j++) 
                  printf "%s%s", a[i,j], (j<NR?OFS:ORS)
              }
        }'
}