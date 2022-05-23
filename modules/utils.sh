#!/bin/bash
#title          :utils
#description    :Some utility functions
#author         :Fabio Cumbo (fabio.cumbo@gmail.com)
#===================================================

DATE="May 17, 2022"
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