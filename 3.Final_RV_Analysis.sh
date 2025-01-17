#!/bin/bash


# =====================================================================
# Script Name:    3.Final_RV_Analysis.sh
# Description:   
# Author:         M. Esperanza Felici
# Usage:          ./3.Final_RV_Analysis.sh <path/to/file>
# =====================================================================

echo '-------- Final polishing --------'
# Check if the file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path/to/file>"
    exit 1
fi

file=$1

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "File not found: $file"
    exit 1
fi

# Read the file line by line
count=0
while IFS= read -r line; do
    count=$((count + 1))
    case $count in
        1) ;;
	2) ;;
	3) ;;  
        4) FinalRes="$line" ;;
        5) ;;
        6) ;;
        7) ;;
        8) ;;
        9) ;;
        10) ;;
        11) ;;
        12) AgProtect="$line" ;;
        13) ;;
        14) ;;
        15) RScripts="$line" ;;
        *) echo "Warning: File contains too many lines. Ignoring extra lines." ;;
    esac
done < "$file"

# Verify that all arguments are set
if [ -z "$AgProtect" ] ||  || [ -z "$FinalRes" ] || [ -z "$RScripts" ] ; then
    echo "Error: The file must contain all the arguments!"
    exit 1
fi

# Print the arguments 
echo 'These are the paths that you defined:'
echo "AgProtect folder: $AgProtect"
echo "Final Results folder: $FinalRes"
echo "3.Final_RV_Analysis_RScripts folder: $RScripts"


#Executing R script
Rscript "$RScripts"/Final_merge.R $FinalRes $AgProtect 

echo "************* End of Analysis *************"

