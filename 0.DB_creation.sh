#!/bin/bash

# =====================================================================
# Script Name:    0.DB_creation.sh
# Description:    Before performing blastp, We must create the databases for all the blastp comparisons.
#                 So, this is what this script does.
#                 For the AvP antigen search, we created the databases for the fasta files stored in all
#                 these folders:
#                    -Bos_taurus
#                    -Felis_silvestris_catus
#                    -Canis_familiaris
#                    -Homo_sapiens
#                    -Ovis_orientalis_aries
#                    -Capra_aegagrus_hircus
#                    -Mus_musculus
#                    -Sus_scrofa
#                    -Equus_ferus_caballus
#                    -Oncorhynchus_mykiss
#                    -DEG_bacteria
#                    -VFDB_full
#                 For more details, see: https://www.ncbi.nlm.nih.gov/books/NBK279690/ 
# Author:         M. Esperanza Felici
# Usage:          ./0.DB_creation.sh
# =====================================================================

####Exit on error, unset variables, or pipe failures
set -euo pipefail


####Functions
usage() {
    echo "Usage: $0 <database_directory>"
    echo
    echo "This script iterates through subdirectories of a specified directory,"
    echo "searches for FASTA (*.faa) files, and creates BLAST protein databases with them."
    echo
    echo "Arguments:"
    echo "  database_directory   The directory containing subdirectories with FASTA files."
    echo
    exit 1
}

#Print a log message with a timestamp
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

#Process subdirectories, check for FASTA files (*.faa), and create a BLAST protein database from them
create_blast_db() {
    local subdir="$1"
    log_message "Processing subdirectory: $subdir"

    pushd "$subdir" > /dev/null

    # Check if there are any .faa files
    if ls *.faa > /dev/null 2>&1; then
        for fasta_file in *.faa; do
            local db_name="${fasta_file%.faa}_db"  # Use the fasta file name as a base for the database name
            log_message "Creating BLAST database for $subdir using $fasta_file..."
            if makeblastdb -in "$fasta_file" -out "$db_name" -dbtype prot -title "$db_name" -parse_seqids; then
                log_message "Successfully created BLAST database: $db_name"
            else
                log_message "Error creating BLAST database for $fasta_file in $subdir"
            fi
        done
    else
        log_message "No FASTA (*.faa) files found in $subdir. Skipping."
    fi

    popd > /dev/null
}

####Main Script
#Check if the user provided the required directory argument
if [[ $# -ne 1 ]]; then
    usage
fi

#Assign the first command-line argument ($1) to the variable DBDIR
DBDIR="$1"

#Directory Existence Check
if [[ ! -d "$DBDIR" ]]; then
    echo "Error: Directory $DBDIR does not exist."
    usage
fi

#Starting
log_message "Starting BLAST database creation in directory: $DBDIR"

#If an item is a subdirectory, it calls the create_blast_db() function to process it and create a BLAST database
for subdir in "$DBDIR"/*; do
    if [[ -d "$subdir" ]]; then
        create_blast_db "$subdir"
    else
        log_message "Skipping non-directory entry: $subdir"
    fi
done

#Ending
log_message "Finished creating BLAST databases."
