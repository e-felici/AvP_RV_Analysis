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



echo "************* BLASTp Databases Creation  *************"

#Define the directory for the databases. This is the directory where each folder has a fasta file for the database creation
DBDIR="/path/to/your/directory"

cd "$DBDIR"

for subdir in *
do 

	echo "Create database for $subdir"
	cd $subdir

	makeblastdb -in *.faa -out "$subdir"_db -dbtype prot -title "$subdir"_db -parse_seqids

	cd "$DBDIR"
done
