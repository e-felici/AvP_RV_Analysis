#!/bin/bash


echo "************* BLASTp Databases Creation  *************"

#Before performing blastp, We must create the databases. See: https://www.ncbi.nlm.nih.gov/books/NBK279690/ 
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
