#!/bin/bash

# =====================================================================
# Script Name:    2.AgProtect_RV_Analysis.sh
# Description:   
# Author:         M. Esperanza Felici
# Usage:          2.AgProtect_RV_Analysis.sh <path/to/file>
# =====================================================================


####Exit on error, unset variables, or pipe failures
set -euo pipefail


####Functions
# Check if required dependencies are installed
check_dependencies() {
    local dependencies=("blastp" "Rscript" "seqkit" "signalp6" "pepstats" "biolib" "rpsblastp" "python" "R")
    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" &>/dev/null; then
            echo "Error: $dep is not installed or not in PATH." >&2
            exit 1
        fi
    done
    echo "All dependencies are installed."
}

#Print a log message with a timestamp
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

#Prepare everything for analysis of a single strain
prepare_analysis() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")

    log_message  "~~~~~~~~~~ Preparation for Analysis of $subdir_name ~~~~~~~~~~"

    # Create subfolders in the strain directory
    log_message  "Creating subfolders for $subdir_name..."
    mkdir -p "$strain_dir/Homology_Analysis_results" \
             "$strain_dir/DEG_results" \
             "$strain_dir/VFDB_full_results" \
             "$strain_dir/Final_results" \
             "$strain_dir/PSORTb_results" \
             "$strain_dir/VaxiJen_results" \
             "$strain_dir/EMBOSS_results" \
             "$strain_dir/SignalP_results" \
             "$strain_dir/SPAAN_results" \
             "$strain_dir/DeepTMHMM_results" \
             "$strain_dir/COG_results"

    # Build the complete list of total protein IDs
    local protein_file="$strain_dir/protein.faa"
    if [[ -f "$protein_file" ]]; then
        grep ">" --binary-files=text "$protein_file" | awk '{print $1}' > "$strain_dir/AllProteinIds-$subdir_name.txt"
        sed -i 's/>//g' "$strain_dir/AllProteinIds-$subdir_name.txt"
    else
        log_message  "WARNING: Protein file not found in $strain_dir. Skipping..."
        return
    fi

    # Run R script for output processing
    log_message  "Running R script for AgProtect..."
    Rscript "$RScripts/Info_AgProtect.R" "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }
}


blast_against_databases() {
    echo '--------BLAST against different databases--------'
    
    # Navigate to the directory containing the databases
    pushd "$DBDIR" || { echo "Error: Could not navigate to $DBDIR"; return 1; }

    # Loop through each folder in the database directory
    for DBfolder in */; do
        # Remove trailing slash
        DBfolder=${DBfolder%/}

        echo "Processing BLAST for database: $DBfolder"

        # Run BLASTP command
        blastp -evalue 100 \
            -query $AgProtect/protein.faa \
            -db $DBfolder/${DBfolder}_db \
            -out $AgProtect/Homology_Analysis_results/$DBfolder.out \
            -outfmt 6
    done
    
    popd
}

characterization() {
    log_message  '--------Protein characterization with pepstats (EMBOSS)--------'
    #executing pepstats (EMBOSS) 
    pepstats \
    	-sequence "$AgProtect"/protein.faa \
    	-outfile "$AgProtect"/EMBOSS_results/EMBOSS-AgProtect 

    # Run R script for output processing
    Rscript "$RScripts"/EMBOSS.R "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }    
}

localization() {
    log_message  '--------Subcellular localization prediction with PSORTb--------'
    
    $Psortb \
 	-i "$AgProtect"/protein.faa \
 	-r "$AgProtect"/PSORTb_results/ \
 	-n \ #Gram negative
 	-o terse \
 	-v
    mv "$AgProtect"/PSORTb_results/*gramneg.txt "$AgProtect"/PSORTb_results/GramNegative

    $Psortb \
 	-i "$AgProtect"/protein.faa \
 	-r "$AgProtect"/PSORTb_results/ \
 	-p \ #Gram positive
 	-o terse \
 	-v
    mv "$AgProtect"/PSORTb_results/*grampos.txt "$AgProtect"/PSORTb_results/GramPositive
    
    mv "$GramAdv"/GramPositiveWithOM "$AgProtect"/PSORTb_results/GramPositiveWithOM
    mv "$GramAdv"/GramNegativeWithoutOM "$AgProtect"/PSORTb_results/GramNegativeWithoutOM
    
    # Run R script for output processing
    Rscript "$RScripts"/PSORTb.R "$AgProtect"  || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }    
}

signal_peptide() {
    log_message  '--------Signal peptides prediction with SignalP--------'
    #Starting analysis with SignalP
    signalp6 \
    --fastafile "$AgProtect"/protein.faa \
    --organism other \
    --output_dir "$AgProtect"/SignalP_results/ \
    --format txt \
    --mode fast
    
    # Run R script for output processing
    Rscript "$RScripts"/SignalP_AgProtect.R "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }    
}

adhesin_identification() {
    log_message  '--------Adhesin identification with SPAAN --------'

    cp "$AgProtect"/protein.faa "$SPAAN"
    pushd "$SPAAN" || { log_message  "Error: Could not change directory to $SPAAN"; return 1; }
    mv protein.faa query.dat 

    # Executing SPAAN
    if [[ ! -x ./askquery ]]; then
        log_message  "Error: SPAAN executable 'askquery' not found or not executable."
        return 1
    fi
    ./askquery

    # Cleaning and copying results
    mv query.out SPAAN-unpolished.txt 
    rm query.dat
    cp SPAAN-unpolished.txt "$AgProtect"/SPAAN_results/

    # Run R script
    if [[ -f "$RScripts/SPAAN_AgProtect.R" ]]; then
        Rscript "$RScripts/SPAAN_AgProtect.R" "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }    
    else
        log_message  "Error: R script $RScripts/SPAAN_AgProtect.R not found."
        return 1
    fi
    
    popd
}

vaxijen() {
    log_message  '--------Search for antigenic proteins with VaxiJen --------' 
    
    # Identify rare amino acids (selenocysteine and unknown residues)
    seqkit locate -d -i -p X "$AgProtect"/protein.faa | awk '{print $1}' > "$AgProtect"/VaxiJen_results/IDs_X
    seqkit locate -p U "$AgProtect"/protein.faa | awk '{print $1}' > "$AgProtect"/VaxiJen_results/IDs_Selenocysteine
    cat "$AgProtect"/VaxiJen_results/IDs_Selenocysteine "$AgProtect"/VaxiJen_results/IDs_X | grep -v "seqID" > "$AgProtect"/VaxiJen_results/IDs_RareAminoAcids
    
    # Remove sequences with rare amino acids
    seqkit grep -i -v -f "$AgProtect"/VaxiJen_results/IDs_RareAminoAcids "$AgProtect"/protein.faa -o "$AgProtect"/VaxiJen_results/AgProtect-protein_filtered.faa
    
    # Split FASTA file into smaller chunks for VaxiJen3 submission
    seqkit split -s 100 -f "$AgProtect"/VaxiJen_results/AgProtect-protein_filtered.faa
    
    # Process each split file
    pushd "$AgProtect"/VaxiJen_results/AgProtect-protein_filtered.faa.split
    for file in *; do
    	log_message  "Processing $file"
    	# Run VaxiJen.py script
    	$VaxiJen "$Chromedriver" "$AgProtect"/VaxiJen_results/AgProtect-protein_filtered.faa.split/ "$file"
    done
    
    #Merge all files
    cat "$AgProtect"/VaxiJen_results/AgProtect-protein_filtered.faa.split/* > "$AgProtect"/VaxiJen_results/VaxiJen3_predictions.csv
    
    # Run R script for output processing
    Rscript "$RScripts"/VaxiJen_AgProtect.R "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    }    
    popd
}

deeptmhmm() {
    log_message  '-------- Transmembrane Topology Prediction and Classification with DeepTMHMM --------'
    
    cp "$AgProtect"/protein.faa "$AgProtect"/DeepTMHMM_results/
    
    pushd "$AgProtect"/DeepTMHMM_results/
    
    # Split the fasta file into smaller chunks
    seqkit split -s 300 -f protein.faa || exit
    
    # Define batch size for processing
    batch_size=4
    files_to_process=("$AgProtect"/DeepTMHMM_results/protein.faa.split/*)
    
    # Initialize file number
    file_number=1
    
    # Process files in batches
    for ((i = 0; i < ${#files_to_process[@]}; i += batch_size)); do
    	batch=("${files_to_process[@]:i:batch_size}")
    	log_message  "Processing batch $((i / batch_size + 1)) out of $(((${#files_to_process[@]} + batch_size - 1) / batch_size))"
	
	# Define maximum retries
	max_retries=3
	
    	# Process each file in the batch
    	for file_to_process in "${batch[@]}"; do
        	
        	attempt=1
        	success=false
        	
        	while [[ $attempt -le $max_retries ]]; do
        		log_message "Attempt $attempt: Running DeepTMHMM on $file_to_process"
        		biolib run DTU/DeepTMHMM --fasta "$file_to_process"

                        if [[ $? -eq 0 ]]; then
            			log_message "Successfully processed $file_to_process on attempt $attempt"
            			success=true
            			break
        		else
            			log_message "ERROR: DeepTMHMM failed on $file_to_process (attempt $attempt)"
            			((attempt++))
        		fi
   		done

    		if [[ $success == false ]]; then
        		log_message "ERROR: DeepTMHMM failed after $max_retries attempts for $file_to_process"
        		# Fail!!
        		exit 1 
    		fi

        	# Move results to the corresponding folder
        	mv biolib_results/predicted_topologies* "$AgProtect"/DeepTMHMM_results/predicted_topologies-$file_number
        	mv biolib_results/TMRs.gff3* "$AgProtect"/DeepTMHMM_results/TMRs.gff3-$file_number
        	
        	# Increment file number and remove processed file
        	((file_number++))
        	rm -f "$file_to_process"
        done
        #Pausing script for 22 hours after processing another batch of files. This is necessary due to the sequence submission limits of DeepTMHMM
        now="$(date)"
        log_message  "Pause script for 22 hours, starting from: $now"
        sleep 22h

     done
     
     #Extracting info about the protein type
     #Combining all predicted_topologies* files to get info about the protein type (Globular/Beta/Alpha with or without SP)
     cat "$AgProtect"/DeepTMHMM_results/predicted_topologies* > "$AgProtect"/DeepTMHMM_results/DeepTMHMM_type_AgProtect.txt
   
     #Extracting info about TMs and their characteristics
     #Create a Temporal Directory
     mkdir -p TEMP
     #Copy TMRs* files
     cp TMRs* TEMP/
     popd
     
     #Move to Temporal directory
     pushd "$AgProtect"/DeepTMHMM_results/TEMP
     # Loop through all files in the current directory
	for file in *
	do
		# Append '//' to the end of each file
    		printf '//
' >> "$file"
    		# Insert '//' at the beginning of each file
    		sed -i '1 i\\/\/' "$file"
    			
    		#Replace older TMRs* files
    		mv "$file" "$AgProtect"/DeepTMHMM_results
	done
     popd
     
     #Move back to subdir
     pushd "$AgProtect"/DeepTMHMM_results
     #Combining and cleaning all TMRs* files
     cat TMRs* | grep -v "##gff-version 3" --binary-files=text > TEMP_DeepTMHMM_TM_AgProtect.txt
     # Separate TM characteristics
     grep -v -e "Length" -e "Number of predicted TMRs" --binary-files=text TEMP_DeepTMHMM_TM_AgProtect.txt > TEMP_DeepTMHMM_TMbreakdown_AgProtect.txt
     #Splitting the file based on //
     csplit -z -s TEMP_DeepTMHMM_TMbreakdown_AgProtect.txt /\/\// '{*}'
     #Moving the split files to Temporal Directory
     mv xx* TEMP
     popd
     
     #Entering the folder with the split files
     pushd TEMP
     #Separating the pieces and counting them
     for file in * ;
     do
		grep -o "signal" $file | wc -l > $file-7.signal
		grep -o "periplasm" $file | wc -l > $file-3.peripl
    		grep -o "inside" $file | wc -l > $file-2.inside
    		grep -o "outside" $file | wc -l > $file-4.outside
   		grep -o "Beta sheet" $file | wc -l > $file-5.beta
    		grep -o "TMhelix" $file | wc -l > $file-6.alfa
    		grep -v "#" $file | grep -v "//" | awk '{print $1;exit;}' > $file-1.name
   		paste $file-* > $file-all
     done
     popd
     
     #Returning to the folder 
     pushd "$AgProtect"/DeepTMHMM_results
     #Merging the pieces back together and cleaning up "proteins" that do not exist (derived from ending and beggining of file)
     cat TEMP/*-all | grep -v --binary-files=text "  0  0  0  0  0  0" > DeepTMHMM_TMbreakdown_AgProtect.txt
     
     # Clean up temporary files
     rm TEMP*
     rm -r TEMP biolib_results
     
     # Run R script for output processing
     Rscript "$RScripts"/DeepTMHMM.R "$AgProtect" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
     } 
     popd   
}

COG() {
    log_message  '--------COG Analysis with RPS-BLAST--------'

    #Running rps-blast
    rpsblastp -evalue 100 \
           -query "$AgProtect"/protein.faa \
           -db "$COG"/Cog_LE/Cog \
           -out "$AgProtect"/COG_results/AgProtect-vs_COG.out \
           -outfmt 6  

    # Run R script for output processing
    Rscript "$RScripts"/COG.R "$AgProtect" "$COG" || {
        log_message  "ERROR: R script failed for AgProtect"
        exit 1
    } 
}

####Main Script

#####Part 1: Getting the paths 
# Perform the dependency check
check_dependencies

# Check if the file argument is provided
if [ "$#" -ne 1 ]; then
    log_message  "Usage: $0 <path/to/file>"
    exit 1
fi

file=$1

# Check if the file exists
if [ ! -f "$file" ]; then
    log_message  "File not found: $file"
    exit 1
fi

# Read the file line by line
count=0
while IFS= read -r line; do
    count=$((count + 1))
    case $count in
        1) WorkDir="$line" ;;
	2) ;;
	3) ;;  
        4) FinalRes="$line" ;;
        5) DBDIR="$line" ;;
        6) ;;
        7) Psortb="$line" ;;
        8) SPAAN="$line" ;;
        9) VaxiJen="$line" ;;
        10) Chromedriver="$line" ;;
        11) COG="$line" ;;
        12) AgProtect="$line" ;;
        13) RScripts="$line" ;;
        14) GramAdv="$line" ;;
	15) RScripts="$line" ;;
        *) log_message  "Warning: File contains too many lines. Ignoring extra lines." ;;
    esac
done < "$file"

# Verify that all arguments are set
if [ -z "$AgProtect" ] || [ -z "$WorkDir" ] || [ -z "$FinalRes" ] || [ -z "$DBDIR" ] || [ -z "$RScripts" ] || [ -z "$Psortb" ] || [ -z "$SPAAN" ] || [ -z "$VaxiJen" ] || [ -z "$Chromedriver" ] || [ -z "$COG" ] || [ -z "$GramAdv" ]; then
    log_message  "Error: The file must contain all the arguments!"
    exit 1
fi

# Print the arguments 
log_message  'These are the paths that you defined:'
log_message  "Path to experimental antigens (AgProtect) folder (inside there is the protein.faa file)): $AgProtect"
log_message  "Path to the final and polished results: $FinalRes"
log_message  "Path to the blastp databases that yoou have created with 0.DB_creation.sh script: $DBDIR"
log_message  "Path to folder with all the R Scripts to analyze AgProtect (2.AgProtect_RV_Analysis_RScripts folder): $RScripts"
log_message  "Path to Perl script wrapper to run a docker run to execute PSORTb inside the docker container: $Psortb"
log_message  "Path to SPAAN folder with all the files to run standalone SPAAN: $SPAAN"
log_message  "Path to VaxiJen.py script: $VaxiJen"
log_message  "Path to Chromedriver binary: $Chromedriver"
log_message  "Path to COG folder with all COG and CDD resources: $COG"
log_message  "Path to strains folder (where each strain is a folder, and inside each folder there is the proteome of that strain (protein.faa file)): $WorkDir"
log_message  "Path to GramPositiveWithOM and GramNegativeWithoutOM (manually downloaded from PSORTb sequence submission page): $GramAdv"


#start analysis
pushd "$AgProtect"

if [[ ! -f protein.faa ]]; then
       	log_message  "Error: File $AgProtect/protein.faa not found!"
       	return 1
fi

if [[ ! -f BacteriaTypes.tsv ]]; then
       	log_message  "Error: File $AgProtect/BacteriaTypes.tsv not found!"
       	return 1
fi

if [[ ! -f HostList.tsv ]]; then
       	log_message  "Error: File $AgProtect/HostList.tsv not found!"
       	return 1
fi

#####Part 2: Pre-process each strain folder
prepare_analysis "$AgProtect" || { log_message  "Error in pre-processing proteins from AgProtect"; return 1; }

#####Part 3: Blastp against several databases
blast_against_databases "$AgProtect" || { log_message  "Error in blastp analysis for proteins in AgProtect"; return 1; }
#Moving some files to their corresponding place
mv "$AgProtect"/Homology_Analysis_results/DEG_bacteria.out "$AgProtect"/Resultados_DEG/AgProtect-vs_DEG.out
mv "$AgProtect"/Homology_Analysis_results/VFDB_full.out "$AgProtect"/Resultados_DEG/AgProtect-vs_VFDB.out

#####Part 4: Homology analysis
Rscript "$RScripts"/Homology_Analysis_AgProtect.R "$AgProtect" || { log_message  "Error in homology analysis for proteins in AgProtect"; return 1; }

#####Part 5: Search for essential proteins with DEG
Rscript "$RScripts"/DEG_AgProtect.R "$AgProtect" || { log_message  "Error in essentiality analysis for proteins in AgProtect"; return 1; }

#####Part 6: Identification of virulence factors with VFDB
Rscript "$RScripts"/VFDB_AgProtect.R "$AgProtect" || { log_message  "Error in virulence analysis for proteins in AgProtect"; return 1; }

#####Part 7: Protein characterization with pepstats (EMBOSS)
characterization || { log_message  "Error in characterization of proteins in AgProtect"; return 1; }

#####Part 8: Subcellular localization prediction with PSORTb
localization || { log_message  "Error in finding subcellular localization of proteins in AgProtect"; return 1; }

#####Part 9: Signal peptides prediction with SignalP
signal_peptide  || { log_message  "Error in finding signal peptides of proteins in AgProtect"; return 1; }
        
#####Part 10: Adhesin Identification with SPAAN
adhesin_identification || { log_message  "Adhesin identification failed for AgProtect"; return 1; }
        
#####Part 11: Search for antigenic proteins with VaxiJen
vaxijen || { log_message  "Antigenic proteins identification failed for AgProtect"; return 1; }
        
#####Part 12: Transmembrane Topology Prediction and Classification with DeepTMHMM
deeptmhmm  || { log_message  "Transmembrane Topology Prediction failed for AgProtect"; return 1; }
        
#####Part 13: COG Analysis with RPS-BLAST
COG  || { log_message  "COG Analysis with RPS-BLAST failed for AgProtect"; return 1; }
        
#####Part 14: Run R script for final output processing
Rscript "$RScripts"/Final_polishing_AgProtect.R "$AgProtect" || {
      	log_message  "ERROR: R script failed for AgProtect"
       	exit 1
    	}
cp "$AgProtect"/Final_results-AgProtect.tsv "$FinalRes"/

        popd 
echo "-------- End of most of the analysis for strain AgProtect --------"
