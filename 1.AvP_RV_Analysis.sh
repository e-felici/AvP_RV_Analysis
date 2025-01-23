#!/bin/bash

# =====================================================================
# Script Name:    1.AvP_RV_Analysis.sh
# Description:   
# Author:         M. Esperanza Felici
# Usage:          ./1.AvP_RV_Analysis.sh <path/to/file>
# =====================================================================


####Exit on error, unset variables, or pipe failures
set -euo pipefail


####Functions
# Check if required dependencies are installed
check_dependencies() {
    local dependencies=("blastp" "Rscript" "seqkit" "mafft" "cd-hit" "signalp6" "pepstats" "biolib" "rpsblastp" "python" "R")
    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" &>/dev/null; then
            echo "Error: $dep is not installed or not in PATH." >&2
            exit 1
        fi
    done
    log_message "All dependencies are installed."
}

#Print a log message with a timestamp
log_message() {
printf "[%s] %s\n" "$(date +'%Y-%m-%d %H:%M:%S')" "$1"
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

    # Prepare fasta for the final step (Conservation among strains)
    local specified_fasta="$strain_dir/SpecifiedStrain-$subdir_name.faa"
    sed "s/^>.*/& $subdir_name/" "$protein_file" > "$specified_fasta"
    mkdir -p "$ConsDir/Sequences/"
    cp "$specified_fasta" "$ConsDir/Sequences/"

    # Run R script for output processing
    log_message  "Running R script for $subdir_name..."
    Rscript "$RScripts/Info.R" "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }
}

homology_analysis() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Homology analysis with BLAST--------'

    # BLASTp: chicken proteome vs proteome of strain $subdir
    blastp -evalue 100 \
           -query "$WorkDir"/"$subdir_name"/protein.faa \
           -db "$DBDIR"/Gallus_gallus/Gallus_gallus_db \
           -out "$WorkDir"/"$subdir_name"/Homology_Analysis_results/"$subdir_name"-vs_chicken.out \
           -outfmt 6  

    # Run R script for output processing
    Rscript "$RScripts"/Homology_Analysis.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    } 
}

essential() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Search for essential proteins with DEG--------'

    # BLASTp: chicken proteome vs proteome of strain $subdir
    blastp -evalue 100 \
           -query "$WorkDir"/"$subdir_name"/protein.faa \
           -db "$DBDIR"/DEG_bacteria/DEG_bacteria_db \
           -out "$WorkDir"/"$subdir_name"/DEG_results/"$subdir_name"-vs_DEG.out \
           -outfmt 6  

    # Run R script for output processing
    Rscript "$RScripts"/DEG.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    } 
}
    
virulence() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Identification of virulence factors with VFDB--------'

    # BLASTp: chicken proteome vs proteome of strain $subdir
    blastp -evalue 100 \
           -query "$WorkDir"/"$subdir_name"/protein.faa \
           -db "$DBDIR"/VFDB_full/VFDB_full_db \
           -out "$WorkDir"/"$subdir_name"/VFDB_results/"$subdir_name"-vs_VFDB.out \
           -outfmt 6  

    # Run R script for output processing
    Rscript "$RScripts"/VFDB.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
}

characterization() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Protein characterization with pepstats (EMBOSS)--------'
    #executing pepstats (EMBOSS) 
    pepstats \
    	-sequence "$WorkDir"/"$subdir_name"/protein.faa \
    	-outfile "$WorkDir"/"$subdir_name"/EMBOSS_results/EMBOSS-$subdir 

    # Run R script for output processing
    Rscript "$RScripts"/EMBOSS.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
}

localization() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Subcellular localization prediction with PSORTb--------'
    #Starting analysis with PSORTb
    $Psortb \
 	-i "$WorkDir"/"$subdir_name"/protein.faa \
 	-r "$WorkDir"/"$subdir_name"/PSORTb_results/ \
 	-n \
 	-o terse \
 	-v

    # Run R script for output processing
    Rscript "$RScripts"/PSORTb.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
}


signal_peptide() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")


    log_message  '--------Signal peptides prediction with SignalP--------'
    #Starting analysis with SignalP
    signalp6 \
    --fastafile "$WorkDir"/"$subdir_name"/protein.faa \
    --organism other \
    --output_dir "$WorkDir"/"$subdir_name"/SignalP_results/ \
    --format txt \
    --mode fast
    
    # Run R script for output processing
    Rscript "$RScripts"/SignalP.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
}

adhesin_identification() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")

    log_message  '--------Adhesin identification with SPAAN --------'

    cp $WorkDir/$subdir_name/protein.faa "$SPAAN"
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
    cp SPAAN-unpolished.txt $WorkDir/$subdir_name/SPAAN_results/

    # Run R script
    if [[ -f "$RScripts/SPAAN.R" ]]; then
        Rscript "$RScripts"/SPAAN.R "$WorkDir" "$subdir_name"  || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
    else
        log_message  "Error: R script $RScripts/SPAAN.R not found."
        return 1
    fi
    
    popd
}

vaxijen() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")
    
    log_message  '--------Search for antigenic proteins with VaxiJen --------' 
    
    # Identify rare amino acids (selenocysteine and unknown residues)
    seqkit locate -d -i -p X "$WorkDir"/"$subdir_name"/protein.faa | awk '{print $1}' > "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_X
    seqkit locate -p U "$WorkDir"/"$subdir_name"/protein.faa | awk '{print $1}' > "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_Selenocysteine
    cat "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_Selenocysteine "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_X | grep -v "seqID" > "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_RareAminoAcids
    
    # Remove sequences with rare amino acids
    seqkit grep -i -v -f "$WorkDir"/"$subdir_name"/VaxiJen_results/IDs_RareAminoAcids "$WorkDir"/"$subdir_name"/protein.faa -o "$WorkDir"/"$subdir_name"/VaxiJen_results/$subdir_name-protein_filtered.faa
    
    # Split FASTA file into smaller chunks for VaxiJen3 submission
    seqkit split -s 100 -f "$WorkDir"/"$subdir_name"/VaxiJen_results/$subdir_name-protein_filtered.faa
    
    # Process each split file
    pushd "$WorkDir"/"$subdir_name"/VaxiJen_results/$subdir_name-protein_filtered.faa.split
    for file in *; do
        # Define maximum retries
        max_retries=3
        attempt=1
        success=false
    
        while [[ $attempt -le $max_retries ]]; do
            # Run VaxiJen.py script
            python "$VaxiJen" "$Chromedriver" "$WorkDir/$subdir_name/VaxiJen_results/$subdir_name-protein_filtered.faa.split/" "$file"
        
            # Check if the output file was created
            output_file="$WorkDir/$subdir_name/VaxiJen_results/$subdir_name-protein_filtered.faa.split/${file}-processed"
            if [[ -f "$output_file" ]]; then
                log_message "Successfully processed $file on attempt $attempt"
                success=true
                break
            else
                log_message "ERROR: Processing failed for $file on attempt $attempt. Retrying..."
                sleep 5
                ((attempt++))
            fi
        done
    
        if [[ $success == false ]]; then
            log_message "ERROR!! Processing failed for $file after $max_retries attempts. Exiting!"
            exit 1
        fi
    done
    
    #Merge all files
    cat "$WorkDir"/"$subdir_name"/VaxiJen_results/$subdir_name-protein_filtered.faa.split/* > "$WorkDir"/"$subdir_name"/VaxiJen_results/VaxiJen3_predictions.csv
    
    # Run R script for output processing
    Rscript "$RScripts"/VaxiJen.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    }    
    popd
}

deeptmhmm() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")
    
    log_message  '-------- Transmembrane Topology Prediction and Classification with DeepTMHMM --------'
    
    cp "$WorkDir"/"$subdir_name"/protein.faa "$WorkDir"/"$subdir_name"/DeepTMHMM_results/
    
    pushd "$WorkDir"/"$subdir_name"/DeepTMHMM_results/
    
    # Split the fasta file into smaller chunks
    seqkit split -s 300 -f protein.faa || exit
    
    # Define batch size for processing
    batch_size=4
    files_to_process=("$WorkDir"/"$subdir_name"/DeepTMHMM_results/protein.faa.split/*)
    
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
        	mv biolib_results/predicted_topologies* "$WorkDir"/"$subdir_name"/DeepTMHMM_results/predicted_topologies-$file_number
        	mv biolib_results/TMRs.gff3* "$WorkDir"/"$subdir_name"/DeepTMHMM_results/TMRs.gff3-$file_number
        	
        	# Increment file number and remove processed file
        	((file_number++))
        	rm -f "$file_to_process"
        done
        #Pausing script for 22 hours after processing another batch of files. This is necessary due to the sequence submission limits of DeepTMHMM
        now="$(date)"
        log_message  "Pause script for 22 hours, starting from: $now"
        sleep 22h

     done
     
        
     # Run R script for output processing
    Rscript "$RScripts"/DeepTMHMM.R "$WorkDir" "$subdir_name" || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    } 
    popd   
}

COG() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")

    log_message  '--------COG Analysis with RPS-BLAST--------'

    #Running rps-blast
    rpsblastp -evalue 100 \
           -query "$WorkDir"/"$subdir_name"/protein.faa \
           -db "$COG"/Cog_LE/Cog \
           -out "$WorkDir"/"$subdir_name"/COG_results/"$subdir_name"-vs_COG.out \
           -outfmt 6  

    # Run R script for output processing
    Rscript "$RScripts"/COG.R "$WorkDir" "$subdir_name" "$COG" || {
        log_message  "ERROR: R script failed for $subdir_name"
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
        2) MultiSeq="$line" ;;
        3) ConsDir="$line" ;;
        4) FinalRes="$line" ;;
        5) DBDIR="$line" ;;
        6) RScripts="$line" ;;
        7) Psortb="$line" ;;
        8) SPAAN="$line" ;;
        9) VaxiJen="$line" ;;
        10) Chromedriver="$line" ;;
        11) COG="$line" ;;
	12) ;;
        13) ;;
        14) ;;
	15) ;;
        *) log_message  "Warning: File contains too many lines. Ignoring extra lines." ;;
    esac
done < "$file"

# Verify that all arguments are set
if [ -z "$WorkDir" ] || [ -z "$MultiSeq" ] || [ -z "$ConsDir" ] || [ -z "$FinalRes" ] || [ -z "$DBDIR" ] || [ -z "$RScripts" ] || [ -z "$Psortb" ] || [ -z "$SPAAN" ] || [ -z "$VaxiJen" ] || [ -z "$Chromedriver" ] || [ -z "$COG" ] ; then
    log_message  "Error: The file must contain all the arguments!"
    exit 1
fi

# Print the arguments 
log_message  'These are the paths that you defined:'
log_message  "Path to strains folder (where each strain is a folder, and inside each folder there is the proteome of that strain (protein.faa file)): $WorkDir"
log_message  "Path to make_multi_seq.pl script from CD_HIT: $MultiSeq"
log_message  "Path to folder in which all the conservation among strains results will be stored: $ConsDir"
log_message  "Path to the final and polished results: $FinalRes"
log_message  "Path to the blastp databases that you have created with 0.DB_creation.sh script: $DBDIR"
log_message  "Path to folder with all the R Scripts to analyze the Av. paragallinarum strains (1.AvP_RV_Analysis_RScripts folder): $RScripts"
log_message  "Path to Perl script wrapper to run a docker run to execute PSORTb inside the docker container: $Psortb"
log_message  "Path to SPAAN folder with all the files to run standalone SPAAN: $SPAAN"
log_message  "Path to VaxiJen.py script: $VaxiJen"
log_message  "Path to Chromedriver binary: $Chromedriver"
log_message  "Path to COG folder with all COG and CDD resources: $COG"

pushd "$WorkDir"
for strain_dir in "$WorkDir"/*; do
    if [[ -d "$strain_dir" ]]; then
        # Checking if fasta exists
    	if [[ ! -f $strain_dir/protein.faa ]]; then
        	log_message  "Error: File $strain_dir/protein.faa not found!"
        	return 1
    	fi

    	#####Part 2: Pre-process each strain folder
        prepare_analysis "$strain_dir" || { log_message  "Error in pre-processing proteins in $strain_dir"; return 1; }

        #####Part 3: Homology analysis
        homology_analysis "$strain_dir" || { log_message  "Error in homology analysis for proteins in $strain_dir"; return 1; }

        #####Part 4: Search for essential proteins with DEG
        essential "$strain_dir" || { log_message  "Error in essentiality analysis for proteins in $strain_dir"; return 1; }

        #####Part 5: Identification of virulence factors with VFDB
        virulence "$strain_dir" || { log_message  "Error in virulence analysis for proteins in $strain_dir"; return 1; }

        #####Part 6: Protein characterization with pepstats (EMBOSS)
        characterization "$strain_dir" || { log_message  "Error in characterization of proteins in $strain_dir"; return 1; }

        #####Part 7: Subcellular localization prediction with PSORTb
        localization "$strain_dir" || { log_message  "Error in finding subcellular localization of proteins in $strain_dir"; return 1; }

        #####Part 8: Signal peptides prediction with SignalP
        signal_peptide "$strain_dir" || { log_message  "Error in finding signal peptides of proteins in $strain_dir"; return 1; }
        
        #####Part 9: Adhesin Identification with SPAAN
        adhesin_identification "$strain_dir" || { log_message  "Adhesin identification failed for $strain_dir"; return 1; }
        
        #####Part 10: Search for antigenic proteins with VaxiJen
        vaxijen "$strain_dir" || { log_message  "Antigenic proteins identification failed for $strain_dir"; return 1; }
        
        #####Part 11: Transmembrane Topology Prediction and Classification with DeepTMHMM
        deeptmhmm "$strain_dir" || { log_message  "Transmembrane Topology Prediction failed for $strain_dir"; return 1; }
        
        #####Part 12: COG Analysis with RPS-BLAST
        COG "$strain_dir" || { log_message  "COG Analysis with RPS-BLAST failed for $strain_dir"; return 1; }
        
        #####Part 13: Run R script for final output processing
        Rscript "$RScripts"/Final_polishing.R "$WorkDir" "$strain_dir" || {
        	log_message  "ERROR: R script failed for $strain_dir"
        	exit 1
    	}
    	cp "$WorkDir"/$strain_dir/Final_results-$strain_dir.tsv "$FinalRes"

        popd || { log_message  "Error: Couldn't return to $WorkDir"; return 1; }
        
    else
        log_message  "Skipping $strain_dir (not a directory)"
    fi
done

#####Part 14: Creating al necessary folders in ConsDir 
log_message  "-------- Creating al necessary folders --------"

pushd $ConsDir
mkdir -p CD-Hit_Results \
         Mafft_Results\
         Other_CD-Hit_Results
popd

pushd $ConsDir/Mafft_Results
mkdir -p Fasta \
         Clustal \
         IDs
popd

pushd $ConsDir/Sequences
cat *faa > All_Sequences.faa

#####Part 15:Clustering the proteins with CD-HIT 
log_message  "-------- Clustering the proteins with CD-HIT --------"
cd-hit -i All_Sequences.faa -o "$ConsDir"/Other_CD-Hit_Results/All_Sequences_clusters.out -d 0 -g 1 -sf 1
# Explanation of flags used in CD-HIT
# -d    0 uses sequence name in FASTA header up to the first whitespace
# -sf   sort FASTA/FASTQ by cluster size (number of sequences), default 0 (no sorting). If set to 1, output sequences by decreasing cluster size
# -g    1 or 0, default 0. By cd-hit's default algorithm, a sequence is clustered to the first cluster meeting the threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster meeting the threshold (accurate but slower mode), but either 1 or 0 won't change the representatives of the final clusters.
#Note: to use make_multi_seq.pl script, "-d 0" and "-g 1" options should be used

#Read the .clstr file and generate a seperate fasta file for each cluster over certain size (in this case, clusters with 1 or more members will be generated)
$MultiSeq All_Sequences.faa "$ConsDir"/Other_CD-Hit_Results/All_Sequences_clusters.out.clstr multi-seq 1    

mv multi-seq/* "$ConsDir"/CD-Hit_Results

#Cleaning
rm -r multi-seq    
popd

#####Part 16: Multiple alignment with Mafft
log_message  "-------- Multiple alignment with Mafft --------"
pushd "$ConsDir"/CD-Hit_Results

for file in *; 
do
	mafft --localpair --maxiterate 1000 --clustalout --anysymbol $file > "$ConsDir"/Mafft_Results/Clustal/$file-aligned.aln
	mafft --localpair --maxiterate 1000 --anysymbol $file > "$ConsDir"/Mafft_Results/Fasta/$file-aligned.faa
done
popd

pushd "$ConsDir"/Mafft_Results/Clustal

#Extracting IDs
for file in *; 
do
	log_message  "Extracting IDs from $file"
	awk '{print $1}' $file | grep -v --binary-files=text "CLUSTAL" | grep -v --binary-files=text "*" | grep -v --binary-files=text ":" | grep -v --binary-files=text "                        " > TEMP-file
	mv TEMP-file "$ConsDir"/Mafft_Results/IDs/$file
done
popd

#####Part 17: Run R script for final conservation output processing
#Executing R script
Rscript "$RScripts"/Conservation.R "$ConsDir" "$FinalRes" "$WorkDir"


log_message  "************* End of Analysis *************"
