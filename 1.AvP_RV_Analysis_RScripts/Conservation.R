#!/usr/bin/Rscript

#load library
library("dplyr")
library("tidyr")
library("readr")
library("purrr")
library("tibble")
library("stringr")
library("forcats")
library("magrittr")
library("seqinr")
library("BALCONY")

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)  # Capture arguments
if (length(args) != 3) {
  stop("Three arguments are required!")
}

ConsDir <- args[1]
FinalRes <- args[2]
WorkDir <- args[3]


cat("Production of entropy (conservation) and strain count results\n")

tryCatch({
#Retrieve all needed files 
#IDs
ID_files = list.files(paste0(ConsDir, "/Mafft_Results/IDs/"), full.names = F)

#Alignments in fasta format
Aln_files = list.files(paste0(ConsDir, "/Mafft_Results/Fasta/"), full.names = FALSE)
Aln_filePath = paste0(ConsDir, "/Mafft_Results/Fasta/", Aln_files)

subfolders <- list.dirs(WorkDir, full.names = TRUE, 
                        recursive = FALSE)

# Count the number of subfolders
num_subfolders <- length(subfolders)

#We assume that 90% of strains is a reasonable amount of strains
num_subfolders <- round(90 * num_subfolders / 100, digits = 0)

#Initialize variable i and vector conservation
i = 1
Tibble <- tibble(Cluster_Number = character(), Mean_Schneider_Entropy = double())

#Entropy analysis for all files
for (i in 1:length(Aln_filePath)) {
      #Reading the alignment data in FASTA format
      alignment = seqinr::read.alignment(file = Aln_filePath[i], 
                                         format = "fasta", 
                                         forceToLower = FALSE)
  
      #Calculating the weights of the sequences based on Henikoff and Henikoff position-based method
      sequences_weights = get_pos_based_seq_weights(alignment) 
                          
      #Calculating the entropy of each of the multiple sequence alignment positions. 
      #(Shannon normalized or Schneider's)
      schneider_cons = schneider_conservativity(alignment = alignment, 
                                                weights = sequences_weights)
      
      #Create tibble
      Temp_Tibble <- tibble(Cluster_Number = Aln_files[i], Mean_Schneider_Entropy = mean(schneider_cons))
      
      #Join
      Tibble <- full_join(Tibble, Temp_Tibble, by = c("Cluster_Number", "Mean_Schneider_Entropy"))
      
      #Clean
      rm(alignment)
      rm(Temp_Tibble)
}

i=1
Tibble2 <- tibble(Cluster_Number = character(), Strain_count = integer())

# Apply the count_proteins function to each FASTA file
for (i in 1:length(Aln_filePath)) {
  # Read the file line by line
  lines <- read_lines(Aln_filePath[i])
  
  # Count lines that start with ">"
  count <- sum(str_starts(lines, ">"))
  
  #Create tibble
  Temp_Tibble <- tibble(Cluster_Number = Aln_files[i], Strain_count = count)
  
  #Join
  Tibble2 <- full_join(Tibble2, Temp_Tibble, by = c("Cluster_Number", "Strain_count"))
  
  #Clean
  rm(count)
  rm(Temp_Tibble)
}

Tibble <- full_join(Tibble, Tibble2, by = "Cluster_Number")

Tibble[[1]] <- gsub("-aligned.faa", "", Tibble[[1]])

# Extract all numeric values from the 'ID_files' object using a regular expression.
IDs = ID_files %>% str_extract(regex("[:digit:]+"))

# Initialize variable i and final_results tibble
i = 1 
final_results <- tibble(ID = character(), 
                        Mean_Schneider_Entropy = double(), 
                        Strain_count = integer(),
                        Cluster_Number = character())

# Loop through each ID in the 'IDs' vector.
for (i in i:length(IDs)) {
  
  # Check if Mean_Schneider_Entropy contains the current ID from the 'IDs' vector.
  position = Tibble$Cluster_Number %in% IDs[i]
  
  # Extract values for that cluster
  entropy = Tibble[position, 2] %>% unlist() %>% unname()
  count = Tibble[position, 3] %>% unlist() %>% unname()
  cluster = Tibble[position, 1] %>% unlist() %>% unname()
  
  # Open the corresponding file and read its contents as a TSV file.
  Temp_Tibble = read_tsv(paste0(ConsDir,"/Mafft_Results/IDs/",ID_files[i]), col_names = FALSE)
  
  colnames(Temp_Tibble)[1] = c("ID")
  
  # Add new columns with extracted numbers
  Temp_Tibble$Mean_Schneider_Entropy = entropy
  Temp_Tibble$Strain_count = count
  Temp_Tibble$Cluster_Number = cluster
  
  final_results <- full_join(final_results, Temp_Tibble, by = c("ID", "Mean_Schneider_Entropy", 
                                                            "Strain_count", "Cluster_Number"))
  
  #Clean
  rm(Temp_Tibble)
  rm(entropy)
  rm(count)
  rm(cluster)
  
  }
final_results$Conservation = 1 - final_results$Mean_Schneider_Entropy

final_results$Conservation_Results <- ifelse(final_results$Strain_count > num_subfolders, 

					     ifelse(final_results$Conservation > 0.95,
						    "Conservation Score > 0.95, more than 90% of the strains",
						    ifelse(final_results$Conservation > 0.8,
							"Conservation Score > 0.80, more than 90% of the strains",
							"Conservation Score < 0.80, more than 90% of the strains"
							  )
						    ),

					     ifelse(final_results$Conservation > 0.95,
						    "Conservation Score > 0.95, less than 90% of the strains",
						    ifelse(final_results$Conservation > 0.8,
							   "Conservation Score > 0.80, less than 90% of the strains",
							   "Conservation Score < 0.80, less than 90% of the strains"
						          )
						   )
					     )



final_results$ID <- final_results$ID %>% str_replace("\\.1", "")

  final_results <- final_results %>% 
  distinct() %>%
  arrange(ID)
  
# Save the final results
write_tsv(final_results, paste0(FinalRes,"/Conservation_results.tsv"))

rm(list = ls())

}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging conservation results: ", conditionMessage(e)))
})
