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

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)  # Capture arguments
if (length(args) != 1) {
  stop("Path to AgProtect argument is required")
}

AgProtect <- args[1]

cat("Filtering, cleaning, and sorting VFDB results\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw <- read_tsv(paste0(AgProtect,"/VFDB_full_results/AgProtect-vs_VFDB.out"), col_names = 1:12)
      #Protein IDs
  ids <- read_tsv(paste0(AgProtect,"/AllProteinIds-AgProtect.txt"), col_names = "ID")

  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < 1e-5) + (raw$X12>50) + (raw$X3 > 25)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble()
  
  # Polish data
  colnames(out)[1:5] = c("ID", "Virulence.Homologue","Virulence.IdentityPercent",
                         "Virulence.Evalue","Virulence.Bitscore")
  
  # Create data frame with proteins that were not classified as VF
  MissingIDs <- ids[!(ids$ID %in% out$ID), ]
  df <- data.frame(matrix("-", nrow = nrow(MissingIDs), ncol = 4))
  colnames(df) <- c("Virulence.Homologue","Virulence.IdentityPercent",
                    "Virulence.Evalue","Virulence.Bitscore")
  MissingIDs <- cbind(MissingIDs, df)
  
  # Combine the results with the missing proteins
  AllResults <- rbind(out, MissingIDs)
  
  # Arrange the combined data frame. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID, Virulence.Bitscore) %>%
    distinct(ID, .keep_all = TRUE)
  
  # Create virulence column
  AllResults$VirulenceFactor = ifelse(AllResults$Virulence.Homologue == "-",
                                      "Non Virulence Factor",
                                      "Probable Virulence Factor")
  
  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/AgProtect-vs_VFDB_final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging VFDB results: ", conditionMessage(e)))
})
