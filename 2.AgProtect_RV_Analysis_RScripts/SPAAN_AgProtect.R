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

cat("Filtering, cleaning, and sorting the results of SPAAN\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw = read_tsv(paste0(AgProtect,"/SPAAN_results/SPAAN-unpolished.txt"))
      #Protein IDs
  ids <- read_tsv(paste0(AgProtect,"/AllProteinIds-AgProtect.txt"), col_names = "ID")

  #Polish IDs
  ids$AdhesinProbability <- "-"
  ids$AdhesinProbability2 <- "Unknown"  
  
  #Polish results
  raw <- select(raw, "Protein name (Annotation)", "Pad-value")
  raw$`Protein name (Annotation)` <- gsub("\\.1.*", ".1", raw$`Protein name (Annotation)`)
  colnames(raw)[1] = "ID"
  colnames(raw)[2] <- "AdhesinProbability"
  raw$ID <- raw$ID %>% str_replace(">", "")
  raw$AdhesinProbability <- as.numeric(raw$AdhesinProbability)
  raw$AdhesinProbability2 = ifelse(raw$AdhesinProbability > 0.7 ,
                                   "Adhesin",
                                   "Non-Adhesin")
  raw$AdhesinProbability <- as.character(raw$AdhesinProbability)
  
  #Merge
  ids <- ids %>% filter(!ID %in% raw$ID)
  raw <- left_join(ids, raw, by=NULL)
  
  # Write the final results to a TSV file
  write_tsv(raw, paste0(AgProtect,"/Final_results/SPAAN-AgProtect-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging SPAAN results: ", conditionMessage(e)))
})
