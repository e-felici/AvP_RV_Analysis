#!/usr/bin/Rscript

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
if (length(args) != 2) {
  stop("Two arguments are required!")
}

FinalRes <- args[1]
AgProtect <- args[2]

tryCatch({
  # Read files 
  #All AvP strains
  #Read files and extract info
  file_list <- list.files(path = FinalRes, 
                          pattern = "Final_results-", 
                          full.names = T)
  # Read all files and assign tibbles named after the files (without extensions)
  file_tibbles <- file_list %>%
    set_names(~ basename(.x) %>% tools::file_path_sans_ext()) %>% # Use file names as names
    map(~ read_tsv(.x, col_names = TRUE))                        # Read each file into a tibble
  
  # Join all tibbles by the "ID" column
  Results <- file_tibbles %>%
    reduce(full_join) 
  
  Conservation <- read_tsv(paste0(FinalRes,"/Conservation_results.tsv"), col_names = TRUE)  
  
  Results <- inner_join(Results, Conservation, by = "ID")
  
  #Adding Experimental Antigens
  AgProtect <- read_tsv(paste0(AgProtect,"/Final_Results-AgProtect.tsv"), col_names = TRUE)
  
  Results <- full_join(Results, AgProtect, by = NULL)
  
  Results$Type_of_Protein2 = ifelse(Results$Type_of_Protein == "BETA",
                                    "Beta Barrel",
                                    ifelse(Results$SignalPeptide == "LIPO",
                                           "Lipoprotein",
                                           ifelse(Results$Type_of_Protein == "TM",
                                                  ifelse(Results$Alpha_helix_TM == 1,
                                                         "1 Alpha helix Transmembrane Step",
                                                         "2 or more Alpha helix Transmembrane Steps"),
                                                  "Secreted")))
  
  write_tsv(Results, paste0(FinalRes, "/AllStrains_AgProtect_Final_results.tsv"))
  
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging definitive results: ", conditionMessage(e)))
})
