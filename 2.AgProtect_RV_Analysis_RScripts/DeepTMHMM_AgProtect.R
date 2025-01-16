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

cat("Filtering, cleaning, and sorting DeepTMHMM results\n")

tryCatch({
# Read results
      #Protein type:
  protType = as.tibble(read_delim(
    paste0(AgProtect,"/DeepTMHMM_Results/DeepTMHMM_type_AgProtect.txt"),
    col_names = F, delim = "|"))
  #Polish tibble
  protType$X1 = str_replace_all(protType$X1, " ", "")
  protType$X2 = str_replace_all(protType$X2, " ", "")
  protType <- protType %>%
    filter(str_detect(X1, ">"))
  protType$X1 <- protType$X1 %>% str_replace(">", "")
  
  #Breakdown of TM topology
  TMbreakdown = as.tibble(read_tsv(
    paste0(AgProtect,"/DeepTMHMM_Results/DeepTMHMM_TMbreakdown_AgProtect.txt"), 
    col_names = F))

 
  # Combine all files
  AllResults <- left_join(protType, TMbreakdown, by = "X1")

#Polish results
colnames(AllResults)[1:9] <- c("ID","Type_of_Protein","Inside_cell_or_cytosol",
                              "Periplasm", "Outside_cell", "Beta_sheet_TM", "Alpha_helix_TM", 
                              "SignalPeptideNumber")

# Arrange the combined data frame by the first column
  AllResults <- arrange(AllResults, ID) %>%
    distinct(ID, .keep_all = TRUE)

# Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/DeepTMHMM-AgProtect-final.tsv"))
  
  rm(list = ls())
  
  }, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging DeppTMHMM results: ", conditionMessage(e)))
})
