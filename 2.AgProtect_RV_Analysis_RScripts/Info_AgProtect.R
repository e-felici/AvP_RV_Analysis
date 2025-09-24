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
if (length(args) != 2) {
  stop("Missing required arguments!")
}

AgProtect <- args[1]
WorkDir <- args[2]

tryCatch({
  #Read files and extract info
  subfolders <- list.dirs(WorkDir, full.names = TRUE, 
                          recursive = FALSE)
  
  # Count the number of subfolders
  num_subfolders <- length(subfolders)
  
  Info <- read_lines(paste0(AgProtect,"/protein.faa")) %>%
    keep(~ str_starts(.x, ">")) %>%
    tibble(Full_Header = .)  %>%
    mutate(
      Full_Header = str_replace_all(Full_Header, ">", ""),     # Remove ">"
      Full_Header = str_replace_all(Full_Header, " ", "_"),      # Replace spaces with underscores (if necessary)
      ID = str_extract(Full_Header, "^[^_]+_[^_]+"),             # Extract the IDs (includes one underscore)
      Description = str_remove(Full_Header, "^[^_]+_[^_]+_?"),
      Strain = "Experimental_Antigens",
      Mean_Schneider_Entropy = "0", 
      Strain_count = num_subfolders, 
      Cluster_Number = "-1",
      Conservation = "1",
      Conservation_Results = "-"
    )
  
  # Arrange the combined data frame by the first column. Keep only distinct
  Info <- arrange(Info, ID) %>%
    distinct(ID, .keep_all = TRUE) %>%
    select(-Full_Header)

  # Write the final results to a TSV file
  write_tsv(Info, paste0(AgProtect,"/Final_results/Info-AgProtect-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging Info results: ", conditionMessage(e)))
})
