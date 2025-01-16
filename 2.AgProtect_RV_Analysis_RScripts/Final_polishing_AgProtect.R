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

tryCatch({
  #Read files and extract info
  file_list <- list.files(paste0(AgProtect,"/Final_results/"), 
                      full.names = T)

  # Read all files and assign tibbles named after the files (without extensions)
  file_tibbles <- file_list %>%
    set_names(~ basename(.x) %>% tools::file_path_sans_ext()) %>% # Use file names as names
    map(~ read_tsv(.x, col_names = TRUE))                        # Read each file into a tibble
  
  # Join all tibbles by the "ID" column
  ALL <- file_tibbles %>%
    reduce(full_join, by = "ID")  

# Write the final results to a TSV file
write_tsv(ALL, paste0(AgProtect,"/Final_Results-AgProtect.tsv"))
rm(list = ls())

}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging final polishing of AgProtect results: ", conditionMessage(e)))
})
