#!/usr/bin/Rscript

library(tidyverse)

# Get the command line arguments
args <- commandArgs(TRUE)

# Get the folder argument
subdir <- args[1]

#Here do not forget to define the home directory
home <- "/path/to/directory/"

tryCatch({
  #Read files and extract info
  file_list <- list.files(paste0(home, "/WorkDir/", subdir,"/Final_results/"), 
                      full.names = T)

  # Read all files and assign tibbles named after the files (without extensions)
  file_tibbles <- file_list %>%
    set_names(~ basename(.x) %>% tools::file_path_sans_ext()) %>% # Use file names as names
    map(~ read_tsv(.x, col_names = TRUE))                        # Read each file into a tibble
  
  # Join all tibbles by the "ID" column
  ALL <- file_tibbles %>%
    reduce(full_join, by = "ID")  # Use "by = 'ID'" to specify the common column

# Write the final results to a TSV file
write_tsv(ALL, paste0(home, "/WorkDir/", subdir, "/Final_Results-", subdir, ".tsv"))
rm(list = ls())

}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging final polishing of this strain results: ", conditionMessage(e)))
})
