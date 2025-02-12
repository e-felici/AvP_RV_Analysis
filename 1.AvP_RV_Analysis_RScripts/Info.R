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
  stop("Two arguments are required!")
}

MAIN <- args[1]
subdir <- args[2]

tryCatch({
  #Read files and extract info
  Info <- read_lines(paste0(MAIN, "/", subdir, "/protein.faa")) %>%
        keep(~ str_starts(.x, ">")) %>%
        tibble(Full_Header = .)  %>%
    mutate(
      Full_Header = str_replace_all(Full_Header, ">", ""),     # Remove ">"
      Full_Header = str_replace_all(Full_Header, " ", "_"),      # Replace spaces with underscores (if necessary)
      ID = str_extract(Full_Header, "^[^_]+_[^_]+"),             # Extract the IDs (includes one underscore)
      Description = str_remove(Full_Header, "^[^_]+_[^_]+_?"),
      Strain = subdir,
      Genus = "Avibacterium",
      Gram_Stain = "Negative",
      Bacteria = "Avibacterium_paragallinarum"
    )
  
  # Arrange the combined data frame by the first column. Keep only distinct
  Info <- arrange(Info, ID) %>%
    distinct(ID, .keep_all = TRUE) %>%
    select(ID, Description, Strain, Genus, Gram_Stain)
  Info$ID <- Info$ID %>% str_replace("\\.1", "")
  
  # Write the final results to a TSV file
  write_tsv(Info, paste0(MAIN, "/", subdir, "/Final_results/Info-", subdir, "-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging Info results: ", conditionMessage(e)))
})
