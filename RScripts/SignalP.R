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

cat("Filtering, cleaning, and sorting the results of the homology analysis\n")

tryCatch({
  #Read raw results
  lines <- as_tibble(read_tsv(paste0(MAIN, "/", subdir, "/SignalP_results/prediction_results.txt"), skip = 1))
  
  #Polish results
  colnames(lines)[1:2] = c("ID", "SignalPeptide")
  lines = select(lines, "ID", "SignalPeptide")
  lines$ID <- lines$ID %>% str_replace("\\.1", "")
  lines <- lines %>% mutate(ID = sub(" .*", "", ID))
  
  # Write the final results to a TSV file
  write_tsv(lines, paste0(MAIN, "/", subdir, "/Final_results/SignalP-", subdir, "-final.tsv"))
  rm(list = ls())
  
  }, error = function(e) {
  warning(paste0("Error arranging SignalP results: ", conditionMessage(e)))
})

