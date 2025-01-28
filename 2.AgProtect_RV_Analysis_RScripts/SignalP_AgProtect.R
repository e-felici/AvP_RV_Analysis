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

cat("Filtering, cleaning, and sorting signalP results\n")

tryCatch({
  #Read raw results
  lines <- as.tibble(read_tsv(paste0(AgProtect,"/SignalP_results/prediction_results.txt"), skip = 1))
  
  #Polish results
  colnames(lines)[1:2] = c("ID", "SignalPeptide")
  lines = select(lines, "ID", "SignalPeptide")
  lines <- lines %>% mutate(ID = sub(" .*", "", ID))
  
  # Write the final results to a TSV file
  write_tsv(lines, paste0(AgProtect,"/Final_results/SignalP-AgProtect-final.tsv"))
  rm(list = ls())

  }, error = function(e) {
  warning(paste0("Error arranging SignalP results: ", conditionMessage(e)))
})

