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

cat("Filtering, cleaning, and sorting VaxiJen results\n")

tryCatch({
  #Read files
      #VaxiJen results
  vaxijen_data <- as_tibble(read_table(paste0(MAIN, "/", subdir, "/VaxiJen_results/VaxiJen3_predictions.csv"),
                             col_names = c("ID", "X2", "X3", "X4", "X5", "AntigenicityResult", "X7", "X8", "X9")))
  #Rare aminoacids
  RareAA <- read_tsv(paste0(MAIN, "/", subdir, "/VaxiJen_results/IDs_RareAminoAcids"), col_names = "ID")
  RareAA$AntigenicityResult <- "-"
  
  # Extracting IDs
ids <- vaxijen_data %>%
  filter(str_detect(ID, "WP_|HBL79")) %>%
  select(ID)

# Extracting result
predicted <- vaxijen_data %>%
  filter(str_detect(AntigenicityResult, "ANTIGEN")) %>%
  select(AntigenicityResult)

# Pasting IDs and refined data
Results <- bind_cols(ids, predicted)
Results <- bind_rows(Results, RareAA)

# Arrange the combined data frame by the first column. Keep only distinct proteins
Results <- arrange(Results, ID) %>%
  distinct(ID, .keep_all = TRUE)
Results$ID <- Results$ID %>% str_replace("\\.1", "")

write_tsv(Results, paste0(MAIN, "/",  subdir,"/Final_results/VaxiJen-", subdir, "-final.tsv"))

rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging VaxiJen results: ", conditionMessage(e)))
})
