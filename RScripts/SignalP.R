#!/usr/bin/Rscript


if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("purrr", quietly = TRUE)) {
  install.packages("purrr")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble")
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  install.packages("magrittr")
}
if (!requireNamespace("forcats", quietly = TRUE)) {
  install.packages("forcats")
}


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
args <- commandArgs(TRUE)

# Get the folder argument
subdir <- args[1]

#Here do not forget to define the home directory
home <- "/path/to/directory/"

cat("Filtering, cleaning, and sorting the results of the homology analysis\n")

tryCatch({
  #Read raw results
  lines <- as_tibble(read_tsv(paste0(home, "/WorkDir/", subdir, "/SignalP_results/prediction_results.txt"), skip = 1))
  
  #Polish results
  colnames(lines)[1:2] = c("ID", "SignalPeptide")
  lines = select(lines, "ID", "SignalPeptide")
  lines$ID <- lines$ID %>% str_replace("\\.1", "")
  lines <- lines %>% mutate(ID = sub(" .*", "", ID))
  
  # Write the final results to a TSV file
  write_tsv(lines, paste0(home, "/WorkDir/", subdir, "/Final_results/SignalP-", subdir, "-final.tsv"))
  rm(list = ls())
  
  }, error = function(e) {
  warning(paste0("Error arranging SignalP results: ", conditionMessage(e)))
})

