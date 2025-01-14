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
home <- "/path/to/directory"

cat("Filtering, cleaning, and sorting the results of SPAAN\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw = read_tsv(paste0(home, "/WorkDir/", subdir, "/SPAAN_results/SPAAN-unpolished.txt"))
      #Protein IDs
  ids <- read_tsv(paste0(home, "/WorkDir/", subdir, "/AllProteinIds-",subdir,".txt"), col_names = "ID")
  colnames(ids)[1] <- "ID"
  
  #Polish IDs
  ids$ID <- ids$ID %>% str_replace("\\.1", "")
  ids$AdhesinProbability <- "-"
  ids$AdhesinProbability2 <- "Unknown"  
 
   #Polish results
  raw <- select(raw, "Protein name (Annotation)", "Pad-value")
  raw$`Protein name (Annotation)` <- gsub("\\.1.*", ".1", raw$`Protein name (Annotation)`)
  colnames(raw)[1] = "ID"
  colnames(raw)[2] <- "AdhesinProbability"
  raw$ID <- raw$ID %>% str_replace("\\.1", "")
  raw$ID <- raw$ID %>% str_replace(">", "")
  raw$AdhesinProbability <- as.numeric(raw$AdhesinProbability)
  raw$AdhesinProbability2 = ifelse(raw$AdhesinProbability > 0.7 ,
                                   "Adhesin",
                                   "Non-Adhesin")
  raw$AdhesinProbability <- as.character(raw$AdhesinProbability)
  
  #Merge
  raw <- full_join(raw, ids, by=NULL)
  
  
  # Write the final results to a TSV file
  write_tsv(raw, paste0(home, "/WorkDir/", subdir, "/Final_results/SPAAN-", subdir, "-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging SPAAN results: ", conditionMessage(e)))
})
