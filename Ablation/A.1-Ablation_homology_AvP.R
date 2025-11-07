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
if (length(args) != 5) {
  stop("Two arguments are required!")
}

MAIN <- args[1]
subdir <- args[2]
evalue <- args[3]
bits<- args[4]
ide <- args[5]

evalue <- as.numeric(evalue)

tryCatch({
  #Read all necessary files:
      #Raw results
  raw = read_tsv(paste0(MAIN,"/", subdir, "/Homology_Analysis_results/", subdir, "-vs_chicken.out"), col_names = 1:12)
      #Protein IDs
  ids <- read_tsv(paste0(MAIN,"/", subdir, "/AllProteinIds-",subdir,".txt"), col_names = "ID")
  colnames(ids)[1] <- "ID"
  
  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < evalue) + (raw$X12>bits) + (raw$X3 > ide)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble() %>%
    distinct(X1, .keep_all = TRUE)
  
  pasaron <- nrow(out)
  
  ids <- read_tsv(paste0(MAIN,"/", subdir, "/AllProteinIds-",subdir,".txt"), col_names = "ID")
  totales <- nrow(ids)
  
  porc <- round(pasaron/totales, 4)
  
  cat(ide, evalue, bits, porc, "\n")
  
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error : ", conditionMessage(e)))
})
