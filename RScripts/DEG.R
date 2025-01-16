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
cat("Filtering, cleaning, and sorting DEG results\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw <- read_tsv(paste0(MAIN, "/", subdir, "/DEG_results/", subdir, "-vs_DEG.out"), col_names = 1:12)
      #Protein IDs
  ids <- read_tsv(paste0(MAIN, "/", subdir, "/AllProteinIds-",subdir,".txt"), col_names = "ID")

  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < 1e-5) + (raw$X12>50) + (raw$X3 > 25)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble()
  
  #Polish data
  colnames(out)[1:5] = c("ID","Essential.Homologue","Essential.IdentityPercent",
                         "Essential.Evalue","Essential.Bitscore")
                        
  # Create data frame with proteins that were not classified as essential
  MissingIDs <- ids[!(ids$ID %in% out$ID), ]
  df <- data.frame(matrix("-", nrow = nrow(MissingIDs), ncol = 4))
  colnames(df) <- c("Essential.Homologue","Essential.IdentityPercent",
                    "Essential.Evalue","Essential.Bitscore")
  MissingIDs <- cbind(MissingIDs, df)
  
  # Combine the raw results with the missing proteins
  AllResults <- rbind(out, MissingIDs)
  
  # Arrange the combined data frame. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID, Essential.Bitscore)%>%
    distinct(ID, .keep_all = TRUE)
  
  #Polish final results
  AllResults$ID <- AllResults$ID %>% str_replace("\\.1", "")
  
  # Create essential protein column
  AllResults$EssentialProtein = ifelse(AllResults$Essential.Homologue == "-",
                                       "Non Essential",
                                       "Essential")
  
  
  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(MAIN, "/", subdir,"/Final_results/DEG-", subdir, "-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging DEG results: ", conditionMessage(e)))
})
