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
if (length(args) != 3) {
  stop("Three arguments are required")
}

AgProtect <- args[1]
DBDir <- args[2]
VF_cat <- args[3]

cat("Filtering, cleaning, and sorting VFDB results\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw <- read_tsv(paste0(AgProtect,"/VFDB_full_results/AgProtect-vs_VFDB.out"), col_names = 1:12)
      #Protein IDs
  ids <- read_tsv(paste0(AgProtect,"/AllProteinIds-AgProtect.txt"), col_names = "ID")

  vf_cat <- read_tsv(VF_cat)
  
  vfdb <- read_lines(paste0(DBDir, "/VFDB_full/VFDB_setB_pro.faa"))
  
  keep <- grep("^>", vfdb)
  
  vfdb <- vfdb[keep]
  
  # Extract up to 3 matches per line
  vf_matches <- str_match_all(vfdb, "\\(?VF\\w+")
  
  # Pad each match to length 3 (with NA if fewer found)
  vf_matrix <- t(sapply(vf_matches, function(x) { length(x) <- 3; x }))
  
  # Build tibble
  vfdb <- as_tibble(vf_matrix, .name_repair = "minimal") %>%
    setNames(paste0("col", 1:3))
  
  rm(vf_matches)
  rm(vf_matrix)
  
  vfdb <- vfdb %>%
    mutate(across(everything(), ~ str_replace_all(., "\\(", "")))
  
  colnames(vfdb) <- c("VF_ID", "delete", "VFCID")
  
  vfdb <- select(vfdb, -delete)
  
  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < 1e-5) + (raw$X12>50) + (raw$X3 > 25)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble()
  
  # Polish data
  colnames(out)[1:5] = c("ID", "Virulence.Homologue","Virulence.IdentityPercent",
                         "Virulence.Evalue","Virulence.Bitscore")
  
  # Create data frame with proteins that were not classified as VF
  MissingIDs <- ids[!(ids$ID %in% out$ID), ]
  df <- data.frame(matrix("-", nrow = nrow(MissingIDs), ncol = 4))
  colnames(df) <- c("Virulence.Homologue","Virulence.IdentityPercent",
                    "Virulence.Evalue","Virulence.Bitscore")
  MissingIDs <- cbind(MissingIDs, df)
  
  # Combine the results with the missing proteins
  AllResults <- rbind(out, MissingIDs)
  
  # Arrange the combined data frame. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID, Virulence.Bitscore) %>%
    distinct(ID, .keep_all = TRUE)
  
  # Create virulence column
  AllResults$VirulenceFactor = ifelse(AllResults$Virulence.Homologue == "-",
                                      "Non Virulence Factor",
                                      "Probable Virulence Factor")
  
  AllResults$VF_ID <- str_replace(AllResults$Virulence.Homologue, "\\(.*", "")
  
  AllResults <- left_join(AllResults, vfdb, by="VF_ID")
  
  AllResults <- left_join(AllResults, vf_cat, by = "VFCID")
  
  AllResults <-  AllResults %>%
    mutate(across(everything(), ~replace_na(., "-")))
  
  AllResults <- arrange(AllResults, ID, Virulence.Bitscore) %>%
    distinct(ID, .keep_all = TRUE)
  
    # Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/AgProtect-vs_VFDB_final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging VFDB results: ", conditionMessage(e)))
})
