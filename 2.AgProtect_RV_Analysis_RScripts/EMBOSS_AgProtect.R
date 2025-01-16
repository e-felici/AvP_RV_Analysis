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
cat("Filtering, cleaning, and sorting EMBOSS results\n")

tryCatch({
  # Read raw results 
  lines <- readLines(paste0(AgProtect,"/EMBOSS_results/EMBOSS-", subdir))
  
  #Extract and polish the IDs
  IDs <- lines[grep("PEPSTATS", lines)] %>% strsplit((" ")[[1]][3]) %>% unlist() %>% as_tibble()
  IDs$value <- IDs$value %>% str_replace("^PEPSTATS of ", "")
  IDs$value <- IDs$value %>% str_replace("from 1 to .*", "")
  IDs$value <- str_replace_all(IDs$value, "\\s+", "")
  colnames(IDs)[1] = "ID"
  
  #Extract and polish info about Molecular weight
  PMs <- lines[grep("Molecula", lines)] %>% strsplit((" ")[[1]][3]) %>% unlist() %>% as_tibble()  
  PMs$value <- PMs$value %>% str_replace("^Molecular weight = ", "")
  PMs$value <- PMs$value %>% str_replace("\tResidues = .*", "") 
  PMs$value <- str_replace_all(PMs$value, "\\s+", "")
  PMs$value <- as.numeric(PMs$value)
  colnames(PMs)[1] = "MW"
  
  #Extract and polish info about isoelectric point
  PIs <- lines[grep("Isoelec", lines)] %>% strsplit((" ")[[1]][3]) %>% unlist() %>% as_tibble()
  PIs$value <- PIs$value %>% str_replace("^Isoelectric Point = ", "")
  PIs$value <- as.numeric(PIs$value)
  colnames(PIs)[1] = "Isoelectric_Point"
  
  #Extract and polish info about probability of expression in inclusion bodies
  InclBods <- lines[grep("bability of", lines, ignore.case = F)] %>% strsplit((" ")[[1]][3]) %>% unlist() %>% as_tibble()
  InclBods$value2 = ifelse(InclBods$value == "^Probability of *",
                           InclBods$value %>% str_replace("Probability of expression in inclusion bodies = ", "") %>% as.numeric(),
                           InclBods$value %>% str_replace("^Improbability of expression in inclusion bodies = ", "") %>% as.numeric())
  InclBods$value2 = ifelse(is.na(InclBods$value2) == T,
                           InclBods$value,
                           1-(InclBods$value2))
  InclBods$value2 <- str_replace_all(InclBods$value2, "\\s+", "")
  InclBods$value2 <- InclBods$value2 %>% str_replace("Probabilityofexpressionininclusionbodies=", "") %>% as.numeric()
  InclBods <- select(InclBods, value2)                      
  colnames(InclBods)[1] = "Inclussion_Bodies_Probability"
  
  #Merge all data
  EMBOSS = cbind(IDs,PMs,PIs,InclBods)
  
  EMBOSS$Isoelectric_Point2 = ifelse(EMBOSS$Isoelectric_Point > 8,
                                     "> 8",
                                     "< 8")
  
  EMBOSS$MW2 = ifelse(EMBOSS$MW > 50000,
                      ifelse(EMBOSS$MW > 100000,
                             "> 100 kDa",
                             "50-100 kDa"),
                      "< 50 kDa")
  
  #Arrange results
  EMBOSS <- arrange(EMBOSS, ID)
  
  #Write the final results to TSV
  write_tsv(EMBOSS, paste0(AgProtect,"/Final_results/EMBOSS-AgProtect-final.tsv"))
  rm(list = ls())
}, error = function(e) {
  warning(paste0("Error arranging results from EMBOSS: ", conditionMessage(e)))
})
