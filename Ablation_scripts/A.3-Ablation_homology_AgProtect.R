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
library("rlang")

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)  # Capture arguments
if (length(args) != 4) {
  stop("Path to AgProtect argument is required")
}

AgProtect <- args[1]
evalue <- args[2]
bits<- args[3]
ide <- args[4]

evalue <- as.numeric(evalue)

tryCatch({
  ##### Read all the files
  #IDs
  ids <- read_tsv(paste0(AgProtect,"/AllProteinIds-AgProtect.txt"), 
                  col_names = "ID")
  totales<- nrow(ids)

  #blastp results for different hosts
  Hom_files = list.files(paste0(AgProtect,"/Homology_Analysis_results/"),
                         full.names = F)
  #List of antigens and hosts
  HostList <- read_tsv(paste0(AgProtect,"/HostList.tsv"), col_names = TRUE)
  HostList  <- HostList %>% 
    mutate(Host = str_replace_all(Host, " ", "")) %>%
    separate_rows(Host, sep = ";")
  
  ######First, we process raw results for all species
  i = 1
  
  Temp <- tibble("ID" = character(),"ide" = double(),
                 "evalue" = double(),"bits" = double(),
                 "Host"= character())
  
  for (i in 1:length(Hom_files)) {
    
    raw <- read_tsv(paste0(AgProtect,"/Homology_Analysis_results/", Hom_files[i]),
                    col_names = F) 
    
    # Define the criteria for filtering the results
    index = ifelse(((raw$X11 < evalue) + (raw$X12>bits) + (raw$X3 > ide)) == 3, TRUE, FALSE)
    
    # Filter the raw data based on the defined criteria
    out = raw[index, c(1, 3, 11, 12)]
    out = out %>% as_tibble() %>%
      distinct(X1, .keep_all = TRUE)
    
    colnames(out) <- c("ID","ide", "evalue","bits")
    
    out$Host <- as.character(Hom_files[i])
    
    Temp <- full_join(out, Temp)
    
  }
   
  ######Then, we filter only the hosts that are relevant for that antigen
  
  Temp <- Temp %>%
    mutate(Host = str_remove(Host, "\\.out$"))
  
  todas <- left_join(HostList, Temp, by=c("ID", "Host") )
  
  todas <- na.omit(todas)
  
  todas <- todas %>% distinct(ID, .keep_all = TRUE)
  
  pasaron <- nrow(todas)
  
  porc <- round(pasaron/totales, 4)
  
  cat(ide, evalue, bits, porc, "\n")
  
  rm(list = ls())  
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging homology analysis results: ", conditionMessage(e)))
})
