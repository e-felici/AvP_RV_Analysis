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
  stop("Path to AgProtect is required")
}

AgProtect <- args[1]

cat("Filtering, cleaning, and sorting COG results\n")

tryCatch({
#Read all files
  Loc_files = list.files(paste0(AgProtect,"/PSORTb_results/"),
                         full.names = F)
  Bact_type <- read_tsv(paste0(AgProtect,"/BacteriaTypes.tsv"), 
                        col_names = TRUE)
  i = 1
  AllResults <- tibble("ID" = character(),
                       "SubcellularLocalization" = character(),
                       "Gram_Stain" = character())
  
  for (i in 1:length(Loc_files)) {
   Temp <- read_tsv(paste0(AgProtect,"/PSORTb_results/", Loc_files[i]),
                  col_names = c("ID","SubcellularLocalization","X")) 
    
    Temp <- Temp %>%
      select(-X) %>%
      mutate(ID = sub(" .*", "", ID))
    
    Temp$Gram_Stain <- Loc_files[i]
    
    AllResults <- full_join(AllResults, Temp, by = c("ID","SubcellularLocalization", "Gram_Stain"))
  }
  
  Bact_type2 <- Bact_type %>% select(-Genus,-Bacteria)
  
  AllResults <- inner_join(AllResults,Bact_type2, by= join_by(ID, Gram_Stain))
  
  AllResults <- AllResults  %>%
    unique() %>% 
    pivot_wider(names_from = Gram_Stain, values_from = SubcellularLocalization) %>%
    mutate(across(everything(), ~ replace_na(.x, "-")))
  
  # Create a single column with the non-"-" values
  AllResults <- AllResults %>%
    rowwise() %>%
    mutate(all_values = list(c_across(-ID)[c_across(-ID) != "-"])) %>% 
    mutate(all_values = if_else(length(all_values) == 1, all_values[[1]], NA_character_)) %>%
    ungroup()
  
  AllResults <- select(AllResults, ID, all_values)
  colnames(AllResults)[2] <- "SubcellularLocalization"
  
  AllResults <- inner_join(AllResults,Bact_type, by = "ID")
  
  # Create Exposition column
  AllResults$Exposition = as.factor(ifelse(AllResults$SubcellularLocalization == "Extracellular" | 
                                          AllResults$SubcellularLocalization == "OuterMembrane",
                                        "Exposed",
                                        ifelse(AllResults$SubcellularLocalization == "Unknown",
                                               "Unknown", 
                                               "Inside the cell")))
  
  AllResults <- arrange(AllResults, ID)

# Write the final merged data frame to a TSV file
write_tsv(AllResults, paste0(AgProtect,"/Final_results/PSORTb-AgProtect-final.tsv"))
rm(list = ls())

}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging COG results: ", conditionMessage(e)))
})
