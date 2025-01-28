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

cat("Filtering, cleaning, and sorting DeepTMHMM results\n")

tryCatch({
  # Step 1: Combine all predicted_topologies* files into one
  type_files <- list.files(paste0(AgProtect,"/DeepTMHMM_results"),
                           pattern = "^predicted_topologies", full.names = TRUE)
  
  #Initialize things
  protType <- tibble(ID = character(), Type_of_Protein = character())
  
  for (file in type_files) {
    #read files
    temp_protType <- as_tibble(read_delim(file,
                                          col_names = FALSE, delim = "|"))
    
    #Polish tibble
    temp_protType$X1 = str_replace_all(temp_protType$X1, " ", "")
    temp_protType$X2 = str_replace_all(temp_protType$X2, " ", "")
    temp_protType <- temp_protType %>%
      filter(str_detect(X1, ">"))
    temp_protType$X1 <- temp_protType$X1 %>% str_replace(">", "")
    colnames(temp_protType)[1:2] <- c("ID", "Type_of_Protein")
    
    #merge with others
    protType <- full_join(protType, temp_protType, by = c("ID", "Type_of_Protein"))
    
  }  
  
  
  
  # Step 2: Copy and process TMRs* files
  TMR_files <- list.files(paste0(AgProtect,"/DeepTMHMM_results"), 
                          pattern = "^TMRs", full.names = TRUE)
  
  #Initialize things
  TMbreakdown <- tibble(ID= character(),
                        Inside_cell_or_cytosol= double(),
                        Periplasm= double(), 
                        Outside_cell= double(), 
                        Beta_sheet_TM= double(), 
                        Alpha_helix_TM= double(), 
                        SignalPeptideNumber= double())
  required_cols <- c("ID","inside","periplasm","outside",
                     "Beta","TMhelix", "signal")
  
  for (File in TMR_files) {
    #read file
    temp_TMbreakdown <- as_tibble(read.table(File,
                                             sep = "", 
                                             comment.char = "#", 
                                             fill = TRUE, 
                                             stringsAsFactors = FALSE,
                                             row.names=NULL,
                                             col.names=c("ID", "Location", "X","Y","Z")))
    #Polish
    temp_TMbreakdown <- temp_TMbreakdown %>% 
      filter(ID != "//") %>%
      select(ID, Location) %>%
      filter(!if_all(everything(), ~ is.na(.) | . == "")) %>%
      filter(str_detect(ID, "AgProtect"))
    
    temp_TMbreakdown <- temp_TMbreakdown  %>%
      group_by(ID, Location) %>%
      summarise(Count = n(), .groups = "drop") %>%
      as_tibble()
    
    temp_TMbreakdown <- temp_TMbreakdown  %>% pivot_wider(
      names_from = Location,
      values_from = Count 
    ) %>%
      mutate(across(everything(), ~ replace_na(.x, 0)))
    
    # Check for missing columns and add them if necessary
    missing_cols <- setdiff(required_cols, colnames(temp_TMbreakdown))
    
    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        temp_TMbreakdown <- temp_TMbreakdown %>%
          mutate(!!col := 0) 
      }
    }
    
    temp_TMbreakdown <- temp_TMbreakdown %>%
      rename(
        Inside_cell_or_cytosol = inside,
        Periplasm = periplasm,
        Outside_cell = outside, 
        Beta_sheet_TM = Beta,
        Alpha_helix_TM = TMhelix, 
        SignalPeptideNumber = signal
      )
    TMbreakdown <- full_join(TMbreakdown, temp_TMbreakdown, by = NULL)
    
  }  
  
  # Combine all files
  AllResults <- left_join(protType, TMbreakdown, by = "ID")
  
  # Arrange the combined data frame by the first column
  AllResults <- arrange(AllResults, ID) %>%
    distinct(ID, .keep_all = TRUE)

  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/DeepTMHMM-AgProtect-final.tsv"))
  
  rm(list = ls())
  
  }, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging DeppTMHMM results: ", conditionMessage(e)))
})
