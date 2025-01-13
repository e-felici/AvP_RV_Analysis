#!/usr/bin/Rscript

#load library
library(tidyverse)

# Get the command line arguments
args <- commandArgs(TRUE)

# Get the folder argument
subdir <- args[1]

#Here do not forget to define the home directory
home <- "/path/to/directory"

# Display message to indicate the start of the process
cat("Filtering, cleaning, and sorting PSORTb results\n")

# Use tryCatch to handle errors during the process
tryCatch({
  # Get the list of files matching the specific pattern in the "PSORTb_results" directory
  lines <- list.files(paste0(home, "/WorkDir/", subdir,"/PSORTb_results/"), 
                      pattern = "_psortb_gramneg.txt", 
                      full.names = F)
  
  # Ensure exactly one file is found; otherwise, handle the error
  if (length(lines) == 1) {
    # Read the raw results
    tibbler <- read.delim(paste0(home, "/WorkDir/", subdir,"/PSORTb_results/",lines[1]),
                          header = TRUE, sep = "\t")
    
    #Polish results
    # Keep only the first string up to the first space
    tibbler <- tibbler %>%
      mutate(SeqID = sub(" .*", "", SeqID))
    tibbler <- select(tibbler, SeqID, Localization)
    colnames(tibbler)[1:2] <- c("ID", "SubcellularLocalization")
    
    # Create Exposition column
    tibbler$Exposition = as.factor(ifelse(tibbler$SubcellularLocalization == "Extracellular" | 
                                            tibbler$SubcellularLocalization == "OuterMembrane",
                                             "Exposed",
                                             ifelse(tibbler$SubcellularLocalization == "Unknown",
                                                    "Unknown", 
                                                    "Inside the cell")))
    tibbler$ID <- tibbler$ID %>% str_replace("\\.1", "")
    tibbler <- arrange(tibbler, ID)
    
    # Write the cleaned and processed data to a new TSV file
    write_tsv(tibbler, paste0(home, "/WorkDir/", subdir, "/Final_results/PSORTb-", subdir, "-final.tsv"))
    rm(list = ls())
    
  } else {
    # Stop the process and show an error if no or multiple files match the pattern
    stop("Error: expected exactly one file matching pattern '_psortb_gramneg.txt'")
  }
  
}, error = function(e) {
  # Handle any unexpected errors and print the error message
  cat("An error occurred:", e$message, "\n")
})
