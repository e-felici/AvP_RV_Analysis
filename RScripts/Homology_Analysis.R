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

cat("Filtering, cleaning, and sorting the results of the homology analysis\n")

tryCatch({
  #Read all necessary files:
      #Raw results
  raw = read_tsv(paste0(home, "/WorkDir/", subdir, "/Homology_Analysis_results/", subdir, "-vs_chicken.out"), col_names = 1:12)
      #Protein IDs
  ids <- read_tsv(paste0(home, "/WorkDir/", subdir, "/AllProteinIds-",subdir,".txt"), col_names = "ID")
  colnames(ids)[1] <- "ID"
  
  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < 1e-5) + (raw$X12>50) + (raw$X3 > 25)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble()
  
  #Polish data
  colnames(out)[1:5] = c("ID", "Host.Homologue","Host.IdentityPercent",
                         "Host.Evalue","Host.Bitscore")
  out$Host <- "Gallus gallus"
  
  # Create data frame with proteins that were not classified as host homologues
  MissingIDs <- ids[!(ids$ID %in% out$ID), ]
  df <- data.frame(matrix("-", nrow = nrow(MissingIDs), ncol = 4))
  df$Host <- "Gallus gallus"
  colnames(df) <- c("Host.Homologue","Host.IdentityPercent",
                    "Host.Evalue","Host.Bitscore", "Host")
  MissingIDs <- cbind(MissingIDs, df)
  
  # Combine the raw results and the missing IDs
  AllResults <- rbind(out, MissingIDs)
  
  # Arrange the combined data frame. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID, Host.Bitscore)
  AllResults <- AllResults %>%
    distinct(ID, .keep_all = TRUE)
  
  #Polish final results
  AllResults$ID <- AllResults$ID %>% str_replace("\\.1", "")
  
  # Create host homologues column
  AllResults$Host_Homologue_Result = ifelse(AllResults$Host.Homologue == "-",
                                       "Non Host Homologue",
                                       "Host Homologue")
  
  
  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(home, "/WorkDir/", subdir, "/Final_results/HostHomology-", subdir, "-final.tsv"))
  rm(list = ls())
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging homology analysis results: ", conditionMessage(e)))
})
