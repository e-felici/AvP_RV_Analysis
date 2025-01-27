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
if (length(args) != 1) {
  stop("Path to AgProtect argument is required")
}

AgProtect <- args[1]

# Print a message indicating the start of filtering, cleaning, and sorting the homology analysis results
cat("Filtering, cleaning, and sorting the results of the homology analysis from AgProtect")


# Attempt to read the TSV file containing the raw homology analysis results
tryCatch({
  ##### Read all the files
  #IDs
  ids <- read_tsv(paste0(AgProtect,"/AllProteinIds-AgProtect.txt"), 
                  col_names = "ID")
  #blastp results for different hosts
  Hom_files = list.files(paste0(AgProtect,"/Homology_Analysis_results/"),
                          full.names = F)
  #List of antigens and hosts
  HostList <- read_tsv(paste0(AgProtect,"/HostList.tsv"), col_names = TRUE)
  HostList  <- HostList %>% mutate(Host = str_replace_all(Host, " ", ""))
 
  ######First, we process raw results for all species
  #Initialize variable i and vector conservation
  i = 1
  AllResults <- tibble("ID" = character(),
                   "Host" = character(),
                   "Host.Homologue" = character(),
                   "Host.IdentityPercent" = double(),
                   "Host.Evalue" = double(),
                   "Host.Bitscore" = double())
  
  #Entropy analysis for all files
  for (i in 1:length(Hom_files)) {
    
    raw <- read_tsv(paste0(AgProtect,"/Homology_Analysis_results/", Hom_files[i]),
                     col_names = F) 
    
  # Define the criteria for filtering the results
  index = ifelse(((raw$X11 < 1e-5) + (raw$X12>50) + (raw$X3 > 25)) == 3, TRUE, FALSE)
  
  # Filter the raw data based on the defined criteria
  out = raw[index, c(1, 2, 3, 11, 12)]
  out = out %>% as_tibble()
  
  #Polish data
  colnames(out)[1:5] = c("ID", "Host.Homologue","Host.IdentityPercent",
                         "Host.Evalue","Host.Bitscore")
  out$Host <- as.character(Hom_files[i])
  
  # Create data frame with proteins that were not classified as host homologues
  MissingIDs <- ids[!(ids$ID %in% out$ID), ]
  df <- data.frame(matrix("-", nrow = nrow(MissingIDs), ncol = 4))
  df$Host <- as.character(Hom_files[i])
  colnames(df) <- c("Host.Homologue","Host.IdentityPercent",
                    "Host.Evalue","Host.Bitscore", "Host")
  MissingIDs <- cbind(MissingIDs, df)
  
  # Combine the raw results and the missing IDs
  Temp <- rbind(out, MissingIDs)
  
  AllResults <- full_join(AllResults, Temp, by = c("ID", "Host", "Host.Homologue",
                                                   "Host.IdentityPercent", "Host.Evalue", 
                                                   "Host.Bitscore"))
  }
  # Arrange the combined data frame by the first column. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID) %>%
    distinct(ID, .keep_all = TRUE)

  #
  write_tsv(AllResults, paste0(AgProtect,"/Homology_Analysis_results/AgProtect-Homology-all.tsv"))
  
  ######Then, we filter only the hosts that are relevant for that antigen
   #Split the Host column into separate rows
  HostList2 <- HostList %>%
    separate_rows(Host, sep = ";")
  
  AllResults$Host_Homologue_Result <- ifelse(AllResults$Host.Homologue == "-",
                                             "Non Host_Homologue",
                                             "Host Homologue")
  
  AllResults <- inner_join(AllResults,HostList2, by= join_by(ID,Host))
  
  rm(HostList2)
  
  AllResults <- AllResults %>%
    select(ID, Host, Host_Homologue_Result) %>%
    pivot_wider(names_from = Host, values_from = Host_Homologue_Result) %>%
    mutate(across(everything(), ~ replace_na(.x, "-")))
  
  # Transform the tibble
  AllResults <- AllResults %>%
    rowwise() %>%
    mutate(
      Host_Homologue_Result = {
        # Get all non-missing valid values (exclude "-") except ID column
        valid_values <- c_across(-ID) %>%  # Exclude ID from c_across
          keep(~ .x != "-")
        
        # Determine the result based on the unique valid values
        unique_values <- unique(valid_values)
        
        case_when(
          all(unique_values == "Host Homologue") ~ "Host Homologue",
          all(unique_values == "Non Host Homologue") ~ "Non Host Homologue",
          length(unique_values) > 1 ~ "Host Homologue and Non Host Homologue, depending on the strain",
          TRUE ~ NA_character_
        )
      }
    ) %>%
    ungroup()
  
  AllResults <- AllResults %>%
    select(ID,Host_Homologue_Result) %>%
    mutate(
      Host.Homologue = "See:AgProtect-Homology-all.tsv",
      Host.IdentityPercent = "See:AgProtect-Homology-all.tsv",
      Host.Evalue = "See:AgProtect-Homology-all.tsv",
      Host.Bitscore= "See:AgProtect-Homology-all.tsv"
    )
  
  AllResults <- inner_join(HostList, AllResults, by = "ID")
  
  rm(HostList)
  
  # Arrange the combined data frame by the first column. Keep only distinct proteins
  AllResults <- arrange(AllResults, ID, Host.Bitscore) %>%
    distinct(ID, .keep_all = TRUE)
  
  
  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/AgProtect-Homology-final.tsv"))
  
  rm(list = ls())
  
  }, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging homology analysis results: ", conditionMessage(e)))
})
