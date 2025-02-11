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
  i = 1
  AllResults <- tibble("ID" = character(),
                   "Host" = character(),
                   "Host.Homologue" = character(),
                   "Host.IdentityPercent" = double(),
                   "Host.Evalue" = double(),
                   "Host.Bitscore" = double())
  
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
  Temp$Host.IdentityPercent <- as.double(Temp$Host.IdentityPercent)
  Temp$Host.Bitscore <- as.double(Temp$Host.Bitscore)
  Temp$Host.Evalue <- as.double(Temp$Host.Evalue)
  
  AllResults <- full_join(AllResults, Temp, by = c("ID", "Host", "Host.Homologue",
                                                   "Host.IdentityPercent", "Host.Evalue", 
                                                   "Host.Bitscore"))
  }
  
  
  write_tsv(AllResults, paste0(AgProtect,"/Homology_Analysis_results/AgProtect-Homology-all.tsv"))
  
  ######Then, we filter only the hosts that are relevant for that antigen
   #Split the Host column into separate rows
  HostList2 <- HostList %>%
    separate_rows(Host, sep = ";")
  
  AllResults$Host_Homologue_Result <- ifelse(AllResults$Host.Homologue == "-",
                                             "Non_Host_Homologue",
                                             "Host_Homologue")
  AllResults <- AllResults %>%
    mutate(Host = str_remove(Host, "\\.out$"))
  
  # Keep distinct strings in 'ID' that maximize Bit Score
  AllResults <- AllResults%>%
    group_by(ID, Host) %>% 
    arrange(desc(Host.Bitscore), desc(Host.IdentityPercent)) %>% # Sort by BitScore and then idpercent (both descending)
    slice(1) %>%                                # Keep the first row in each group
    ungroup()       
  
  AllResults <- inner_join(AllResults,HostList2, by= join_by(ID,Host))
  
  rm(HostList2)
  
  data <- AllResults %>%
    group_by(ID) %>%
    mutate(
      Host_Homologue_Result_All = case_when(
        all(Host_Homologue_Result == "Host_Homologue") ~ "Host Homologue",
        all(Host_Homologue_Result == "Non_Host_Homologue") ~ "Non Host Homologue",
        TRUE ~ "Host Homologue and Non Host Homologue, depending on the Host"
      )
    ) %>%
    ungroup()
  
  compressed_data <- data %>%
    group_by(ID) %>%
    summarize(
      Host = paste(Host, collapse = ";"),
      Host.Homologue = paste(Host.Homologue, collapse = ";"),
      Host.IdentityPercent = paste(Host.IdentityPercent, collapse = ";"),
      Host.Evalue = paste(Host.Evalue, collapse = ";"),
      Host.Bitscore = paste(Host.Bitscore, collapse = ";"),
      .groups = "drop"
    )
  
  data <- data%>%
    select(ID, Host_Homologue_Result_All) %>%
    distinct(ID, .keep_all = TRUE)
  
AllResults <- full_join(compressed_data,data, by="ID")

AllResults  <- AllResults %>% mutate(Host.IdentityPercent = str_replace_all(Host.IdentityPercent, "NA", "-"))
AllResults  <- AllResults %>% mutate(Host.Evalue = str_replace_all(Host.Evalue, "NA", "-"))
AllResults  <- AllResults %>% mutate(Host.Bitscore = str_replace_all(Host.Bitscore, "NA", "-"))

  # Write the final results to a TSV file
  write_tsv(AllResults, paste0(AgProtect,"/Final_results/AgProtect-Homology-final.tsv"))
  
  rm(list = ls())  
  }, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging homology analysis results: ", conditionMessage(e)))
})
