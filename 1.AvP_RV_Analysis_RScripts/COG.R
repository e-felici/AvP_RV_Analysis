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
  stop("Three arguments are required!")
}

MAIN <- args[1]
subdir <- args[2]
COG <- args[3]

MAIN <-"~/Busqueda_antigenos/ALL"
subdir <- "GCF_000442905.1"
COG <- "~/rpsbproc"

cat("Filtering, cleaning, and sorting COG results\n")
tryCatch({
  ####### Read all necessary files####### 
  ##Raw results
  lines <- readLines(paste0(MAIN, "/", subdir, "/COG_results/", subdir, "-multi_protein.out"))
  
  ##Protein Ids
  AllProteinIds = read_tsv(paste0(MAIN, "/", subdir, "/AllProteinIds-", subdir, ".txt"),
                           col_names = "ID")
  AllProteinIds$ID <- AllProteinIds$ID %>% str_replace("\\.1", "")
  
  ##Summary info about the CD models
  cddid <- as.tibble(read_tsv(paste0(COG, "/data/cddid.tbl"), 
                              col_names = F))
  colnames(cddid)[1:2] = c("PSSM_ID","COG_ID")
  cddid$PSSM_ID <- as.character(cddid$PSSM_ID)
  
  ## COG definitions: Info about COG descriptions
  cogdef <- as.tibble(read_tsv(paste0(COG, "/data/cog-24.def.tab"), 
                               col_names = F))
  cogdef <- select(cogdef, X1, X2)
  colnames(cogdef)[1:2] = c("COG_ID", "COG_categories")
  
  ## Descriptions of COG functional categories
  cogdet <- as.tibble(read_tsv(paste0(COG, "/data/cog-24.fun.tab"), 
                               col_names = F, comment = "#"))
  colnames(cogdet)[1:3] = c("COG", "General_Functional_Category", "COG_category_description")
  
  ##Superfamily info: Descriptions of COG functional categories
  superfamily <- as.tibble(read_tsv(paste0(COG, "/data/family_superfamily_links"), 
                                    col_names = F))
  colnames(superfamily)[1] <- "target"
  colnames(superfamily)[3] <- "COG_ID"
  superfamily <- superfamily %>%
    filter(if_any(everything(), ~ str_detect(as.character(.), "COG")))
  superfamily <- superfamily %>% select(target, COG_ID)
  
  ##Manually curated COGs (old)
  COG_olds <- as.tibble(read_tsv(paste0(COG, "/data/cog_anot_manual.tbl"), 
                                 col_names = T))
  COG_olds$PSSM_ID <-  as.character(COG_olds$PSSM_ID) 
  
  
  
  
  ###### Raw results Processing ######
  # Find all SITES and ENDSITES indices
  sites_starts <- grep("^SITES", lines)
  sites_ends <- grep("^ENDSITES", lines)
  
  # Build a vector of all lines to remove
  to_remove <- unlist(mapply(function(start, end) start:end, sites_starts, sites_ends))
  to_remove2 <- grep("#", lines)
  to_remove3 <- grep("DATA", lines)
  to_remove4 <- grep("SESSION", lines)
  
  to_remove <- c(to_remove, to_remove2,to_remove3,to_remove4)
  
  # Keep only lines not in to_remove
  lines_clean <- lines[-to_remove]
  
  #clean
  rm(list = ls(pattern = "^to_remove"))
  rm(list = ls(pattern = "^sites"))
  
  # Find indices of QUERY and DOMAINS lines
  query_idx <- grep("^QUERY", lines_clean)
  domain_idx <- grep("^DOMAINS", lines_clean)
  
  # Extract 5th column from QUERY lines: IDs
  query_ids <- sapply(strsplit(lines_clean[query_idx], "\\s+"), `[`, 5)
  
  # For each DOMAINS, extract the 9th column from the next line (the domain data line): PSSM_IDs
  domain_9th <- sapply(domain_idx, function(i) {
    domain_line <- lines_clean[i + 1]
    strsplit(domain_line, "\t")[[1]][4]
  })
  
  # Combine into a tibble
  result <- tibble(
    ID = query_ids,
    PSSM_ID= domain_9th
  )
  
  result <- inner_join(result, cddid, by = "PSSM_ID")
  
  #clean
  rm(cddid)
  rm(list = ls(pattern = "^domain"))
  rm(list = ls(pattern = "^query"))
  rm(list = ls(pattern = "^lines"))
  
  ##From COG IDs obtain COG Functional Categories 
  result <- left_join(result, cogdef, by = "COG_ID")
  result <- result %>% select(ID, PSSM_ID, COG_ID, COG_categories)
  
  # Extract only the first functional category
  result <- result %>% mutate(COG = substr(COG_categories, 1, 1))
  
  ##Add little explanation about COG first functional categories
  result <- left_join(result, cogdet, by = "COG")
  
  #Some PSSM IDs are superfamily cluster records (i.e: start with "cl", check: 
  # https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CDSequenceCluster)
  #And as a consequence, they could not be processed with previous code (NA rows).
  #So, to "translate" the to COG categories:
  result_superfamily <- result %>% 
    filter(if_any(everything(), is.na))
  
  result_superfamily <- inner_join(result_superfamily, superfamily, by = "COG_ID")
  
  rm(superfamily)
  
  colnames(result_superfamily)[3] <- "cl_ID"
  colnames(result_superfamily)[8] <- "COG_ID"
  
  result_superfamily <- select(result_superfamily, ID, PSSM_ID, COG_ID)
  
  result_superfamily <- left_join(result_superfamily, cogdef, by = "COG_ID")
  
  result_superfamily <- result_superfamily %>%
    group_by(ID) %>%
    summarise(
      COG_ID = str_c(COG_ID, collapse = ", "),
      PSSM_ID = first(PSSM_ID),
      COG_categories = first(COG_categories)
    )
  
  #clean  
  rm(cogdef)
  
  #On the other hand, some COG IDs are old, so here we use the manually curated file:    
  result_COG_old <- inner_join(result_superfamily, COG_olds, by =c("COG_ID", "PSSM_ID"))
  
  rm(COG_olds)
  
  result_COG_old <- result_COG_old %>% select(-COG_categories.x)
  
  colnames(result_COG_old)[4] <- "COG_categories"
  
  result_superfamily <- result_superfamily %>% drop_na()
  
  result_superfamily <- result_superfamily %>% mutate(COG = substr(COG_categories, 1, 1))
  
  result_superfamily <- left_join(result_superfamily, cogdet, by = "COG")
  
  rm(cogdet)
  
  result <- full_join(result, result_COG_old)
  result <- full_join(result, result_superfamily)
  
  rm(result_COG_old)
  rm(result_superfamily)
  
  result$ID <- result$ID %>% str_replace("\\.1", "")
  
  # Merge data with IDs to add the proteins that do not fit into any of the COG categories
  result <- full_join(result, AllProteinIds, by = "ID")
  
  #Polish data 
  replace_spaces <- function(column) {
    column %>% str_replace_all(" ", "_")
  }
  result <- result %>%
    mutate(across(everything(), replace_spaces))
  
  result <- result %>%
    mutate(across(everything(), ~replace_na(., "-")))
  
  # Arrange the combined data frame by the first column. Keep only distinct
  result <- arrange(result, ID) %>%
    distinct(ID, .keep_all = TRUE)
  
  # Write the final results to a TSV file
  write_tsv(result, paste0(MAIN, "/", subdir, "/Final_results/COG-", subdir, "-final.tsv"),
            na = "-")
  
  rm(list = ls())
  
}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging COG results: ", conditionMessage(e)))
})
