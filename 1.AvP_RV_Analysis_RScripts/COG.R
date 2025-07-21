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

cat("Filtering, cleaning, and sorting COG results\n")
tryCatch({
  ####### Read all necessary files####### 
  ##Raw results
  lines <- readLines(paste0(MAIN, "/", subdir, "/COG_results/", subdir, "-multi_protein.out"))
  
  ##Protein Ids
  AllProteinIds = read_tsv(paste0(MAIN, "/", subdir, "/AllProteinIds-", subdir, ".txt"),
                           col_names = "ID")
  AllProteinIds$ID <- AllProteinIds$ID %>% str_replace("\\.1", "")
  
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
  to_remove <- grep("^#", lines)
  
  lines_clean <- lines[-to_remove]
  
  #clean
  rm(to_remove)
  
  # Find indices of QUERY and DOMAINS lines
  query_idx <- grep("^QUERY", lines_clean)
  domain_idx <- grep("^DOMAINS", lines_clean)
  
  query_ids <- as.tibble(lines_clean[query_idx])
  domain_9 <- tibble(lines_clean[domain_idx + 1])
  
  split_matrix_query <- str_split_fixed(query_ids$value, "\\s+", 12)
  split_matrix_dom <- str_split_fixed(domain_9$`lines_clean[domain_idx + 1]`, "\\s+", 12)
  
  query_ids <- as_tibble(split_matrix_query, .name_repair = "minimal")
  domain_9 <- as_tibble(split_matrix_dom, .name_repair = "minimal")
  
  colnames(query_ids) <- paste0("V", 1:12)
  colnames(domain_9) <- paste0("V", 1:12)
  
  query_ids <- query_ids %>%
    distinct(V2, .keep_all = TRUE)
  domain_9 <- domain_9 %>%
    distinct(V2, .keep_all = TRUE)
  
  query_ids <- query_ids %>% select(V2, V5)
  domain_9 <- domain_9 %>% select(V2, V9)
  
  result  <- full_join(query_ids, domain_9)
  
  rm(list = ls(pattern = "^domain"))
  rm(list = ls(pattern = "^query"))
  rm(list = ls(pattern = "^split"))
  rm(list = ls(pattern = "^lines"))
  
  result <- result %>% select(-V2)
  colnames(result) <- c("ID", "COG_ID")
  

  
  ##From COG IDs obtain COG Functional Categories 
  result <- left_join(result, cogdef, by = "COG_ID")
  result <- result %>% select(ID, COG_ID, COG_categories)
  
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
  
  colnames(result_superfamily)[2] <- "cl_ID"
  colnames(result_superfamily)[7] <- "COG_ID"
  
  result_superfamily <- select(result_superfamily, ID, COG_ID)
  
  result_superfamily <- left_join(result_superfamily, cogdef, by = "COG_ID")
  
  result_superfamily <- result_superfamily %>%
    group_by(ID) %>%
    summarise(
      COG_ID = str_c(COG_ID, collapse = ", "),
      COG_categories = first(COG_categories)
    )
  
  #clean  
  rm(cogdef)
  
  #On the other hand, some COG IDs are old, so here we use the manually curated file:    
  result_COG_old <- inner_join(result_superfamily, COG_olds, by ="COG_ID")
  
  rm(COG_olds)
  
  result_COG_old <- result_COG_old %>% select(-COG_categories.x, -PSSM_ID)
  
  colnames(result_COG_old)[3] <- "COG_categories"
  
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
