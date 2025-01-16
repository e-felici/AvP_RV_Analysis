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
if (length(args) != 2) {
  stop("Two arguments are required")
}

AgProtect <- args[1]
COG <- args[2]

cat("Filtering, cleaning, and sorting COG results\n")

tryCatch({
#Read all necessary files
    ##Raw rps-blast results
AllResults = read_tsv(paste0(AgProtect, "/COG_Results/AgProtect-vs_COG.out"), 
                      col_names = 1:12)

    ##Protein Ids
AllProteinIds = read_tsv(paste0(AgProtect, "/AllProteinIds-AgProtect.txt"),
                         col_names = "ID")

    ##COG and CDD resources:
      #Summary info about the CD models
cddid <- as.tibble(read_tsv(paste0(COG, "/cddid.tbl"), 
                            col_names = F))
      #Polish tibble
colnames(cddid)[1] = "CDDID"
cddid$CDDID <- as.character(cddid$CDDID)

      #Info about COG descriptions
cogdef <- as.tibble(read_tsv(paste0(COG,"/cog-20.def.tab"), 
                             col_names = F))
      #Polish tibble
cogdef <- select(cogdef, X1, X2)
colnames(cogdef)[1] = "COGID"

      #Descriptions of COG functional categories
cogdet <- as.tibble(read_tsv(paste0(COG,"/fun-20.tab"), 
                             col_names = F))
      #Polish tibble
colnames(cogdet)[1] = "COG"


##From Best-hit results obtain CDD ID
#Filter rows and subset data 
index = ifelse(((AllResults$X11 < 1e-5) + (AllResults$X12 > 50) + (AllResults$X3 > 25)) == 3, TRUE, FALSE)
AllResults = AllResults[index, c(1, 2, 3, 11, 12)]

# Polish data
AllResults = AllResults %>% 
  as_tibble() %>% 
  group_by(X1) %>% 
  filter(row_number() == 1)
colnames(AllResults)[1] = "ID"

# Merge data with IDs to add the proteins that do not fit into any of the COG categories
AllResults <- full_join(AllResults, AllProteinIds, by = "ID")

# Polish data
AllResults$X2 <- gsub("CDD:", "", AllResults$X2)
AllResults$X2 <- as.character(AllResults$X2)
AllResults$X3 <- as.character(AllResults$X3)
AllResults$X11 <- as.character(AllResults$X11)
AllResults$X12 <- as.character(AllResults$X12)
colnames(AllResults)[2] = "CDDID"

##From CDD IDs obtain COG IDs
# Merge
AllResults <- left_join(AllResults, cddid, by = "CDDID")

# Polish data
AllResults <- select(AllResults, ID, CDDID, X3.y, X4, X3.x, X11, X12, X2)
colnames(AllResults)[8] = "COGID"

##From COG IDs obtain COG Functional Categories 
# Merge
AllResults <- left_join(AllResults, cogdef, by = "COGID")

# Extract only the first functional category
AllResults <- AllResults %>% mutate(COG = substr(X2, 1, 1))

##Add little explanation about COG first functional categories
# Merge
AllResults <- left_join(AllResults, cogdet, by = "COG")

#Polish data 
replace_spaces <- function(column) {
  column %>% str_replace_all(" ", "_")
}
AllResults <- AllResults %>%
  mutate(across(everything(), replace_spaces))

colnames(AllResults)[1:12] <- c("ID", "PSSM-Id", "CD_short_name", "CD_description", 
                            "COG_IdentityPercent", "COG.Evalue", "COG.Bitscore", "COG",
                            "COG_category", "COG_category_u", "COG_color", 
                            "COG_category_description")


# Arrange the combined data frame by the first column. Keep only distinct
AllResults <- arrange(AllResults, ID) %>%
  distinct(ID, .keep_all = TRUE)

# Write the final results to a TSV file
write_tsv(AllResults, paste0(AgProtect,"/Final_results/COG-AgProtect-final.tsv"),
          na = "-")
rm(list = ls())

}, 
# Catch any errors that may occur during execution
error = function(e) {
  # Print a warning message if an error occurs
  warning(paste0("Error arranging COG results: ", conditionMessage(e)))
})
