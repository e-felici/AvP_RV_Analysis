#!/usr/bin/Rscript

#Catch Balcony tar path
args <- commandArgs(trailingOnly = TRUE)  
if (length(args) != 1) {
  stop("Please, do not forget to define the path to BALCONY_0.2.10.tar.gz!")
}
BALCONY_path <- args[1]

#install packages
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
if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}
if (!requireNamespace("rlang", quietly = TRUE)) {
  install.packages("rlang")
}
if (!requireNamespace("Rpdb", quietly = TRUE)) {
  install.packages("Rpdb")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!requireNamespace("BALCONY", quietly = TRUE)) {
  install.packages(BALCONY_path, repos = NULL, type="source")
}
cat("Everything installed!")

