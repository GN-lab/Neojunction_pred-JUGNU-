#!/usr/bin/env Rscript
# Step 14a: Generate input CSV file for running MHCFlurry 2.0
# Modified to read HLA alleles directly from OptiType results per sample
# November 24, 2025 | Gaurav Raichand | The Institute of Cancer Research
#     The input CSV file is expected to contain columns 
#     "allele", "peptide", and, optionally, "n_flank", and "c_flank".
#
# To run on command line: $ mhcflurry-predict INPUT.csv –out RESULT.csv
# Prerequisite: Source config.sh to set INPUT_DIR and OPTITYPE_OUTPUT_DIR

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(data.table)  # For faster reading and operations

# Load Directories from Env (sourced from config.sh)
input_dir    <- Sys.getenv("INPUT_DIR")          # /data/rds/DMP/UCEC/EVOLIMMU/graichand/Neojuction_pred/SSNIP/0_Input_Files
directory_12 <- Sys.getenv("STEP12_OUTPUT_DIR")  # == OUTPUT_DIR from config
directory_14 <- Sys.getenv("STEP14_OUTPUT_DIR")  # == OUTPUT_DIR from config

# Path to OptiType results (configure in config.sh or use default)
optitype_dir <- Sys.getenv("OPTITYPE_OUTPUT_DIR")

# Input filenames (n-mers from Step 12 outputs in directory_12)
nmers_08_file <- "2023_0812_hlathenalist_msic_08mers.tsv"
nmers_09_file <- "2023_0812_hlathenalist_msic_09mers.tsv"
nmers_10_file <- "2023_0812_hlathenalist_msic_10mers.tsv"
nmers_11_file <- "2023_0812_hlathenalist_msic_11mers.tsv"

# Samples file in INPUT_DIR
samples_file <- file.path(input_dir, "samples.txt")  # Adjust if samples.txt is elsewhere

# Hardcode run date for outputs
run_date <- "2025_1124"

# Load Samples List
user_samples <- fread(
  samples_file,
  header   = FALSE,
  col.names = "sample",
  fill     = TRUE,
  na.strings = c("", "NA")
)$sample

user_samples <- unique(user_samples[!is.na(user_samples)])
cat("Loaded", length(user_samples), "user samples from", samples_file, "\n")

###########################################################################
# Extract HLA Alleles from OptiType Results -------------------------
###########################################################################

extract_hla_from_optitype <- function(optitype_dir, user_samples) {
  # Extract HLA alleles from OptiType *_result.tsv files per sample
  # Returns a data.table with columns: sample, alleles_str (pipe-separated, HLA format)
  # OptiType output has columns: A1, A2, B1, B2, C1, C2, Reads, Objective
  
  hla_results_list <- list()
  
  for (sample in user_samples) {
    result_file <- file.path(optitype_dir, sample, paste0(sample, "_result.tsv"))
    
    if (!file.exists(result_file)) {
      cat("  WARNING: OptiType result file not found for", sample, "\n")
      next
    }
    
    # Read OptiType result
    optitype_result <- tryCatch(
      fread(result_file, header = TRUE, na.strings = c("", "NA")),
      error = function(e) {
        cat("  ERROR reading", result_file, ":", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(optitype_result) || nrow(optitype_result) == 0) {
      cat("  WARNING: No data in OptiType result for", sample, "\n")
      next
    }
    
    alleles <- c()
    col_names <- tolower(colnames(optitype_result))
    
    # helper to normalize OptiType allele (e.g. "A*24:02") to "HLA-A*24:02"
    norm_allele <- function(raw, locus) {
      if (is.na(raw) || raw == "") return(NA_character_)
      raw <- trimws(raw)
      # remove leading "A*", "B*", "C*" from OptiType output
      core <- sub("^[ABC]\\*", "", raw)
      paste0("HLA-", locus, "*", core)
    }
    
    # Extract A alleles (A1, A2)
    if ("a1" %in% col_names && "a2" %in% col_names) {
      a1_raw <- as.character(optitype_result[[which(col_names == "a1")[1]]][1])
      a2_raw <- as.character(optitype_result[[which(col_names == "a2")[1]]][1])
      a1 <- norm_allele(a1_raw, "A")
      a2 <- norm_allele(a2_raw, "A")
      if (!is.na(a1) && a1 != "") alleles <- c(alleles, a1)
      if (!is.na(a2) && a2 != "") alleles <- c(alleles, a2)
    }
    
    # Extract B alleles (B1, B2)
    if ("b1" %in% col_names && "b2" %in% col_names) {
      b1_raw <- as.character(optitype_result[[which(col_names == "b1")[1]]][1])
      b2_raw <- as.character(optitype_result[[which(col_names == "b2")[1]]][1])
      b1 <- norm_allele(b1_raw, "B")
      b2 <- norm_allele(b2_raw, "B")
      if (!is.na(b1) && b1 != "") alleles <- c(alleles, b1)
      if (!is.na(b2) && b2 != "") alleles <- c(alleles, b2)
    }
    
    # Extract C alleles (C1, C2)
    if ("c1" %in% col_names && "c2" %in% col_names) {
      c1_raw <- as.character(optitype_result[[which(col_names == "c1")[1]]][1])
      c2_raw <- as.character(optitype_result[[which(col_names == "c2")[1]]][1])
      c1 <- norm_allele(c1_raw, "C")
      c2 <- norm_allele(c2_raw, "C")
      if (!is.na(c1) && c1 != "") alleles <- c(alleles, c1)
      if (!is.na(c2) && c2 != "") alleles <- c(alleles, c2)
    }
    
    if (length(alleles) > 0) {
      alleles_unique <- unique(alleles)
      alleles_str <- paste(alleles_unique, collapse = "|")
      hla_results_list[[sample]] <- data.table(
        alleles_str = alleles_str,
        sample      = sample
      )
      cat("  Extracted", length(alleles_unique), "unique alleles for", sample, "\n")
    } else {
      cat("  WARNING: No valid alleles extracted for", sample, "\n")
    }
  }
  
  if (length(hla_results_list) == 0) {
    stop("No HLA alleles successfully extracted from OptiType results.")
  }
  
  hla_data_combined <- rbindlist(hla_results_list)
  return(hla_data_combined)
}

# Extract HLA alleles from OptiType per sample
cat("\nExtracting HLA alleles from OptiType results directory:\n", optitype_dir, "\n\n")
hla_data_filtered <- extract_hla_from_optitype(optitype_dir, user_samples)
cat("\nSuccessfully extracted HLA data for", nrow(hla_data_filtered), "samples.\n\n")

###########################################################################
#  ORIGINAL: Parse Alleles ------------------------------------------------
###########################################################################

# Parse alleles: Split by "|", extract unique class I (HLA-A/B/C)
hla_alleles_list <- strsplit(hla_data_filtered$alleles_str, "\\|")
valid_lists      <- hla_alleles_list[lengths(hla_alleles_list) > 0]
all_alleles      <- unique(unlist(valid_lists))
class_i_alleles  <- all_alleles[grepl("^HLA-[ABC]", all_alleles)]

alleles_tib <- tibble(allele = class_i_alleles)
cat("Extracted", nrow(alleles_tib), "unique class I HLA alleles from OptiType results.\n")
cat("Alleles:", paste(class_i_alleles, collapse = ", "), "\n\n")

# Load Nmer Files (from Step 12 in directory_12)
setwd(directory_12)
nmers_08 <- fread(nmers_08_file, na.strings = c("", "NA"))
nmers_09 <- fread(nmers_09_file, na.strings = c("", "NA"))
nmers_10 <- fread(nmers_10_file, na.strings = c("", "NA"))
nmers_11 <- fread(nmers_11_file, na.strings = c("", "NA"))

###########################################################################
#  Step 1: Edit Dataframes ------------------------------------------------
###########################################################################

# Remove TPM if present
if ("TPM" %in% colnames(nmers_08)) nmers_08[, TPM := NULL]
if ("TPM" %in% colnames(nmers_09)) nmers_09[, TPM := NULL]
if ("TPM" %in% colnames(nmers_10)) nmers_10[, TPM := NULL]
if ("TPM" %in% colnames(nmers_11)) nmers_11[, TPM := NULL]

# Rename columns for MHCflurry (n_mer -> peptide, ctex_up -> n_flank, ctex_dn -> c_flank)
setnames(nmers_08, old = colnames(nmers_08), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_09, old = colnames(nmers_09), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_10, old = colnames(nmers_10), new = c("peptide", "n_flank", "c_flank"))
setnames(nmers_11, old = colnames(nmers_11), new = c("peptide", "n_flank", "c_flank"))

# Clean flanks: Remove "-" and ensure character
nmers_08[, `:=`(
  n_flank = gsub("-", "", as.character(n_flank)),
  c_flank = gsub("-", "", as.character(c_flank))
)]
nmers_09[, `:=`(
  n_flank = gsub("-", "", as.character(n_flank)),
  c_flank = gsub("-", "", as.character(c_flank))
)]
nmers_10[, `:=`(
  n_flank = gsub("-", "", as.character(n_flank)),
  c_flank = gsub("-", "", as.character(c_flank))
)]
nmers_11[, `:=`(
  n_flank = gsub("-", "", as.character(n_flank)),
  c_flank = gsub("-", "", as.character(c_flank))
)]

# Replicate peptides for unique alleles (efficient data.table way)
alleles_dt <- as.data.table(alleles_tib)

nmer_dfs     <- list(nmers_08, nmers_09, nmers_10, nmers_11)
nmer_lengths <- c(8, 9, 10, 11)

for (i in 1:4) {
  nmer_i <- nmer_dfs[[i]]
  if (nrow(nmer_i) == 0 || nrow(alleles_dt) == 0) {
    cat("Skipping", nmer_lengths[i], "mers: empty data.\n")
    next
  }

  nmer_expanded    <- nmer_i[rep(1:.N, times = nrow(alleles_dt))]
  alleles_repeated <- alleles_dt[rep(1:.N, each = nrow(nmer_i))]
  nmer_all         <- cbind(alleles_repeated, nmer_expanded)

  if ("peptide.1" %in% names(nmer_all)) {
    setnames(nmer_all, "peptide.1", "peptide")
  }

  assign(paste0(nmer_lengths[i], "mer_final"), nmer_all)
  cat("Processed", nmer_lengths[i], "mers:", nrow(nmer_all), "rows (peptides × alleles).\n")
}

# Write to STEP14_OUTPUT_DIR
setwd(directory_14)
suffix <- paste0(run_date, ".csv")  # Distinguish from original

fwrite(get("8mer_final"),
       paste0("08mer_mhcflurry_input_", suffix),
       sep = ",", na = "NA", col.names = TRUE, quote = TRUE)

fwrite(get("9mer_final"),
       paste0("09mer_mhcflurry_input_", suffix),
       sep = ",", na = "NA", col.names = TRUE, quote = TRUE)

fwrite(get("10mer_final"),
       paste0("10mer_mhcflurry_input_", suffix),
       sep = ",", na = "NA", col.names = TRUE, quote = TRUE)

fwrite(get("11mer_final"),
       paste0("11mer_mhcflurry_input_", suffix),
       sep = ",", na = "NA", col.names = TRUE, quote = TRUE)

print("Step 14a: MHCflurry input CSVs generated with alleles extracted from OptiType results per sample.")
