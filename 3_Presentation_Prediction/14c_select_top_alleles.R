#!/usr/bin/env Rscript
# Step 14c: Select TOP HLA-allele presentation score for each n-mer (MHCFlurry 2.0)
# Now processes BOTH ALT and WT/self sides, mirroring the same logic.
# Updated to read HLA per sample directly from OptiType results
# July 03, 2026 | Gaurav Raichand | The Institute of Cancer Research
# Prerequisite: Source config.sh to set INPUT_DIR, OUTPUT_DIR, STEP14_OUTPUT_DIR, OPTITYPE_OUTPUT_DIR

###########################################################################
# Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(data.table)
library(parallel)

# Load Directories from environment
input_dir      <- Sys.getenv("INPUT_DIR")
directory_14_in  <- Sys.getenv("OUTPUT_DIR")       # Where MHCflurry outputs are (Step 14b)
directory_14_out <- Sys.getenv("STEP14_OUTPUT_DIR")
optitype_dir     <- Sys.getenv("OPTITYPE_OUTPUT_DIR")

# Auto-detect the newest matching file instead of a hardcoded date (NEW).
# The old hardcoded RUN_DATE <- "2025_1124" no longer matched what 14b.sh
# actually writes (RUN_DATE="2026_1124" there) -- meaning this script was
# silently reading stale November 2025 predictions instead of the current
# run's output, with no error to warn you. This can't go stale the same
# way again, since it always picks whatever file is actually newest.
latest_file <- function(pattern, dir = ".") {
  files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files match pattern ", pattern, " in ", dir)
  files[which.max(file.info(files)$mtime)]
}

###########################################################################
# NEW: Load samples and extract HLA from OptiType -------------------------
###########################################################################

# Load sample IDs
samples_list <- fread(
  file.path(input_dir, "samples.txt"),
  header = FALSE, col.names = "sample"
)$sample %>%
  unique()

cat("Loaded", length(samples_list), "samples from samples.txt\n")

# Function copied / adapted from Step 14a
extract_hla_from_optitype <- function(optitype_dir, user_samples) {
  # Returns data.table with columns: sample, alleles_str (pipe-separated HLA-A/B/C)
  hla_results_list <- list()
  
  for (sample in user_samples) {
    result_file <- file.path(optitype_dir, sample, paste0(sample, "_result.tsv"))
    
    if (!file.exists(result_file)) {
      cat("  WARNING: OptiType result file not found for", sample, "\n")
      next
    }
    
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
      core <- sub("^[ABC]\\*", "", raw)
      paste0("HLA-", locus, "*", core)
    }
    
    # A1, A2
    if ("a1" %in% col_names && "a2" %in% col_names) {
      a1_raw <- as.character(optitype_result[[which(col_names == "a1")[1]]][1])
      a2_raw <- as.character(optitype_result[[which(col_names == "a2")[1]]][1])
      a1 <- norm_allele(a1_raw, "A")
      a2 <- norm_allele(a2_raw, "A")
      if (!is.na(a1) && a1 != "") alleles <- c(alleles, a1)
      if (!is.na(a2) && a2 != "") alleles <- c(alleles, a2)
    }
    
    # B1, B2
    if ("b1" %in% col_names && "b2" %in% col_names) {
      b1_raw <- as.character(optitype_result[[which(col_names == "b1")[1]]][1])
      b2_raw <- as.character(optitype_result[[which(col_names == "b2")[1]]][1])
      b1 <- norm_allele(b1_raw, "B")
      b2 <- norm_allele(b2_raw, "B")
      if (!is.na(b1) && b1 != "") alleles <- c(alleles, b1)
      if (!is.na(b2) && b2 != "") alleles <- c(alleles, b2)
    }
    
    # C1, C2
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
        sample      = sample,
        alleles_str = alleles_str
      )
      cat("  Extracted", length(alleles_unique), "unique alleles for", sample, "\n")
    } else {
      cat("  WARNING: No valid alleles extracted for", sample, "\n")
    }
  }
  
  if (length(hla_results_list) == 0) {
    stop("No HLA alleles successfully extracted from OptiType results.")
  }
  
  rbindlist(hla_results_list)
}

cat("\nExtracting HLA alleles from OptiType results directory:\n", optitype_dir, "\n\n")
hla_data_filtered <- extract_hla_from_optitype(optitype_dir, samples_list)
cat("\nSuccessfully extracted HLA data for", nrow(hla_data_filtered), "samples (rows in hla_data_filtered).\n\n")

# Expand pipe-separated alleles into per-sample, per-allele table
# hla_parsed: sample, allele (class I only)
hla_parsed <- as.data.table(hla_data_filtered)[
  , .(allele = unlist(strsplit(alleles_str, "\\|"))), by = sample
][
  grepl("^HLA-[ABC]", allele)
][
  , .(sample, allele)
][
  !is.na(allele) & allele != ""
][
  , .SD[!duplicated(.SD)], .SDcols = c("sample", "allele")
]

cat("Prepared HLA data for", nrow(hla_parsed), "sample-allele pairs from OptiType.\n")

###########################################################################
# Load MHCflurry Outputs -- ALT side ---------------------------------------
###########################################################################

setwd(directory_14_in)

# Standardise column names to what the rest of the script expects
standardize_cols <- function(dt) {
  if ("mhcflurry_affinity" %in% names(dt) && !("mhcflurry_binding_affinity" %in% names(dt))) {
    data.table::setnames(dt, "mhcflurry_affinity", "mhcflurry_binding_affinity")
  }
  dt
}

alt_file_08 <- latest_file("^08mers_flank_mhcflurry_[0-9_]+\\.csv$", directory_14_in)
alt_file_09 <- latest_file("^09mers_flank_mhcflurry_[0-9_]+\\.csv$", directory_14_in)
alt_file_10 <- latest_file("^10mers_flank_mhcflurry_[0-9_]+\\.csv$", directory_14_in)
alt_file_11 <- latest_file("^11mers_flank_mhcflurry_[0-9_]+\\.csv$", directory_14_in)
cat("ALT MHCflurry files (newest match):\n",
    alt_file_08, "\n", alt_file_09, "\n", alt_file_10, "\n", alt_file_11, "\n\n")

mhc_08 <- standardize_cols(fread(alt_file_08, na = c("", "NA")))
mhc_09 <- standardize_cols(fread(alt_file_09, na = c("", "NA")))
mhc_10 <- standardize_cols(fread(alt_file_10, na = c("", "NA")))
mhc_11 <- standardize_cols(fread(alt_file_11, na = c("", "NA")))

###########################################################################
# Load MHCflurry Outputs -- WT/self side (NEW) -----------------------------
###########################################################################

wt_file_08 <- latest_file("^08mers_flank_mhcflurry_wt_[0-9_]+\\.csv$", directory_14_in)
wt_file_09 <- latest_file("^09mers_flank_mhcflurry_wt_[0-9_]+\\.csv$", directory_14_in)
wt_file_10 <- latest_file("^10mers_flank_mhcflurry_wt_[0-9_]+\\.csv$", directory_14_in)
wt_file_11 <- latest_file("^11mers_flank_mhcflurry_wt_[0-9_]+\\.csv$", directory_14_in)
cat("WT MHCflurry files (newest match):\n",
    wt_file_08, "\n", wt_file_09, "\n", wt_file_10, "\n", wt_file_11, "\n\n")

###########################################################################
# Step 1: Cohort-Wide Top per Peptide -- ALT ------------------------------
###########################################################################

for (dt in list(mhc_08, mhc_09, mhc_10, mhc_11)) {
  dt[is.na(dt)] <- ""
}

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1))
cat("Running with", num_cores, "cores (data.table will use these natively)\n")

process_top <- function(dt, score_col = "mhcflurry_presentation_score") {
  if (nrow(dt) == 0) return(tibble())
  # Use data.table's native grouped operation instead of splitting into
  # parallel chunks -- foreach copies data into each worker's memory,
  # multiplying peak usage by num_cores. data.table does this in-place.
  setDT(dt)
  result <- dt[, .SD[which.max(get(score_col))],
               by = .(peptide, n_flank, c_flank),
               .SDcols = c("allele", score_col, "mhcflurry_binding_affinity")]
  as_tibble(result)
}

per_sample_max <- function(top_df, hla_df) {
  if (nrow(top_df) == 0 || nrow(hla_df) == 0) {
    return(tibble())
  }
  top_df %>%
    left_join(hla_df, by = "allele") %>%
    group_by(peptide, n_flank, c_flank, sample) %>%
    summarise(
      max_presentation      = max(mhcflurry_presentation_score, na.rm = TRUE),
      best_allele_for_sample = allele[which.max(mhcflurry_presentation_score)],
      binding_affinity      = mhcflurry_binding_affinity[which.max(mhcflurry_presentation_score)],
      .groups = "drop"
    ) %>%
    filter(!is.na(sample))
}

mhc_08_top <- process_top(mhc_08)
mhc_09_top <- process_top(mhc_09)
mhc_10_top <- process_top(mhc_10)
mhc_11_top <- process_top(mhc_11)

mhc_08_sample <- per_sample_max(mhc_08_top, hla_parsed)
mhc_09_sample <- per_sample_max(mhc_09_top, hla_parsed)
mhc_10_sample <- per_sample_max(mhc_10_top, hla_parsed)
mhc_11_sample <- per_sample_max(mhc_11_top, hla_parsed)

setwd(directory_14_out)
current_date <- format(Sys.Date(), "%Y%m%d")

fwrite(mhc_08_top, paste0("mhcflurry_08mer_top_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_09_top, paste0("mhcflurry_09mer_top_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_10_top, paste0("mhcflurry_10mer_top_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_11_top, paste0("mhcflurry_11mer_top_", current_date, ".tsv"), sep = "\t")

fwrite(mhc_08_sample, paste0("mhcflurry_08mer_per_sample_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_09_sample, paste0("mhcflurry_09mer_per_sample_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_10_sample, paste0("mhcflurry_10mer_per_sample_", current_date, ".tsv"), sep = "\t")
fwrite(mhc_11_sample, paste0("mhcflurry_11mer_per_sample_", current_date, ".tsv"), sep = "\t")

cat("ALT tops exported. Starting WT (much larger files) -- one length at a time to keep peak memory down...\n")

# Free the ALT objects now that they're written -- not needed anymore, and
# every bit of memory freed here is memory available for the WT files.
rm(mhc_08, mhc_09, mhc_10, mhc_11,
   mhc_08_top, mhc_09_top, mhc_10_top, mhc_11_top,
   mhc_08_sample, mhc_09_sample, mhc_10_sample, mhc_11_sample)
gc()

###########################################################################
# Step 2: WT/self side -- ONE LENGTH AT A TIME (NEW) -----------------------
###########################################################################
# Previously all 4 WT files (17-18GB each) were loaded into memory
# simultaneously, which is what caused repeated out-of-memory failures --
# no amount of requesting more cluster memory was going to fix that cleanly,
# since actual usage (with data.table/tibble copies during processing) is
# well above the raw file size, multiple times over. Processing one length
# at a time, writing its output immediately, and freeing it before loading
# the next cuts peak memory by roughly 4x. This is the real fix -- the
# SLURM resource numbers matter far less once peak memory is this much
# smaller.

wt_files <- list(
  "08" = wt_file_08,
  "09" = wt_file_09,
  "10" = wt_file_10,
  "11" = wt_file_11
)

for (len in names(wt_files)) {
  cat("\n--- WT", len, "mers ---\n")
  wt_path <- wt_files[[len]]

  mhc_wt <- standardize_cols(fread(wt_path, na = c("", "NA")))
  mhc_wt[is.na(mhc_wt)] <- ""
  cat("  Loaded", nrow(mhc_wt), "rows from", basename(wt_path), "\n")

  mhc_wt_top <- process_top(mhc_wt)
  rm(mhc_wt); gc()
  cat("  Cohort-wide top computed:", nrow(mhc_wt_top), "rows\n")

  mhc_wt_sample <- per_sample_max(mhc_wt_top, hla_parsed)
  cat("  Per-sample top computed:", nrow(mhc_wt_sample), "rows\n")

  fwrite(mhc_wt_top, paste0("mhcflurry_", len, "mer_wt_top_", current_date, ".tsv"), sep = "\t")
  fwrite(mhc_wt_sample, paste0("mhcflurry_", len, "mer_wt_per_sample_", current_date, ".tsv"), sep = "\t")
  cat("  Wrote WT", len, "mer output files.\n")

  rm(mhc_wt_top, mhc_wt_sample); gc()
}

print("Step 14c: Tops exported for BOTH ALT and WT, using HLA alleles extracted from OptiType per sample.")# Replace NA with "" (vectorized for efficiency)
nmer_08[is.na(nmer_08)] <- ""
nmer_09[is.na(nmer_09)] <- ""
nmer_10[is.na(nmer_10)] <- ""
nmer_11[is.na(nmer_11)] <- ""

# Set up parallel processing (use SLURM cores if available)
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1))
num_cores <- max(1, num_cores)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Function to process one n-mer dataset in parallel
process_nmer <- function(dt) {
  # Split into chunks for parallel processing
  chunk_size <- ceiling(nrow(dt) / num_cores)
  indices <- split(1:nrow(dt), ceiling(seq_along(1:nrow(dt)) / chunk_size))
  
  # Parallel foreach over chunks
  result_list <- foreach(idx = indices, .combine = rbind, .packages = c("data.table")) %dopar% {
    dt_chunk <- dt[idx, ]
    dt_chunk[, .SD[which.max(mhcflurry_presentation_score)], by = .(peptide, n_flank, c_flank)]
  }
  
  result_list
}

# Process each n-mer in parallel
nmer_08_edit <- process_nmer(nmer_08)
nmer_09_edit <- process_nmer(nmer_09)
nmer_10_edit <- process_nmer(nmer_10)
nmer_11_edit <- process_nmer(nmer_11)

# Stop the cluster
stopCluster(cl)

# Export Files ------------------------------------------------------------
setwd(directory_14_out)
current_date <- format(Sys.Date(), "%Y%m%d")  # Use today's date for output naming
fwrite(nmer_08_edit, paste0("mhcflurry_08mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_09_edit, paste0("mhcflurry_09mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_10_edit, paste0("mhcflurry_10mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
fwrite(nmer_11_edit, paste0("mhcflurry_11mer_selected_alleles_", current_date, ".tsv"), sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)

print("Step 14c completed: Selected top HLA-alleles for each n-mer and saved output files.")
