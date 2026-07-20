#!/usr/bin/env Rscript
# Step 13b: Select TOP HLA-allele presentation score for each n-mer (NetMHCpan)
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: For each n-mer predicted by NetMHCpan 4.2, keep only strong binder
#          (SB) rows, then select the best HLA allele per peptide using the
#          lowest EL rank (= best presentation likelihood).
#
# Input:   netmhcpan_XXmer_YYYY_MMDD.tsv         (ALT, from Step 13a)
#          netmhcpan_XXmer_wt_YYYY_MMDD.tsv      (WT,  from Step 13a)
# Output:  netmhcpan_XXmer_selected_alleles_YYYYMMDD.tsv
#          netmhcpan_XXmer_wt_selected_alleles_YYYYMMDD.tsv
#
# Key columns in input:
#   allele | peptide | core | netmhcpan_EL_score | netmhcpan_EL_rank |
#   netmhcpan_BA_score | netmhcpan_BA_rank | binder
#
# binder values: "<= SB" (strong), "<= WB" (weak), "NB" (non-binder)
# Selection logic:
#   1. Keep rows where binder == "<= SB"
#   2. For each peptide, pick the row with the lowest netmhcpan_EL_rank
#      (ties broken by lowest EL rank value -> highest EL score)
#   No additional rank cutoff applied here -- intersection with MHCflurry
#   results in a downstream step will further tighten the candidate list.

###########################################################################
#  Step 0: Load packages and config ---------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(data.table)

# Directories from environment (set in config.sh / job submission)
#
# IMPORTANT: 13a writes its TSV output to OUTPUT_DIR (main results directory).
#            13b reads from there and writes its selected-allele files to
#            STEP13_OUTPUT_DIR (which may be the same or a sub-directory).
#            If STEP13_OUTPUT_DIR is not set, fall back to OUTPUT_DIR for both.
input_dir    <- Sys.getenv("OUTPUT_DIR")       # where 13a wrote its files
directory_13 <- Sys.getenv("STEP13_OUTPUT_DIR") # where 13b writes its output

if (nchar(input_dir) == 0)    stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_13) == 0) {
  warning("[WARN] STEP13_OUTPUT_DIR not set -- writing output to OUTPUT_DIR")
  directory_13 <- input_dir
}

dir.create(directory_13, showWarnings = FALSE, recursive = TRUE)

# Detect run_date from the most recently written 09mer ALT file in input_dir.
# 13a uses RUN_DATE=$(date +%Y_%m%d) e.g. "2026_0716" -- this avoids a date
# mismatch when 13b runs on a different calendar day from 13a.
nmp_files <- list.files(input_dir, pattern = "^netmhcpan_09mer_[0-9]{4}_[0-9]{4}\\.tsv$")
if (length(nmp_files) == 0) {
  stop("[ERROR] No netmhcpan_09mer_YYYY_MMDD.tsv files found in: ", input_dir,
       "\n  Has Step 13a finished running?")
}
# Pick the most recently modified file and extract its date stamp
nmp_mtimes  <- file.info(file.path(input_dir, nmp_files))$mtime
run_date    <- sub("^netmhcpan_09mer_(.+)\\.tsv$", "\\1",
                   nmp_files[which.max(nmp_mtimes)])
# out_date for output filenames: strip underscores -> YYYYMMDD
out_date    <- gsub("_", "", run_date)

cat("[INFO] Reading 13a files from:     ", input_dir,    "\n")
cat("[INFO] Writing 13b files to:       ", directory_13, "\n")
cat("[INFO] Detected 13a run_date:      ", run_date,     "\n")
cat("[INFO] Output files will be dated: ", out_date,     "\n")

###########################################################################
#  Step 1: Helper -- load, filter SB, select best allele per peptide -----
###########################################################################

select_top_allele <- function(length_label, type = "alt") {

  # Build input filename
  if (type == "wt") {
    infile  <- file.path(input_dir,    paste0("netmhcpan_", length_label, "mer_wt_", run_date, ".tsv"))
    outfile <- file.path(directory_13, paste0("netmhcpan_", length_label, "mer_wt_selected_alleles_", out_date, ".tsv"))
  } else {
    infile  <- file.path(input_dir,    paste0("netmhcpan_", length_label, "mer_", run_date, ".tsv"))
    outfile <- file.path(directory_13, paste0("netmhcpan_", length_label, "mer_selected_alleles_", out_date, ".tsv"))
  }

  if (!file.exists(infile)) {
    warning("[WARN] File not found, skipping: ", infile)
    return(invisible(NULL))
  }

  dt <- fread(infile, na.strings = c("", "NA"))

  cat(sprintf("[INFO] %s %smer: %d rows loaded\n", type, length_label, nrow(dt)))

  # Convert EL rank to numeric (should already be, but guard against parse issues)
  dt[, netmhcpan_EL_rank := as.numeric(netmhcpan_EL_rank)]

  # Step 1: Keep strong binders only
  sb <- dt[binder %in% c("SB", "<= SB")]
  cat(sprintf("[INFO] %s %smer: %d strong binder rows (binder == 'SB' or '<= SB')\n",
              type, length_label, nrow(sb)))

  if (nrow(sb) == 0) {
    warning("[WARN] No strong binders found for ", type, " ", length_label, "mer -- output will be empty")
    fwrite(sb, outfile, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
    return(invisible(sb))
  }

  # Step 2: For each unique peptide, keep the row with the lowest EL rank
  # (i.e. the allele that presents it best)
  # If two alleles tie on EL rank, break ties by highest EL score
  setorder(sb, peptide, netmhcpan_EL_rank, -netmhcpan_EL_score)
  best <- sb[, .SD[1], by = peptide]

  cat(sprintf("[INFO] %s %smer: %d unique peptides after selecting best allele\n",
              type, length_label, nrow(best)))

  fwrite(best, outfile, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
  cat(sprintf("[INFO] Written: %s\n", outfile))

  return(invisible(best))
}

###########################################################################
#  Step 2: Run for ALT (cancer-specific) peptides ------------------------
###########################################################################

cat("\n=== ALT (cancer-specific) ===\n")
nmer_08_alt <- select_top_allele("08", "alt")
nmer_09_alt <- select_top_allele("09", "alt")
nmer_10_alt <- select_top_allele("10", "alt")
nmer_11_alt <- select_top_allele("11", "alt")

###########################################################################
#  Step 3: Run for WT (native self-peptide) predictions ------------------
###########################################################################

cat("\n=== WT (native self-peptides) ===\n")
nmer_08_wt <- select_top_allele("08", "wt")
nmer_09_wt <- select_top_allele("09", "wt")
nmer_10_wt <- select_top_allele("10", "wt")
nmer_11_wt <- select_top_allele("11", "wt")

###########################################################################
#  Step 4: Summary --------------------------------------------------------
###########################################################################

cat("\n=== 13b Summary ===\n")
for (len in c("08", "09", "10", "11")) {
  for (tp in c("alt", "wt")) {
    if (tp == "wt") {
      f <- file.path(directory_13, paste0("netmhcpan_", len, "mer_wt_selected_alleles_", out_date, ".tsv"))
    } else {
      f <- file.path(directory_13, paste0("netmhcpan_", len, "mer_selected_alleles_", out_date, ".tsv"))
    }
    if (file.exists(f)) {
      n <- nrow(fread(f, nrows = Inf, select = 1L))
      cat(sprintf("  %s %smer: %d peptides -> %s\n", tp, len, n, f))
    } else {
      cat(sprintf("  %s %smer: MISSING\n", tp, len))
    }
  }
}

cat("\n[DONE] Step 13b complete.\n")