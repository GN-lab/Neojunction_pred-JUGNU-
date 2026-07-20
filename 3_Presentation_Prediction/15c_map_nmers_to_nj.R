#!/usr/bin/env Rscript
# Step 15c: Map n-mers back to neojunctions (ALT + WT)
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: For each peptide in the concordance tiers (from 15a), find which
#          neojunction it came from by searching it as a substring in aa.seq.alt
#          from 2023_0812_complete_list_all_mers.tsv.
#
# Key insight from old pipeline:
#   junc.id is obtained by str_detect(aa.seq.alt, fixed(peptide)) against
#   complete_list_all_mers.tsv -- NOT from Res_AA_Prediction_Confirmed.
#   complete_list_all_mers.tsv has ALL n-mers with junction metadata already.
#
# Input (15a): alt_concordance_all_YYYYMMDD.tsv
#              wt_concordance_tier1_YYYYMMDD.tsv
# Input (results): 2023_0812_complete_list_all_mers.tsv
#
# Output:
#   alt_neoA_to_neoJ_map_YYYYMMDD.tsv     -- peptide x junction (all ALT tiers)
#   alt_immunogenic_njs_YYYYMMDD.tsv       -- per-junction summary
#   wt_neoA_to_neoJ_map_YYYYMMDD.tsv      -- peptide x junction (WT Tier 1)
#   wt_immunogenic_njs_YYYYMMDD.tsv        -- per-junction summary

###########################################################################
#  Step 0: Packages and config
###########################################################################

rm(list = ls(all.names = TRUE))
library(data.table)
library(stringr)

output_dir   <- Sys.getenv("OUTPUT_DIR")
directory_15 <- Sys.getenv("STEP15_OUTPUT_DIR")

if (nchar(output_dir)   == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_15) == 0) directory_15 <- output_dir

###########################################################################
#  Step 1: Detect dates
###########################################################################

alt_scan <- list.files(directory_15,
                        pattern = "^alt_concordance_all_[0-9]{8}\\.tsv$")
if (length(alt_scan) == 0)
  stop("[ERROR] No alt_concordance_all_YYYYMMDD.tsv in: ", directory_15)
date_15a <- sub("^alt_concordance_all_([0-9]{8})\\.tsv$", "\\1",
                alt_scan[which.max(
                  file.info(file.path(directory_15, alt_scan))$mtime)])

current_date <- format(Sys.Date(), "%Y%m%d")
cat("[INFO] 15a date:    ", date_15a,    "\n")
cat("[INFO] Output date: ", current_date, "\n")

###########################################################################
#  Step 2: Load complete_list_all_mers.tsv -- the junction reference
#  This file has junc.id, symbol, type, fs, aa.seq.alt for every junction
###########################################################################

cat("\n[STEP 2] Loading complete_list_all_mers...\n")
mers_file <- file.path(output_dir, "2023_0812_complete_list_all_mers.tsv")
if (!file.exists(mers_file))
  stop("[ERROR] File not found: ", mers_file)

mers <- fread(mers_file, na.strings = c("", "NA"), quote = "",
              select = c("junc.id", "symbol", "type",
                         "aa.change", "ln.diff", "aa.seq.alt"))
cat(sprintf("[INFO] complete_list_all_mers: %d rows\n", nrow(mers)))

# Derive fs from aa.change and ln.diff (same logic as old script)
mers[, fs := fifelse(
  grepl("shift|fs", aa.change, ignore.case = TRUE) | (ln.diff %% 3 != 0),
  "fs", "in-frame"
)]

# Keep one row per junc.id (some junctions have multiple transcripts)
# Keep the longest aa.seq.alt for best substring coverage
mers[, seq_len := nchar(as.character(aa.seq.alt))]
setorder(mers, junc.id, -seq_len)
mers_dedup <- mers[, .SD[1], by = junc.id]
cat(sprintf("[INFO] Unique junctions: %d\n", nrow(mers_dedup)))

###########################################################################
#  Helper: map peptides to junctions via substring search
###########################################################################

map_peptides_to_junctions <- function(candidates, label) {

  if (nrow(candidates) == 0) {
    cat(sprintf("[WARN] No candidates for: %s\n", label))
    return(list(map = data.table(), nj_summary = data.table()))
  }

  # Get unique peptides to avoid redundant searches
  unique_peps <- unique(candidates$peptide)
  cat(sprintf("[INFO] %s: %d unique peptides to map against %d junctions\n",
              label, length(unique_peps), nrow(mers_dedup)))

  # For each unique peptide, find which junctions contain it in aa.seq.alt
  cat(sprintf("[INFO] Running substring search for %s...\n", label))

  pep_to_junc <- rbindlist(lapply(unique_peps, function(pep) {
    hits <- mers_dedup[str_detect(aa.seq.alt, fixed(pep))]
    if (nrow(hits) == 0) return(NULL)
    hits[, peptide := pep]
    hits[, .(peptide, junc.id, symbol, type, fs)]
  }), fill = TRUE)

  if (nrow(pep_to_junc) == 0) {
    cat(sprintf("[WARN] No junction hits found for: %s\n", label))
    return(list(map = data.table(), nj_summary = data.table()))
  }
  cat(sprintf("[INFO] %s: %d peptide-junction hits\n", label, nrow(pep_to_junc)))

  # Join back all allele/tier/score info from candidates
  score_cols <- intersect(
    c("peptide", "allele", "concordance_tier",
      "netmhcpan_EL_score", "netmhcpan_EL_rank",
      "mhcflurry_affinity", "mhcflurry_presentation_score",
      "combined_score"),
    colnames(candidates)
  )
  cand_meta <- unique(candidates[, score_cols, with = FALSE])

  # Join peptide-junction hits with candidate scores
  mapped <- merge(pep_to_junc, cand_meta,
                  by = "peptide", all.x = TRUE, allow.cartesian = TRUE)

  # Build neo_id: peptide|junc.id|allele
  mapped[, neo_id := paste(peptide, junc.id, allele, sep = "|")]

  cat(sprintf("[INFO] %s: %d total rows after joining scores\n",
              label, nrow(mapped)))

  # Per-junction summary
  nj_summary <- mapped[, .(
    n_tier1_peptides    = sum(concordance_tier == "Tier1_HighConfidence",  na.rm=TRUE),
    n_tier2_peptides    = sum(concordance_tier == "Tier2_MediumConfidence", na.rm=TRUE),
    n_tier3_peptides    = sum(grepl("Tier3", concordance_tier), na.rm=TRUE),
    n_total_peptides    = .N,
    n_alleles           = uniqueN(allele),
    best_combined_score = max(combined_score,       na.rm=TRUE),
    best_EL_rank        = min(netmhcpan_EL_rank,    na.rm=TRUE),
    best_mhc_affinity   = min(mhcflurry_affinity,   na.rm=TRUE)
  ), by = .(junc.id, symbol, type, fs)]

  nj_summary[, top_tier := fcase(
    n_tier1_peptides > 0, "Tier1_HighConfidence",
    n_tier2_peptides > 0, "Tier2_MediumConfidence",
    n_tier3_peptides > 0, "Tier3_Discordant",
    default = "None"
  )]

  setorder(nj_summary, -n_tier1_peptides, -best_combined_score)

  list(map = mapped, nj_summary = nj_summary)
}

###########################################################################
#  Step 3: Map ALT candidates (all tiers)
###########################################################################

cat("\n[STEP 3] Loading and mapping ALT concordance results...\n")
alt_all <- fread(file.path(directory_15,
                            paste0("alt_concordance_all_", date_15a, ".tsv")),
                 na.strings = c("", "NA"), quote = "")

# Strip quotes from string columns if present
for (col in c("allele","peptide","binder","nmp_call","mhc_call","concordance_tier")) {
  if (col %in% colnames(alt_all))
    set(alt_all, j = col, value = gsub('^"|"$', "", alt_all[[col]]))
}

alt_all <- alt_all[concordance_tier != "Excluded"]
cat(sprintf("[INFO] ALT candidates: %d rows\n", nrow(alt_all)))
print(alt_all[, .N, by = concordance_tier][order(concordance_tier)])

alt_result <- map_peptides_to_junctions(alt_all, "ALT")

###########################################################################
#  Step 4: Map WT Tier 1
###########################################################################

cat("\n[STEP 4] Loading and mapping WT Tier 1...\n")
wt_file <- file.path(directory_15,
                      paste0("wt_concordance_tier1_", date_15a, ".tsv"))
if (file.exists(wt_file)) {
  wt_t1 <- fread(wt_file, na.strings = c("", "NA"), quote = "")
  for (col in c("allele","peptide","binder","nmp_call","mhc_call","concordance_tier")) {
    if (col %in% colnames(wt_t1))
      set(wt_t1, j = col, value = gsub('^"|"$', "", wt_t1[[col]]))
  }
  cat(sprintf("[INFO] WT Tier 1 candidates: %d rows\n", nrow(wt_t1)))
  wt_result <- map_peptides_to_junctions(wt_t1, "WT Tier1")
} else {
  cat("[WARN] WT Tier 1 file not found -- skipping\n")
  wt_result <- list(map = data.table(), nj_summary = data.table())
}

###########################################################################
#  Step 5: Write outputs
###########################################################################

cat("\n[STEP 5] Writing outputs...\n")
setwd(directory_15)

if (nrow(alt_result$map) > 0) {
  fwrite(alt_result$map,
         paste0("alt_neoA_to_neoJ_map_", current_date, ".tsv"),
         sep = "\t", na = "NA", quote = FALSE)
  fwrite(alt_result$nj_summary,
         paste0("alt_immunogenic_njs_", current_date, ".tsv"),
         sep = "\t", na = "NA", quote = FALSE)
  cat(sprintf("[INFO] ALT map: %d rows -> alt_neoA_to_neoJ_map_%s.tsv\n",
              nrow(alt_result$map), current_date))
  cat(sprintf("[INFO] ALT NJs: %d -> alt_immunogenic_njs_%s.tsv\n",
              nrow(alt_result$nj_summary), current_date))
}

if (nrow(wt_result$map) > 0) {
  fwrite(wt_result$map,
         paste0("wt_neoA_to_neoJ_map_", current_date, ".tsv"),
         sep = "\t", na = "NA", quote = FALSE)
  fwrite(wt_result$nj_summary,
         paste0("wt_immunogenic_njs_", current_date, ".tsv"),
         sep = "\t", na = "NA", quote = FALSE)
  cat(sprintf("[INFO] WT map:  %d rows -> wt_neoA_to_neoJ_map_%s.tsv\n",
              nrow(wt_result$map), current_date))
  cat(sprintf("[INFO] WT NJs:  %d -> wt_immunogenic_njs_%s.tsv\n",
              nrow(wt_result$nj_summary), current_date))
}

###########################################################################
#  Summary
###########################################################################

cat("\n=== 15c Summary ===\n")
cat(sprintf("  ALT mapped rows:            %d\n", nrow(alt_result$map)))
cat(sprintf("  ALT immunogenic junctions:  %d\n", nrow(alt_result$nj_summary)))
if (nrow(alt_result$nj_summary) > 0) {
  cat("  ALT top tier breakdown:\n")
  print(alt_result$nj_summary[, .N, by = top_tier][order(top_tier)])
}
cat(sprintf("  WT Tier1 mapped rows:       %d\n", nrow(wt_result$map)))
cat(sprintf("  WT immunogenic junctions:   %d\n", nrow(wt_result$nj_summary)))
cat("\n[DONE] Step 15c complete.\n")