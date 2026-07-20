#!/usr/bin/env Rscript
# Step 15a: Cross-analysis -- NetMHCpan 4.2 x MHCflurry 2.0 Concordance
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# TWO-STEP LOGIC:
#
# STEP 1 -- Clean ALT immunopeptidome:
#   Load ALL ALT predictions (NetMHCpan + MHCflurry, full universe).
#   Remove ANY peptide sequence found in WT predictions from EITHER tool
#   at ANY binding level (NB, WB, SB -- doesn't matter).
#   What remains is the clean tumour-specific peptide universe.
#
# STEP 2 -- Concordance tiering (applied separately to ALT and WT):
#   Tier 1 (High confidence):   NMP EL rank < 0.5% AND MHCflurry < 500nM
#   Tier 2 (Medium confidence): NMP EL rank < 2.0% AND MHCflurry < 500nM
#   Tier 3 (Discordant):        One tool binder, other NB -- flagged
#   Excluded:                   Both NB -- dropped
#
# Output files:
#   ALT (tumour-specific, WT-filtered):
#     alt_concordance_tier1_YYYYMMDD.tsv   -- both tools strong
#     alt_concordance_tier2_YYYYMMDD.tsv   -- moderate agreement
#     alt_concordance_tier3_YYYYMMDD.tsv   -- discordant
#     alt_concordance_all_YYYYMMDD.tsv     -- all tiers combined
#   WT (native immunopeptidome reference):
#     wt_concordance_tier1_YYYYMMDD.tsv
#     wt_concordance_tier2_YYYYMMDD.tsv
#     wt_concordance_tier3_YYYYMMDD.tsv
#     wt_concordance_all_YYYYMMDD.tsv
#   Shared:
#     cross_alg_all_nmers_YYYYMMDD.tsv    -- clean ALT full join (for 15b/c/d)
#     wt_exclusion_peptides_YYYYMMDD.txt  -- list of excluded WT peptides

###########################################################################
#  Step 0: Packages and config
###########################################################################

rm(list = ls(all.names = TRUE))
library(data.table)

output_dir   <- Sys.getenv("OUTPUT_DIR")
directory_13 <- Sys.getenv("STEP13_OUTPUT_DIR")
directory_14 <- Sys.getenv("STEP14_OUTPUT_DIR")
directory_15 <- Sys.getenv("STEP15_OUTPUT_DIR")

if (nchar(output_dir)   == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_13) == 0) directory_13 <- output_dir
if (nchar(directory_14) == 0) directory_14 <- output_dir
if (nchar(directory_15) == 0) directory_15 <- output_dir

dir.create(directory_15, showWarnings = FALSE, recursive = TRUE)

# Binding thresholds
NMP_SB_RANK <- 0.5    # NetMHCpan EL rank strictly < this = Strong Binder
NMP_WB_RANK <- 2.0    # NetMHCpan EL rank strictly < this = Weak Binder
MHC_AFF_NM  <- 500    # MHCflurry affinity strictly < this nM = binder

###########################################################################
#  Step 1: Detect file dates from existing files
###########################################################################

# 13a raw date: netmhcpan_09mer_YYYY_MMDD.tsv
nmp_raw_scan <- list.files(directory_13,
                            pattern = "^netmhcpan_09mer_[0-9]{4}_[0-9]{4}\\.tsv$")
if (length(nmp_raw_scan) == 0)
  stop("[ERROR] No netmhcpan_09mer_YYYY_MMDD.tsv in: ", directory_13)
nmp_run_date <- sub("^netmhcpan_09mer_([0-9]{4}_[0-9]{4})\\.tsv$", "\\1",
                    nmp_raw_scan[which.max(
                      file.info(file.path(directory_13, nmp_raw_scan))$mtime)])

# 14b ALT date: 09mers_flank_mhcflurry_YYYY_MMDD.csv
mhc_scan <- list.files(directory_14,
                        pattern = "^09mers_flank_mhcflurry_[0-9]{4}_[0-9]{4}\\.csv$")
if (length(mhc_scan) == 0)
  stop("[ERROR] No 09mers_flank_mhcflurry_YYYY_MMDD.csv in: ", directory_14)
mhc_date <- sub("^09mers_flank_mhcflurry_([0-9]{4}_[0-9]{4})\\.csv$", "\\1",
                mhc_scan[which.max(
                  file.info(file.path(directory_14, mhc_scan))$mtime)])

# 14b WT date: 09mers_flank_mhcflurry_wt_YYYY_MMDD.csv
mhc_wt_scan <- list.files(directory_14,
                           pattern = "^09mers_flank_mhcflurry_wt_[0-9]{4}_[0-9]{4}\\.csv$")
if (length(mhc_wt_scan) == 0)
  stop("[ERROR] No 09mers_flank_mhcflurry_wt_YYYY_MMDD.csv in: ", directory_14)
mhc_wt_date <- sub("^09mers_flank_mhcflurry_wt_([0-9]{4}_[0-9]{4})\\.csv$", "\\1",
                   mhc_wt_scan[which.max(
                     file.info(file.path(directory_14, mhc_wt_scan))$mtime)])

current_date <- format(Sys.Date(), "%Y%m%d")

cat("=== File dates detected ===\n")
cat("[INFO] NetMHCpan raw (13a):      ", nmp_run_date, "\n")
cat("[INFO] MHCFlurry ALT (14b):      ", mhc_date,     "\n")
cat("[INFO] MHCFlurry WT  (14b):      ", mhc_wt_date,  "\n")
cat("[INFO] Output date:              ", current_date,  "\n")
cat("[INFO] Tier 1: NMP EL rank <",  NMP_SB_RANK, "% AND MHCflurry <", MHC_AFF_NM, "nM\n")
cat("[INFO] Tier 2: NMP EL rank <",  NMP_WB_RANK, "% AND MHCflurry <", MHC_AFF_NM, "nM\n")

###########################################################################
#  Helper: load NetMHCpan raw TSV (all binding levels)
###########################################################################

load_nmp_raw <- function(len, type = "alt") {
  if (type == "wt") {
    f <- file.path(directory_13,
                   paste0("netmhcpan_", len, "mer_wt_", nmp_run_date, ".tsv"))
  } else {
    f <- file.path(directory_13,
                   paste0("netmhcpan_", len, "mer_", nmp_run_date, ".tsv"))
  }
  if (!file.exists(f)) { warning("[WARN] Missing: ", f); return(data.table()) }
  cat(sprintf("[INFO]   Loading NMP %s %smer...\n", type, len))
  dt <- fread(f, na.strings = c("", "NA"),
              select = c("allele", "peptide", "netmhcpan_EL_score",
                         "netmhcpan_EL_rank", "netmhcpan_BA_score",
                         "netmhcpan_BA_rank", "binder"))
  dt[, netmhcpan_EL_rank  := as.numeric(netmhcpan_EL_rank)]
  dt[, netmhcpan_EL_score := as.numeric(netmhcpan_EL_score)]
  dt
}

###########################################################################
#  Helper: load MHCflurry raw flank CSV (full prediction universe)
###########################################################################

load_mhc_raw <- function(len, type = "alt") {
  if (type == "wt") {
    f <- file.path(directory_14,
                   paste0(len, "mers_flank_mhcflurry_wt_", mhc_wt_date, ".csv"))
  } else {
    f <- file.path(directory_14,
                   paste0(len, "mers_flank_mhcflurry_", mhc_date, ".csv"))
  }
  if (!file.exists(f)) { warning("[WARN] Missing: ", f); return(data.table()) }
  cat(sprintf("[INFO]   Loading MHC %s %smer...\n", type, len))
  dt <- fread(f, na.strings = c("", "NA"),
              select = c("allele", "peptide", "n_flank", "c_flank",
                         "mhcflurry_affinity", "mhcflurry_presentation_score",
                         "mhcflurry_presentation_percentile"))
  dt[, mhcflurry_affinity           := as.numeric(mhcflurry_affinity)]
  dt[, mhcflurry_presentation_score := as.numeric(mhcflurry_presentation_score)]
  dt
}

###########################################################################
#  Helper: assign concordance tiers to a joined data.table
###########################################################################

assign_tiers <- function(dt) {
  dt[, nmp_EL_rank := as.numeric(netmhcpan_EL_rank)]
  dt[, mhc_aff     := as.numeric(mhcflurry_affinity)]

  # NMP call -- strict < thresholds
  dt[, nmp_call := fcase(
    nmp_EL_rank < NMP_SB_RANK, "SB",
    nmp_EL_rank < NMP_WB_RANK, "WB",
    default = "NB"
  )]

  # MHCflurry call -- strict < 500nM
  dt[, mhc_call := fcase(
    mhc_aff < MHC_AFF_NM, "binder",
    default = "NB"
  )]

  # Concordance tier
  dt[, concordance_tier := fcase(
    nmp_call == "SB" & mhc_call == "binder",              "Tier1_HighConfidence",
    nmp_call == "WB" & mhc_call == "binder",              "Tier2_MediumConfidence",
    nmp_call == "SB" & mhc_call == "NB",                  "Tier3_Discordant_NMPstrong",
    nmp_call %in% c("WB","NB") & mhc_call == "binder",   "Tier3_Discordant_MHCstrong",
    default = "Excluded"
  )]

  # Combined score: mean of NMP EL score + MHCflurry presentation score (both 0-1)
  dt[, combined_score := rowMeans(
    cbind(as.numeric(netmhcpan_EL_score),
          as.numeric(mhcflurry_presentation_score)),
    na.rm = TRUE
  )]
  dt
}

###########################################################################
#  Helper: join NMP + MHCflurry on allele + peptide, keep best row each
###########################################################################

join_tools <- function(nmp_dt, mhc_dt) {
  # Best row per allele+peptide in each tool
  nmp_best <- nmp_dt[order(netmhcpan_EL_rank)][, .SD[1], by = .(allele, peptide)]
  mhc_best <- mhc_dt[order(mhcflurry_affinity)][, .SD[1], by = .(allele, peptide)]
  merge(nmp_best, mhc_best, by = c("allele", "peptide"), all = TRUE)
}

###########################################################################
#  Helper: split tiered table into files and write
###########################################################################

write_tiers <- function(dt, prefix) {
  setwd(directory_15)

  tier1 <- dt[concordance_tier == "Tier1_HighConfidence"]
  tier2 <- dt[concordance_tier == "Tier2_MediumConfidence"]
  tier3 <- dt[grepl("Tier3", concordance_tier)]
  all_t <- dt[concordance_tier != "Excluded"]

  setorder(tier1, nmp_EL_rank)
  setorder(tier2, nmp_EL_rank)
  setorder(tier3, concordance_tier, nmp_EL_rank)
  setorder(all_t, concordance_tier, nmp_EL_rank)

  fwrite(tier1, paste0(prefix, "_tier1_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(tier2, paste0(prefix, "_tier2_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(tier3, paste0(prefix, "_tier3_", current_date, ".tsv"), sep="\t", na="NA")
  fwrite(all_t, paste0(prefix, "_all_",   current_date, ".tsv"), sep="\t", na="NA")

  cat(sprintf("\n  [%s] Tier 1 (High confidence):   %d rows\n", prefix, nrow(tier1)))
  cat(sprintf("  [%s] Tier 2 (Medium confidence): %d rows\n", prefix, nrow(tier2)))
  cat(sprintf("  [%s] Tier 3 (Discordant):        %d rows\n", prefix, nrow(tier3)))
  cat(sprintf("  [%s] Excluded (both NB):         %d rows\n", prefix,
              nrow(dt[concordance_tier == "Excluded"])))

  list(tier1=tier1, tier2=tier2, tier3=tier3, all=all_t)
}

###########################################################################
#  ======================== STEP 1: BUILD WT EXCLUSION SET ===============
#  Load ALL WT predictions from BOTH tools at ALL binding levels.
#  Every peptide sequence seen in WT is added to the exclusion set.
#  This is done BEFORE any ALT processing.
###########################################################################

cat("\n========== STEP 1: Building WT exclusion set ==========\n")

cat("[INFO] Loading NetMHCpan WT (all levels)...\n")
nmp_wt_all <- rbindlist(
  lapply(c("08","09","10","11"), load_nmp_raw, type = "wt"), fill = TRUE)
cat(sprintf("[INFO] NetMHCpan WT total rows: %d\n", nrow(nmp_wt_all)))

cat("[INFO] Loading MHCflurry WT (all levels)...\n")
mhc_wt_all <- rbindlist(
  lapply(c("08","09","10","11"), load_mhc_raw, type = "wt"), fill = TRUE)
cat(sprintf("[INFO] MHCflurry WT total rows: %d\n", nrow(mhc_wt_all)))

# Union of ALL WT peptide sequences from both tools
wt_exclusion <- unique(c(
  unique(nmp_wt_all$peptide),
  unique(mhc_wt_all$peptide)
))
cat(sprintf("\n[INFO] Total unique WT peptides to exclude: %d\n", length(wt_exclusion)))

# Write exclusion list for traceability
setwd(directory_15)
writeLines(wt_exclusion,
           paste0("wt_exclusion_peptides_", current_date, ".txt"))
cat(sprintf("[INFO] WT exclusion list written: wt_exclusion_peptides_%s.txt\n",
            current_date))

###########################################################################
#  ======================== STEP 2a: CLEAN ALT IMMUNOPEPTIDOME ===========
#  Load ALL ALT predictions, remove WT peptides, then tier what remains.
###########################################################################

cat("\n========== STEP 2a: Clean ALT immunopeptidome ==========\n")

cat("[INFO] Loading NetMHCpan ALT (all levels)...\n")
nmp_alt_all <- rbindlist(
  lapply(c("08","09","10","11"), load_nmp_raw, type = "alt"), fill = TRUE)
cat(sprintf("[INFO] NetMHCpan ALT total rows: %d\n", nrow(nmp_alt_all)))

cat("[INFO] Loading MHCflurry ALT (all levels)...\n")
mhc_alt_all <- rbindlist(
  lapply(c("08","09","10","11"), load_mhc_raw, type = "alt"), fill = TRUE)
cat(sprintf("[INFO] MHCflurry ALT total rows: %d\n", nrow(mhc_alt_all)))

# Apply WT exclusion filter FIRST -- before any tiering
nmp_alt_clean <- nmp_alt_all[!peptide %in% wt_exclusion]
mhc_alt_clean <- mhc_alt_all[!peptide %in% wt_exclusion]

cat(sprintf("\n[INFO] NMP ALT after WT exclusion: %d rows (removed %d)\n",
            nrow(nmp_alt_clean), nrow(nmp_alt_all) - nrow(nmp_alt_clean)))
cat(sprintf("[INFO] MHC ALT after WT exclusion: %d rows (removed %d)\n",
            nrow(mhc_alt_clean), nrow(mhc_alt_all) - nrow(mhc_alt_clean)))

# Now join and tier the clean ALT set
cat("\n[INFO] Joining clean ALT NMP + MHCflurry...\n")
alt_joined <- join_tools(nmp_alt_clean, mhc_alt_clean)
cat(sprintf("[INFO] Clean ALT combined rows: %d\n", nrow(alt_joined)))

alt_tiered <- assign_tiers(alt_joined)

cat("\n--- ALT concordance counts ---\n")
print(alt_tiered[, .N, by = concordance_tier][order(-N)])

cat("\n[INFO] Writing ALT tier files...\n")
alt_results <- write_tiers(alt_tiered, "alt_concordance")

# Write full joined table for downstream 15b/c/d compatibility
setwd(directory_15)
fwrite(alt_tiered[concordance_tier != "Excluded"],
       paste0("cross_alg_all_nmers_", current_date, ".tsv"),
       sep="\t", na="NA")

###########################################################################
#  ======================== STEP 2b: WT NATIVE IMMUNOPEPTIDOME ===========
#  Tier the WT predictions independently as a reference set.
###########################################################################

cat("\n========== STEP 2b: WT native immunopeptidome ==========\n")

cat("\n[INFO] Joining WT NMP + MHCflurry...\n")
wt_joined  <- join_tools(nmp_wt_all, mhc_wt_all)
cat(sprintf("[INFO] WT combined rows: %d\n", nrow(wt_joined)))

wt_tiered  <- assign_tiers(wt_joined)

cat("\n--- WT concordance counts ---\n")
print(wt_tiered[, .N, by = concordance_tier][order(-N)])

cat("\n[INFO] Writing WT tier files...\n")
wt_results <- write_tiers(wt_tiered, "wt_concordance")

###########################################################################
#  Summary
###########################################################################

cat("\n=== 15a Summary ===\n")
cat("\nALT (tumour-specific, WT-filtered):\n")
cat(sprintf("  Tier 1 -- High confidence (NMP SB + MHC <500nM):  %d\n",
            nrow(alt_results$tier1)))
cat(sprintf("  Tier 2 -- Medium confidence (NMP WB + MHC <500nM): %d\n",
            nrow(alt_results$tier2)))
cat(sprintf("  Tier 3 -- Discordant (tools disagree):             %d\n",
            nrow(alt_results$tier3)))

cat("\nWT native immunopeptidome (reference):\n")
cat(sprintf("  Tier 1 -- High confidence:   %d\n", nrow(wt_results$tier1)))
cat(sprintf("  Tier 2 -- Medium confidence: %d\n", nrow(wt_results$tier2)))
cat(sprintf("  Tier 3 -- Discordant:        %d\n", nrow(wt_results$tier3)))

cat(sprintf("\n  WT peptides excluded from ALT: %d\n", length(wt_exclusion)))
cat(sprintf("  All files written to: %s\n", directory_15))
cat("\n[DONE] Step 15a complete.\n")