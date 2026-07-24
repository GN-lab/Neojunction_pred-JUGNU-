#!/usr/bin/env Rscript
# Title: "Step 11b: WT Native Proteome Peptide Generation"
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Generate 8-11mer peptides from the full expressed human proteome
#          (UniProt canonical + isoforms, combined FASTA) to replace the
#          junction-only WT peptides that Step 12 currently writes.
#
#          Output files are written with EXACTLY the same filenames that
#          Steps 13a and 14a expect -- so nothing downstream needs to change:
#
#            2023_0812_hlathenalist_msic_08mers_wt.tsv   (-> 13a NetMHCpan)
#            2023_0812_hlathenalist_msic_09mers_wt.tsv
#            2023_0812_hlathenalist_msic_10mers_wt.tsv
#            2023_0812_hlathenalist_msic_11mers_wt.tsv
#
#          This script runs AFTER Step 12 and BEFORE Step 13a.
#          Step 12 still writes its own junction-level _wt files first;
#          this script overwrites them with the full-proteome versions.
#
# Config variables used (all set in config.sh):
#   PATH_TO_FASTA_UNIPROT   -- directory containing combined FASTA
#   OUTPUT_DIR              -- where Step 12 wrote its files (same dir 13a reads)
#   SLURM_CPUS_PER_TASK     -- parallelism (set by SLURM automatically)

###########################################################################
# 0. Packages
###########################################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(stringr)
  library(doParallel)
  library(foreach)
})

start_time <- proc.time()

###########################################################################
# 1. Paths from config.sh env vars
###########################################################################

fasta_dir  <- Sys.getenv("PATH_TO_FASTA_UNIPROT")   # INPUT_DIR in config
output_dir <- Sys.getenv("OUTPUT_DIR")               # results/

if (fasta_dir == "") stop("[ERROR] PATH_TO_FASTA_UNIPROT is not set. Source config.sh first.")
if (output_dir == "") stop("[ERROR] OUTPUT_DIR is not set. Source config.sh first.")

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
num_cores <- max(1L, num_cores - 1L)   # leave one core free for OS

cat("[CONFIG] PATH_TO_FASTA_UNIPROT:", fasta_dir, "\n")
cat("[CONFIG] OUTPUT_DIR:           ", output_dir, "\n")
cat("[CONFIG] Cores to use:         ", num_cores, "\n\n")

###########################################################################
# 2. Load combined UniProt FASTA (canonical + isoforms)
###########################################################################

fasta_path <- file.path(fasta_dir, "UP000005640_9606_combined.fasta")

# Friendly error if user hasn't built the combined file yet
if (!file.exists(fasta_path)) {
  stop(
    "[ERROR] Combined FASTA not found at: ", fasta_path, "\n",
    "  Build it once with:\n",
    "    cd ", fasta_dir, "\n",
    "    cat UP000005640_9606.fasta UP000005640_9606_additional.fasta \\\n",
    "        > UP000005640_9606_combined.fasta\n"
  )
}

cat("[INFO] Loading UniProt combined FASTA:", fasta_path, "\n")
fasta_aa   <- readAAStringSet(fasta_path)
fasta_seqs <- as.character(fasta_aa)
cat("[INFO] Sequences loaded:", length(fasta_seqs), "\n")

# Quality filter: remove sequences with ambiguous/stop characters
before <- length(fasta_seqs)
fasta_seqs <- fasta_seqs[!grepl("[*BJOUXZ]", fasta_seqs)]
cat("[INFO] After quality filter:", length(fasta_seqs),
    "(removed", before - length(fasta_seqs), ")\n\n")

###########################################################################
# 3. Tile all sequences into 8-11mers in parallel
#
#    Per length, parallelise across proteins (not across lengths) so each
#    worker handles one protein at a time and memory stays predictable.
#    After tiling, deduplicate on n_mer within that length -- many proteins
#    share peptides, and we only need one representative row per unique
#    peptide since downstream tools (NetMHCpan / MHCflurry) score by
#    sequence, not by protein of origin.
###########################################################################

cat("[INFO] Starting peptide tiling (", num_cores, "cores)...\n")

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export the sequence vector to workers once (avoids re-serialising per task)
clusterExport(cl, varlist = "fasta_seqs", envir = environment())

all_mers_list <- foreach(
  h         = 8:11,
  .combine  = "c",          # collect a list of data.tables
  .packages = c("data.table", "stringr")
) %dopar% {

  tile_one <- function(seq, h) {
    len <- nchar(seq)
    if (is.na(len) || len < h) return(NULL)
    starts <- seq_len(len - h + 1L)
    data.table(
      n_mer   = substring(seq, starts, starts + h - 1L),
      n_flank = substring(seq, pmax(1L, starts - 30L), starts - 1L),
      c_flank = substring(seq, starts + h, pmin(len, starts + h + 29L))
    )
  }

  chunks <- lapply(fasta_seqs, tile_one, h = h)
  dt     <- rbindlist(chunks, use.names = FALSE, fill = FALSE)

  # Deduplicate: keep first occurrence of each unique peptide
  dt <- dt[!duplicated(dt$n_mer)]

  # Drop any remaining non-standard AAs that slipped through the earlier filter
  dt <- dt[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", n_mer)]

  list(dt)   # wrap in list() so foreach .combine = "c" gives a plain list
}

stopCluster(cl)

cat("[INFO] Tiling complete.\n\n")

###########################################################################
# 4. Write output files in exactly the format Steps 13a and 14a expect
#
#    Format (matches what Step 12 writes):
#      Column 1: n_mer      -- the peptide sequence
#      Column 2: ctex_up    -- 30-char left-padded N-flank (pad char = "-")
#      Column 3: ctex_dn    -- 30-char right-padded C-flank (pad char = "-")
#      Column 4: TPM        -- blank (no TPM available at proteome level)
#
#    14a then reads these, strips the "-" padding, and expands by allele.
###########################################################################

mer_lengths <- 8:11

for (idx in seq_along(mer_lengths)) {

  h      <- mer_lengths[idx]
  lenpad <- sprintf("%02d", h)
  dt_h   <- all_mers_list[[idx]]

  n_raw <- nrow(dt_h)
  cat(sprintf("[INFO] %smer: %d unique peptides before final filters\n", lenpad, n_raw))

  # Pad flanks to exactly 30 characters (left-pad N-flank, right-pad C-flank)
  dt_h[, ctex_up := str_pad(n_flank, width = 30L, side = "left",  pad = "-")]
  dt_h[, ctex_dn := str_pad(c_flank, width = 30L, side = "right", pad = "-")]
  dt_h[, TPM     := ""]

  out_dt   <- dt_h[, .(n_mer, ctex_up, ctex_dn, TPM)]
  out_file <- file.path(output_dir,
                        paste0("2023_0812_hlathenalist_msic_", lenpad, "mers_wt.tsv"))

  fwrite(out_dt, out_file, sep = "\t", na = "NA",
         col.names = TRUE, quote = FALSE)

  cat(sprintf("[INFO] Written: %s (%d rows)\n", out_file, nrow(out_dt)))
}

###########################################################################
# 5. Done
###########################################################################

runtime <- proc.time() - start_time
cat("\n[DONE] Step 11b complete.\n")
cat(sprintf("[DONE] Total runtime: %.1f seconds (%.1f minutes)\n",
            runtime[3], runtime[3] / 60))
cat("[DONE] The four _wt.tsv files in", output_dir,
    "now reflect the full human proteome.\n")
cat("[DONE] Steps 13a and 14a will pick them up automatically -- no changes needed.\n")