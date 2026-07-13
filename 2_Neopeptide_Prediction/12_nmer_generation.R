#!/usr/bin/env Rscript
# Title: "Step 12: Generate AA n-mers from Predicted AA Sequences"
# December 04, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: The previous step generated all of the predicted amino acid sequences that would be generated
#          from the neojunctions. This step will generate all of the predicted n-mers (from 8-mers to 11-mers),
#          and return only the n-mers that cannot be generated in the respective WT amino acid sequence.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))
library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
library(stringr)  # For str_pad in flank padding (if needed elsewhere; otherwise optional)
library(tibble)

# Establish Directories ---------------------------------------------------
directory_00 <- Sys.getenv("META_FILES_PATH")
directory_01 <- Sys.getenv("STEP01_OUTPUT_DIR")
directory_11 <- Sys.getenv("STEP11_OUTPUT_DIR")
directory_12 <- Sys.getenv("OUTPUT_DIR")

# Load Files --------------------------------------------------------------

# From External: TPM Values (Hartwig transcript TPM matrix)
setwd(directory_00)
filename_tpm = "transcript_tpm_matrix.tsv"
# Specify col_types to suppress parsing warnings (adjust based on your data; this assumes 'tx' is character and others are double)
tpm = read_tsv(filename_tpm, na = c("", "NA"), col_types = cols(tx = col_character(), .default = col_double()))

# Inspect for initial duplicates (for debugging; remove once fixed)
duplicated_cols <- colnames(tpm)[duplicated(colnames(tpm))]
if (length(duplicated_cols) > 0) {
  print("Initial duplicates found:")
  print(duplicated_cols)
}

# Ensure unique column names to avoid filter error
colnames(tpm) <- make.unique(colnames(tpm))

# Rename 'tx' to 'enst' (original expects 'enst')—check if 'tx' exists first
if ("tx" %in% colnames(tpm)) {
  tpm <- tpm %>% rename(enst = tx)
} else {
  warning("'tx' column not found—check your TPM file headers.")
}

# From TPM header: Get ACTUAL sample names
tpm_colnames <- colnames(tpm)[-(1:2)]  # Skip 'tx' and 'gene_id' columns
message("[INFO] Found ", length(tpm_colnames), " samples in TPM matrix")
message("[INFO] First 5 samples: ", paste(head(tpm_colnames, 5), collapse = ", "))

# Create tcga dataframe from TPM samples
tcga <- data.frame(case = tpm_colnames)
message("[INFO] Using ", nrow(tcga), " samples from TPM file")
  
# From Step 11: List of Predicted AA Sequences
setwd(directory_11)
filename_aa = list.files(pattern = "^Res_AA_Prediction_Confirmed_[0-9]{8}\\.tsv$") %>% sort(decreasing = TRUE) %>% .[1]
aa = read_tsv(filename_aa, na = c("", "NA"))

###########################################################################
#  Step 1: Edit the TPM dataframe -----------------------------------------
###########################################################################

# The following steps are taken from Step 3 (03_tpm_filter_10.R)
# For greater detail on each step, refer to that step.

colnames(tpm) = substr(colnames(tpm), 1, 15)
colnames(tpm) <- make.unique(colnames(tpm))
print(length(unique(colnames(tpm)))) 
# Inspect post-shortening duplicates (for debugging; remove once fixed)
duplicated_cols_after <- colnames(tpm)[duplicated(colnames(tpm))]
if (length(duplicated_cols_after) > 0) {
  print("Duplicates after shortening:")
  print(duplicated_cols_after)
} else {
  print("No duplicates after fix.")
}

tpm.edit = tpm %>% dplyr::filter(substr(enst, 1, 4) == "ENST")
tpm.edit = tpm.edit %>% dplyr::select(enst, colnames(tpm)[is.element(colnames(tpm), tcga$case)])
tpm.edit = tpm.edit %>% gather("case", "log2tpm", 2:ncol(tpm.edit))
tpm.edit = tpm.edit %>% mutate(tpm = round((2^as.numeric(log2tpm)) - 0.001, 4))
tpm.edit = tpm.edit %>% mutate(tpm = ifelse(is.na(tpm), 0, tpm))
tpm.edit = tpm.edit %>% mutate(tpm = ifelse(tpm < 0, 0, tpm))
tpm.edit = tpm.edit %>% dplyr::select( - log2tpm)
tpm.edit = tpm.edit %>% spread(case, tpm)
tpm.edit = tpm.edit %>%  mutate(enst.1 = substr(enst, 1, 15))

# Select out only the median TPM for all of the genes
tpm.edit = tpm.edit %>% 
  mutate(TPM = apply(tpm.edit[, -1], 1, median)) %>% 
  dplyr::select(enst.1, TPM)


###########################################################################
#  Step 2: Edit and Filter AA Sequences -----------------------------------
###########################################################################

# Remove NA sequences or any sequences that do not start with a start codon (M)

aa_valid = aa %>% 
  dplyr::select(junc.id, symbol, enst, enst.model, aa.change, aa.seq.wt, aa.seq.alt, ln.wt, ln.alt, ln.diff, sc, type) %>% 
  dplyr::filter(!is.na(aa.seq.wt) & !is.na(aa.seq.alt)) %>%              # Remove any sequences with NA in either the WT or Mut variants
  dplyr::filter(aa.seq.wt != aa.seq.alt) %>%                             # Ensure that the WT and Mut variants are not the same
  dplyr::filter(grepl("^M", aa.seq.wt) & grepl("^M", aa.seq.alt)) %>%    # Ensure that both AA sequences start with M
  mutate(aa.seq.alt = gsub("\\*.*", "*", aa.seq.alt))                    # Remove the first stop codon (*) and any amino acids following it


# aa %>% dim() %>% print()       --> 467  17
# aa_valid %>% dim() %>% print() --> 390  12

# Add the TPM values as the final column
aa_valid = aa_valid %>% 
  left_join(tpm.edit, by = c("enst" = "enst.1"))

# Turn all NA values under TPM into blanks
aa_valid$TPM[is.na(aa_valid$TPM)] = ""

# Remove duplicate aa.alts
aa_valid.1 = aa_valid[!duplicated(aa_valid$aa.seq.alt), ]

# Optional smoke-test toggle (NEW, same env var as Step 11): set
# TEST_N_JUNCTIONS before running to limit this run to just the first N
# AA predictions from Step 11, e.g. `export TEST_N_JUNCTIONS=5`. Unset
# (the default) means a normal full run.
test_n <- Sys.getenv("TEST_N_JUNCTIONS", unset = "")
if (test_n != "") {
  test_n <- as.integer(test_n)
  cat("[TEST MODE] TEST_N_JUNCTIONS is set -- limiting to the first", test_n,
      "of", nrow(aa_valid.1), "AA predictions.\n")
  aa_valid.1 <- aa_valid.1 %>% dplyr::slice(1:min(test_n, nrow(aa_valid.1)))
}
 
# ###########################################################################
# #  Step 3: Edit the AA sequences to only Include 8-11 mers ----------------
# ###########################################################################
# # Generate list of n-mers for MUTANT peptides (aa.seq.alt)
# setwd(directory_12)
# for (h in 8:11){
#   summary_nmer_h = NULL
#   for (i in 1:nrow(aa_valid.1)){
#     aa_i = aa_valid.1 %>% dplyr::slice(i)
#     aa_sequence = aa_valid.1 %>% dplyr::slice(i) %>% pull(aa.seq.alt)
#     
#     # Iterate through the entire length of the peptide
#     for (j in 1:nchar(aa_sequence)){
#       print(paste(h, round((i/nrow(aa_valid.1)*100), digits = 2), round((j/nchar(aa_sequence)*100), digits = 2)))
#       
#       # Generate the n-mer and their flanks
#       nmer_h = substr(aa_sequence, start = j, stop = j+h-1)
#       flank_10n = substr(aa_sequence, start = j-10, stop = j-1)
#       flank_10c = substr(aa_sequence, start = j+h, stop = j+h+10-1)
#       flank_30n = substr(aa_sequence, start = j-30, stop = j-1)
#       flank_30c = substr(aa_sequence, start = j+h, stop = j+h+30-1)
#       
#       append_row = aa_i %>% 
#         mutate(n_flank = flank_30n) %>%
#         mutate(n_mer = nmer_h) %>%
#         mutate(c_flank = flank_30c)
#       
#       if(nchar(nmer_h) == h){summary_nmer_h = rbind(summary_nmer_h, append_row)}
#     }
#   }
#   if (h == 8){nmer_8_mut = summary_nmer_h ; write_tsv(nmer_8_mut, "2023_0802_all_iterations_alt_list_08mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 9){nmer_9_mut = summary_nmer_h ; write_tsv(nmer_9_mut, "2023_0802_all_iterations_alt_list_09mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 10){nmer_10_mut = summary_nmer_h ; write_tsv(nmer_10_mut, "2023_0802_all_iterations_alt_list_10mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 11){nmer_11_mut = summary_nmer_h ; write_tsv(nmer_11_mut, "2023_0802_all_iterations_alt_list_11mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
# }

# Use SLURM cores if provided; leave one free
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1))
num_cores <- max(1, num_cores - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

generate_nmers <- function(aa_i, aa_sequence, h) {
  seq_len <- nchar(aa_sequence)
  if (is.na(seq_len) || seq_len < h) return(NULL)
  starts <- 1:(seq_len - h + 1)
  nmer_h    <- substring(aa_sequence, starts, starts + h - 1)
  flank_30n <- substring(aa_sequence, pmax(starts - 30, 1), starts - 1)
  flank_30c <- substring(aa_sequence, starts + h, pmin(starts + h + 29, seq_len))
  dt <- data.table(aa_i)[rep(1, length(nmer_h)), ]
  # aa_start/aa_end (NEW): the n-mer's 1-based position within aa_sequence --
  # this was already being computed as `starts` above, just never exported.
  # This IS the peptide coordinate (e.g. for a SOPRANO bed file): no exon or
  # genomic information needed, just this index pair plus enst.model.
  dt[, `:=`(n_flank = flank_30n, n_mer = nmer_h, c_flank = flank_30c,
            aa_start = starts, aa_end = starts + h - 1)]
  dt
}

setwd(directory_12)
for (type in c("alt", "wt")) {
  for (h in 8:11) {
    # Return a list() per task, and use rbindlist on the driver
    chunks <- foreach(
      i = 1:nrow(aa_valid.1),
      .combine = "c",
      .packages = c("data.table")
    ) %dopar% {
      aa_i <- aa_valid.1[i, ]
      aa_sequence <- aa_i[[paste0("aa.seq.", type)]]
      out <- generate_nmers(aa_i, aa_sequence, h)
      list(out)  # IMPORTANT: list() so the combiner gets a plain list
    }
    # Stack safely on the driver
    summary_nmer_h <- rbindlist(chunks, use.names = TRUE, fill = FALSE)
    # Ensure a data.frame/tibble for write_tsv
    summary_nmer_h <- as_tibble(summary_nmer_h)
    # Write exact filenames as in the original
    if (h == 8)  write_tsv(summary_nmer_h, paste0("2023_0802_all_iterations_", ifelse(type=="alt","alt","wt"), "_list_08mers.tsv"), na="NA", col_names=TRUE, escape="double")
    if (h == 9)  write_tsv(summary_nmer_h, paste0("2023_0802_all_iterations_", ifelse(type=="alt","alt","wt"), "_list_09mers.tsv"), na="NA", col_names=TRUE, escape="double")
    if (h == 10) write_tsv(summary_nmer_h, paste0("2023_0802_all_iterations_", ifelse(type=="alt","alt","wt"), "_list_10mers.tsv"), na="NA", col_names=TRUE, escape="double")
    if (h == 11) write_tsv(summary_nmer_h, paste0("2023_0802_all_iterations_", ifelse(type=="alt","alt","wt"), "_list_11mers.tsv"), na="NA", col_names=TRUE, escape="double")
  }
}
stopCluster(cl)

###########################################################################
#  Step 4: Remove aa.seq.alt n-mers plus flanks found in aa.seq.wt --------
###########################################################################
setwd(directory_12)
nmer08_mut = read_tsv("2023_0802_all_iterations_alt_list_08mers.tsv", col_names = T)
nmer09_mut = read_tsv("2023_0802_all_iterations_alt_list_09mers.tsv", col_names = T)
nmer10_mut = read_tsv("2023_0802_all_iterations_alt_list_10mers.tsv", col_names = T)
nmer11_mut = read_tsv("2023_0802_all_iterations_alt_list_11mers.tsv", col_names = T)

nmer08_wt = read_tsv("2023_0802_all_iterations_wt_list_08mers.tsv", col_names = T)
nmer09_wt = read_tsv("2023_0802_all_iterations_wt_list_09mers.tsv", col_names = T)
nmer10_wt = read_tsv("2023_0802_all_iterations_wt_list_10mers.tsv", col_names = T)
nmer11_wt = read_tsv("2023_0802_all_iterations_wt_list_11mers.tsv", col_names = T)

nmer_08_neo <- anti_join(nmer08_mut, nmer08_wt, by = c("n_flank", "n_mer", "c_flank")) # 50563 neoantigens
nmer_09_neo <- anti_join(nmer09_mut, nmer09_wt, by = c("n_flank", "n_mer", "c_flank")) # 50801 neoantigens
nmer_10_neo <- anti_join(nmer10_mut, nmer10_wt, by = c("n_flank", "n_mer", "c_flank")) # 51039 neoantigens
nmer_11_neo <- anti_join(nmer11_mut, nmer11_wt, by = c("n_flank", "n_mer", "c_flank")) # 51275 neoantigens

nmer_all = rbind((rbind(rbind(nmer_08_neo, nmer_09_neo), nmer_10_neo)), nmer_11_neo)

###########################################################################
#  Step 4b (NEW): Flag n-mers spanning an exon-exon boundary --------------
###########################################################################
# Uses the ALT exon -> AA coordinate map written by Step 11
# (Exon_CDS_AA_Map_*.tsv, only rows that passed CDS-reconstruction
# validation) to check whether each retained n-mer's [aa_start, aa_end]
# window is drawn from a single exon or spans a junction. This is the
# highest-risk case for a coordinate mistake, so it's flagged explicitly
# rather than left implicit.

setwd(directory_11)
exon_map_candidates <- list.files(pattern = "^Exon_CDS_AA_Map_[0-9]{8}\\.tsv$")

if (length(exon_map_candidates) > 0) {
  exon_map_file <- sort(exon_map_candidates, decreasing = TRUE)[1]
  exon_map_dt <- fread(exon_map_file, na.strings = c("", "NA"))
  data.table::setnames(exon_map_dt, c("aa_start", "aa_end"), c("exon_aa_start", "exon_aa_end"))

  nmer_dt <- as.data.table(nmer_all)
  nmer_dt[, row_id := .I]

  flag_result <- tryCatch({
    data.table::setkey(exon_map_dt, enst.model, exon_aa_start, exon_aa_end)

    overlaps <- data.table::foverlaps(
      nmer_dt[, .(row_id, enst.model, aa_start, aa_end)],
      exon_map_dt[, .(enst.model, exon_aa_start, exon_aa_end, exon_idx)],
      by.x = c("enst.model", "aa_start", "aa_end"),
      by.y = c("enst.model", "exon_aa_start", "exon_aa_end"),
      type = "any",
      nomatch = NULL
    )

    span_summary <- overlaps[, .(n_exons_spanned = data.table::uniqueN(exon_idx)), by = row_id]

    merged <- merge(nmer_dt, span_summary, by = "row_id", all.x = TRUE)
    merged[is.na(n_exons_spanned), n_exons_spanned := 0L]
    merged[, spans_exon_boundary := n_exons_spanned > 1]
    merged[, row_id := NULL]
    merged
  }, error = function(e) {
    cat("[WARNING] Exon-boundary flagging failed:", conditionMessage(e), "\n")
    nmer_dt[, row_id := NULL]
    nmer_dt[, `:=`(n_exons_spanned = NA_integer_, spans_exon_boundary = NA)]
    nmer_dt
  })

  nmer_all <- as_tibble(flag_result)
  cat("Exon-boundary flagging: ", sum(nmer_all$spans_exon_boundary, na.rm = TRUE),
      " of ", nrow(nmer_all), " retained n-mers span an exon-exon junction (0 unflaggable: ",
      sum(is.na(nmer_all$spans_exon_boundary)), ")\n", sep = "")
} else {
  cat("[WARNING] No Exon_CDS_AA_Map file found in", directory_11,
      "-- run the updated Step 11 first. Skipping exon-boundary flagging.\n")
  nmer_all <- nmer_all %>% dplyr::mutate(n_exons_spanned = NA_integer_, spans_exon_boundary = NA)
}

setwd(directory_12)

###########################################################################
#  Step 5: Export Complete Files ------------------------------------------
###########################################################################
setwd(directory_12)
# Export the Complete n-mer files
filename_8mer  = "2023_0812_complete_list_08mers.tsv"
filename_9mer  = "2023_0812_complete_list_09mers.tsv"
filename_10mer = "2023_0812_complete_list_10mers.tsv"
filename_11mer = "2023_0812_complete_list_11mers.tsv"
filename_all   = "2023_0812_complete_list_all_mers.tsv"

write_tsv(nmer_08_neo, filename_8mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_09_neo, filename_9mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_10_neo, filename_10mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_11_neo, filename_11mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_all, filename_all, na = "NA", col_names = T, quote_escape = "double") ;

###########################################################################
#  Step 6: Filter for cancer-specific n-mers ------------------------------
###########################################################################

library(Biostrings)  # For readAAStringSet

setwd(directory_12)
nmer_8  <- read_tsv(filename_8mer, na = c("", "NA"))
nmer_9  <- read_tsv(filename_9mer, na = c("", "NA"))
nmer_10 <- read_tsv(filename_10mer, na = c("", "NA"))
nmer_11 <- read_tsv(filename_11mer, na = c("", "NA"))

# Load canonical + isoform FASTAs (NEW: both combined into one exclusion set)
# UP000005640_9606.fasta          = canonical sequences (one per gene)
# UP000005640_9606_additional.fasta = known alternative isoforms
# Filtering against both means any peptide from ANY annotated human protein
# sequence (canonical or isoform) is removed before MHCflurry runs --
# no point predicting binding for something that already exists in a normal
# human protein. Saves compute in 14a/14b and keeps the neoantigen list clean.
fasta_dir <- Sys.getenv("PATH_TO_FASTA_UNIPROT")
setwd(fasta_dir)

filename_fasta_canonical <- "UP000005640_9606.fasta"
filename_fasta_isoform   <- "UP000005640_9606_additional.fasta"

fastaFile_canonical <- readAAStringSet(filename_fasta_canonical)
cat("Canonical FASTA:", length(fastaFile_canonical), "sequences\n")

if (file.exists(filename_fasta_isoform)) {
  fastaFile_isoform <- readAAStringSet(filename_fasta_isoform)
  cat("Isoform FASTA:", length(fastaFile_isoform), "sequences\n")
  # Combine both into one list for the k-mer extraction loop below
  fasta_list <- as.list(as.character(c(fastaFile_canonical, fastaFile_isoform)))
  cat("Combined:", length(fasta_list), "sequences total (canonical + isoform)\n")
} else {
  warning("Isoform FASTA not found at: ", filename_fasta_isoform,
          " -- filtering against canonical sequences only. Download with:\n",
          "wget -P ", fasta_dir, " https://ftp.uniprot.org/pub/databases/uniprot/",
          "current_release/knowledgebase/reference_proteomes/Eukaryota/",
          "UP000005640/UP000005640_9606_additional.fasta.gz && ",
          "gunzip ", fasta_dir, "/UP000005640_9606_additional.fasta.gz")
  fasta_list <- as.list(as.character(fastaFile_canonical))
}
rm(fastaFile_canonical)
if (exists("fastaFile_isoform")) rm(fastaFile_isoform)
gc()

# Function to extract all h-mers from a single sequence
extract_mers <- function(seq, h) {
  len <- nchar(seq)
  if (len < h) return(character(0))
  starts <- 1:(len - h + 1)
  substring(seq, starts, starts + h - 1)
}

# Set up cluster and export
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 1))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary objects and functions
clusterExport(cl, varlist = c("fasta_list", "extract_mers"), envir = environment())

# Now your foreach loop (adapted from previous suggestion)
for (h in 8:11) {
  # Select df_h (as before)
  if (h == 8) df_h <- nmer_8
  if (h == 9) df_h <- nmer_9
  if (h == 10) df_h <- nmer_10
  if (h == 11) df_h <- nmer_11
  
  # Parallel extraction of h-mers
  mers_list <- foreach(i = 1:length(fasta_list), .combine = 'c', .packages = "base") %dopar% {
    extract_mers(fasta_list[[i]], h)  # Use the exported list
  }
  
  # Unique proteome h-mers
  all_mers <- unique(mers_list)
  
  # Filter
  keep_idx <- ! (df_h$n_mer %in% all_mers)
  df_valid <- df_h[keep_idx, ]
  
  # Assign (as before)
  if (h == 8) df_valid_08 <- df_valid
  if (h == 9) df_valid_09 <- df_valid
  if (h == 10) df_valid_10 <- df_valid
  if (h == 11) df_valid_11 <- df_valid
}

# Stop cluster
stopCluster(cl)

setwd(directory_12)
write_tsv(df_valid_08, "2023_0812_cancer_specific_08mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_09, "2023_0812_cancer_specific_09mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_10, "2023_0812_cancer_specific_10mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_11, "2023_0812_cancer_specific_11mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;

# NOTE: After changing to readAA rather than readDNA, the total changed from 21489 --> 21787

setwd(directory_12)
df_valid_08 = read_tsv("2023_0812_cancer_specific_08mers.tsv", col_names = T); df_valid_08$c_flank <- ifelse(is.na(df_valid_08$c_flank), "", df_valid_08$c_flank); df_valid_08$n_flank <- ifelse(is.na(df_valid_08$n_flank), "", df_valid_08$n_flank); df_valid_08$n_flank <- gsub("\\*", "", df_valid_08$n_flank); df_valid_08$c_flank <- gsub("\\*", "", df_valid_08$c_flank); df_valid_08 <- df_valid_08[!grepl("\\*", df_valid_08$n_mer), ]
df_valid_09 = read_tsv("2023_0812_cancer_specific_09mers.tsv", col_names = T); df_valid_09$c_flank <- ifelse(is.na(df_valid_09$c_flank), "", df_valid_09$c_flank); df_valid_09$n_flank <- ifelse(is.na(df_valid_09$n_flank), "", df_valid_09$n_flank); df_valid_09$n_flank <- gsub("\\*", "", df_valid_09$n_flank); df_valid_09$c_flank <- gsub("\\*", "", df_valid_09$c_flank); df_valid_09 <- df_valid_09[!grepl("\\*", df_valid_09$n_mer), ]
df_valid_10 = read_tsv("2023_0812_cancer_specific_10mers.tsv", col_names = T); df_valid_10$c_flank <- ifelse(is.na(df_valid_10$c_flank), "", df_valid_10$c_flank); df_valid_10$n_flank <- ifelse(is.na(df_valid_10$n_flank), "", df_valid_10$n_flank); df_valid_10$n_flank <- gsub("\\*", "", df_valid_10$n_flank); df_valid_10$c_flank <- gsub("\\*", "", df_valid_10$c_flank); df_valid_10 <- df_valid_10[!grepl("\\*", df_valid_10$n_mer), ]
df_valid_11 = read_tsv("2023_0812_cancer_specific_11mers.tsv", col_names = T); df_valid_11$c_flank <- ifelse(is.na(df_valid_11$c_flank), "", df_valid_11$c_flank); df_valid_11$n_flank <- ifelse(is.na(df_valid_11$n_flank), "", df_valid_11$n_flank); df_valid_11$n_flank <- gsub("\\*", "", df_valid_11$n_flank); df_valid_11$c_flank <- gsub("\\*", "", df_valid_11$c_flank); df_valid_11 <- df_valid_11[!grepl("\\*", df_valid_11$n_mer), ]

# Function to pad strings with hyphens to reach a total length of 30 characters
pad_to_30_front <- function(s) {
  str_pad(s, width = 30, side = "left", pad = "-")
}

pad_to_30_back <- function(s) {
  str_pad(s, width = 30, side = "right", pad = "-")
}

df_valid_08$ctex_up <- sapply(df_valid_08$n_flank, pad_to_30_front); df_valid_08$ctex_dn <- sapply(df_valid_08$c_flank, pad_to_30_back)
df_valid_09$ctex_up <- sapply(df_valid_09$n_flank, pad_to_30_front); df_valid_09$ctex_dn <- sapply(df_valid_09$c_flank, pad_to_30_back)
df_valid_10$ctex_up <- sapply(df_valid_10$n_flank, pad_to_30_front); df_valid_10$ctex_dn <- sapply(df_valid_10$c_flank, pad_to_30_back)
df_valid_11$ctex_up <- sapply(df_valid_11$n_flank, pad_to_30_front); df_valid_11$ctex_dn <- sapply(df_valid_11$c_flank, pad_to_30_back)

# Generate the HLAthena Compatible n-mer files for each of the algorithms (MSi, MSiC, and MSiCE)
# Files for the MSiCE algorithm
MSiCE_8  = df_valid_08 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM") 
MSiCE_9  = df_valid_09 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")
MSiCE_10 = df_valid_10 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")
MSiCE_11 = df_valid_11 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")

# Files for the MSiC algorithm
MSiC_8  = df_valid_08 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_9  = df_valid_09 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_10 = df_valid_10 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_11 = df_valid_11 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")

# Files for the MSiC algorithm
MSi_8 = df_valid_08 %>% dplyr::select("n_mer") %>%
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "") 
MSi_9 = df_valid_09 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")
MSi_10 = df_valid_10 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")
MSi_11 = df_valid_11 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")

# Export the HLAthena Compatible n-mer files
filename_MSi_8  = "2023_0812_hlathenalist_msi_08mers.tsv"
filename_MSi_9  = "2023_0812_hlathenalist_msi_09mers.tsv"
filename_MSi_10 = "2023_0812_hlathenalist_msi_10mers.tsv"
filename_MSi_11 = "2023_0812_hlathenalist_msi_11mers.tsv"

filename_MSiC_8  = "2023_0812_hlathenalist_msic_08mers.tsv"
filename_MSiC_9  = "2023_0812_hlathenalist_msic_09mers.tsv"
filename_MSiC_10 = "2023_0812_hlathenalist_msic_10mers.tsv"
filename_MSiC_11 = "2023_0812_hlathenalist_msic_11mers.tsv"

filename_MSiCE_8  = "2023_0812_hlathenalist_msice_08mers.tsv"
filename_MSiCE_9  = "2023_0812_hlathenalist_msice_09mers.tsv"
filename_MSiCE_10 = "2023_0812_hlathenalist_msice_10mers.tsv"
filename_MSiCE_11 = "2023_0812_hlathenalist_msice_11mers.tsv"

write_tsv(MSi_8, filename_MSi_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSi_9, filename_MSi_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSi_10, filename_MSi_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSi_11, filename_MSi_11, na = "NA", col_names = T, quote_escape = "double")

write_tsv(MSiC_8, filename_MSiC_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiC_9, filename_MSiC_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiC_10, filename_MSiC_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSiC_11, filename_MSiC_11, na = "NA", col_names = T, quote_escape = "double")

write_tsv(MSiCE_8, filename_MSiCE_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiCE_9, filename_MSiCE_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiCE_10, filename_MSiCE_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSiCE_11, filename_MSiCE_11, na = "NA", col_names = T, quote_escape = "double")

###########################################################################
#  Step 6b (NEW): Peptide -> coordinate map for ALT n-mers ----------------
###########################################################################
# The MSi/MSiC/MSiCE files above intentionally keep only what MHCflurry
# needs (n_mer, ctex_up, ctex_dn, TPM) -- everything else gets dropped by
# that select(). This companion file keeps the SAME peptide+flank identity
# (so it can be joined back onto MHCflurry's own output later, since
# MHCflurry echoes back whatever peptide/flank values it was given) PLUS
# junc.id, enst.model, aa_start, aa_end -- needed downstream (Step 15e) to
# build a SOPRANO-style peptide coordinate bed file.

coord_map_08 <- df_valid_08 %>% dplyr::select(junc.id, enst.model, n_mer, ctex_up, ctex_dn, aa_start, aa_end)
coord_map_09 <- df_valid_09 %>% dplyr::select(junc.id, enst.model, n_mer, ctex_up, ctex_dn, aa_start, aa_end)
coord_map_10 <- df_valid_10 %>% dplyr::select(junc.id, enst.model, n_mer, ctex_up, ctex_dn, aa_start, aa_end)
coord_map_11 <- df_valid_11 %>% dplyr::select(junc.id, enst.model, n_mer, ctex_up, ctex_dn, aa_start, aa_end)

write_tsv(coord_map_08, "2023_0812_peptide_coordinate_map_08mers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(coord_map_09, "2023_0812_peptide_coordinate_map_09mers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(coord_map_10, "2023_0812_peptide_coordinate_map_10mers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(coord_map_11, "2023_0812_peptide_coordinate_map_11mers.tsv", na = "NA", col_names = T, quote_escape = "double")

cat("Wrote peptide coordinate map files (junc.id, enst.model, aa_start, aa_end) -- used by Step 15e to build the SOPRANO bed file.\n")

print("step 12 has run successfullly, please find your result files in the output directory")

###########################################################################
#  Step 7 (NEW): Export WT n-mers as a "self-peptide" reference set -------
###########################################################################
# Purpose: Step 6 above applies a UniProt-proteome "cancer-specific" filter
# to the ALT-derived n-mers only -- that filter is intentionally NOT applied
# here, because WT-derived n-mers are *expected* to occur in the normal
# proteome; that's the whole point of using them as the self-peptide
# reference. These WT n-mers are formatted (same MSiC layout: n_mer,
# ctex_up, ctex_dn, TPM) so Step 14a/14b can run them through MHCflurry
# exactly like the ALT side. The resulting WT presentation/affinity
# predictions are then used in Step 15e to remove any ALT "neoantigen"
# (peptide + HLA allele) that also arises from the WT/self sequence --
# something the exact n_flank/n_mer/c_flank match in Step 4 above can miss
# when flanking context differs slightly but the core presented peptide is
# effectively the same self-epitope.

setwd(directory_12)

# Re-use the WT n-mer tables already loaded in Step 4 above (nmer08_wt ... nmer11_wt)
wt_raw_list <- list("08" = nmer08_wt, "09" = nmer09_wt, "10" = nmer10_wt, "11" = nmer11_wt)

clean_wt_nmer <- function(df) {
  df$c_flank <- ifelse(is.na(df$c_flank), "", df$c_flank)
  df$n_flank <- ifelse(is.na(df$n_flank), "", df$n_flank)
  df$n_flank <- gsub("\\*", "", df$n_flank)
  df$c_flank <- gsub("\\*", "", df$c_flank)
  df <- df[!grepl("\\*", df$n_mer), ]
  df
}

for (h in names(wt_raw_list)) {
  df_wt_h <- clean_wt_nmer(wt_raw_list[[h]])
  df_wt_h$ctex_up <- sapply(df_wt_h$n_flank, pad_to_30_front)
  df_wt_h$ctex_dn <- sapply(df_wt_h$c_flank, pad_to_30_back)

  MSiC_wt_h <- df_wt_h %>%
    dplyr::select("n_mer", "ctex_up", "ctex_dn") %>%
    mutate(TPM = "")

  filename_wt_h <- paste0("2023_0812_hlathenalist_msic_", h, "mers_wt.tsv")
  write_tsv(MSiC_wt_h, filename_wt_h, na = "NA", col_names = T, quote_escape = "double")
  cat("Wrote WT reference n-mer file:", filename_wt_h, "(", nrow(MSiC_wt_h), "rows )\n")
}

print("step 12 (WT reference n-mer export) has run successfully -- see *_mers_wt.tsv files, to be consumed by Step 14a/14b (WT side) and filtered against in Step 15e")
