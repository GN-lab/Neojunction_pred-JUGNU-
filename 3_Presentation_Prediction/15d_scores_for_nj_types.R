#!/usr/bin/env Rscript
# Step 15d: Matrix + Bed files + Figures (Tier 1 only)
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Using Tier 1 (high confidence, both tools agree) candidates
#          mapped to neojunctions in 15c:
#   1. Generate bed files (ALT tumour-specific + WT native)
#   2. Build sample x neoantigen matrix (Tier 1 only)
#   3. Generate HLA and sample summary tables
#   4. Generate figures (NJ type, frameshift, HLA, sample)
#
# Input (15c): alt_neoA_to_neoJ_map_YYYYMMDD.tsv
#              wt_neoA_to_neoJ_map_YYYYMMDD.tsv
# Input (results): Count_Neojunctions_YYYYMMDD.tsv
#                  Patient_List_Post_TumorPurity_Filter_0.60.txt
#                  2023_0812_peptide_coordinate_map_XXmers.tsv

###########################################################################
#  Step 0: Packages and config
###########################################################################

rm(list = ls(all.names = TRUE))
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

output_dir        <- Sys.getenv("OUTPUT_DIR")
directory_15      <- Sys.getenv("STEP15_OUTPUT_DIR")
directory_figures <- Sys.getenv("STEP15_FIGURES_DIR")

if (nchar(output_dir)   == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_15) == 0) directory_15 <- output_dir
if (nchar(directory_figures) == 0)
  directory_figures <- file.path(output_dir, "figures", "step15")

dir.create(directory_figures, showWarnings = FALSE, recursive = TRUE)

# Tier colours (consistent with 15b)
TIER_COLOURS <- c(
  "Tier1_HighConfidence"       = "#2dc653",
  "Tier2_MediumConfidence"     = "#f4a100",
  "Tier3_Discordant_NMPstrong" = "#e63946",
  "Tier3_Discordant_MHCstrong" = "#e07a5f"
)

# NJ type colours
TYPE_COLOURS <- c(
  "A3.gain"            = "#f94144", "A3.loss"          = "#f3722c",
  "A5.gain"            = "#f8961e", "A5.loss"          = "#f9c74f",
  "ES"                 = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
  "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS"           = "#577590"
)

###########################################################################
#  Step 1: Detect dates from existing files
###########################################################################

# 15c output date
map_scan <- list.files(directory_15,
                        pattern = "^alt_neoA_to_neoJ_map_[0-9]{8}\\.tsv$")
if (length(map_scan) == 0)
  stop("[ERROR] No alt_neoA_to_neoJ_map_YYYYMMDD.tsv in: ", directory_15,
       "\n  Has Step 15c finished?")
date_15c <- sub("^alt_neoA_to_neoJ_map_([0-9]{8})\\.tsv$", "\\1",
                map_scan[which.max(
                  file.info(file.path(directory_15, map_scan))$mtime)])

# Count NJ date
count_scan <- list.files(output_dir,
                          pattern = "^Count_Neojunctions_[0-9]{8}\\.tsv$")
date_count <- if (length(count_scan) > 0) {
  sub("^Count_Neojunctions_([0-9]{8})\\.tsv$", "\\1",
      count_scan[which.max(
        file.info(file.path(output_dir, count_scan))$mtime)])
} else { NULL }

current_date <- format(Sys.Date(), "%Y%m%d")

cat("[INFO] 15c output date:  ", date_15c,    "\n")
cat("[INFO] Count NJ date:    ", date_count,   "\n")
cat("[INFO] Output date:      ", current_date, "\n")

###########################################################################
#  Step 2: Load 15c mapped outputs
###########################################################################

cat("\n[STEP 2] Loading 15c junction maps...\n")

alt_map <- fread(file.path(directory_15,
                            paste0("alt_neoA_to_neoJ_map_", date_15c, ".tsv")),
                 na.strings = c("", "NA"), quote = "")
cat(sprintf("[INFO] ALT map rows: %d\n", nrow(alt_map)))

# Filter to Tier 1 only for matrix and bed
alt_t1 <- alt_map[concordance_tier == "Tier1_HighConfidence"]
cat(sprintf("[INFO] ALT Tier 1 rows: %d\n", nrow(alt_t1)))
cat(sprintf("[INFO] ALT Tier 1 unique peptides:   %d\n", uniqueN(alt_t1$peptide)))
cat(sprintf("[INFO] ALT Tier 1 unique junctions:  %d\n", uniqueN(alt_t1$junc.id)))

# WT map (optional)
wt_map_file <- file.path(directory_15,
                          paste0("wt_neoA_to_neoJ_map_", date_15c, ".tsv"))
wt_t1 <- if (file.exists(wt_map_file)) {
  dt <- fread(wt_map_file, na.strings = c("", "NA"), quote = "")
  dt[concordance_tier == "Tier1_HighConfidence"]
} else {
  cat("[WARN] WT map file not found -- skipping WT bed\n")
  data.table()
}
cat(sprintf("[INFO] WT Tier 1 rows: %d\n", nrow(wt_t1)))

###########################################################################
#  Step 3: Generate bed files
###########################################################################

cat("\n[STEP 3] Generating bed files...\n")

# Load coordinate maps
coord_files <- file.path(output_dir, paste0(
  "2023_0812_peptide_coordinate_map_", c("08","09","10","11"), "mers.tsv"))
coord_files_exist <- coord_files[file.exists(coord_files)]

if (length(coord_files_exist) == 0) {
  warning("[WARN] No coordinate map files found -- skipping bed generation")
} else {
  coord_map <- rbindlist(lapply(coord_files_exist, fread,
                                na.strings = c("", "NA")),
                         use.names = TRUE, fill = TRUE)

  # Standardise column names
  if ("ctex_up" %in% colnames(coord_map))
    setnames(coord_map, "ctex_up", "n_flank")
  if ("ctex_dn" %in% colnames(coord_map))
    setnames(coord_map, "ctex_dn", "c_flank")
  if ("n_mer" %in% colnames(coord_map))
    setnames(coord_map, "n_mer", "peptide")

  cat(sprintf("[INFO] Coordinate map: %d rows\n", nrow(coord_map)))

  make_bed <- function(candidates, label, outfile) {
    if (nrow(candidates) == 0) {
      warning(sprintf("[WARN] No candidates for %s bed", label))
      return(invisible(NULL))
    }

    # Join on peptide + flanks if flanks present, otherwise peptide only.
    # WT map from 15c does not carry flanks so we fall back to peptide only.
    join_by <- intersect(c("peptide","n_flank","c_flank"), colnames(candidates))
    if (length(join_by) == 0) join_by <- "peptide"
    cat(sprintf("[INFO] %s bed: joining on: %s\n", label,
                paste(join_by, collapse = "+")))

    joined  <- merge(candidates, coord_map,
                     by = join_by, all.x = TRUE,
                     allow.cartesian = TRUE)
    joined  <- joined[!is.na(enst.model)]

    if (nrow(joined) == 0) {
      # Try fallback: join on peptide only regardless
      cat(sprintf("[INFO] %s bed: retrying with peptide-only join...\n", label))
      joined <- merge(candidates[, .(peptide, allele)],
                      coord_map, by = "peptide",
                      all.x = TRUE, allow.cartesian = TRUE)
      joined <- joined[!is.na(enst.model)]
    }

    if (nrow(joined) == 0) {
      warning(sprintf("[WARN] No coordinate matches for %s bed", label))
      return(invisible(NULL))
    }

    # Collapse all HLA alleles per transcript + position + peptide
    bed <- joined[, .(
      HLA_ALLELES = paste(sort(unique(allele)), collapse = ",")
    ), by = .(
      ENST_ID  = enst.model,
      AA_START = aa_start,
      AA_END   = aa_end,
      PEPTIDE  = peptide
    )]
    setorder(bed, ENST_ID, AA_START)

    fwrite(bed, file.path(output_dir, outfile),
           sep = "\t", col.names = TRUE, quote = FALSE)
    cat(sprintf("[INFO] %s bed: %d rows | %d peptides | %d transcripts -> %s\n",
                label, nrow(bed), uniqueN(bed$PEPTIDE),
                uniqueN(bed$ENST_ID), outfile))
    invisible(bed)
  }

  # ALT Tier 1 bed -- tumour-specific immunogenic peptides
  alt_bed <- make_bed(
    alt_t1, "ALT Tier1",
    paste0("immunogenic_peptides_", current_date, ".bed")
  )

  # WT Tier 1 bed -- native immunopeptidome reference
  # The coordinate map only contains ALT peptides so we cannot use make_bed.
  # Instead: load complete_list_all_mers to get aa.seq.wt and enst.model,
  # then find each WT peptide's position within aa.seq.wt to get AA coords.
  if (nrow(wt_t1) > 0) {
    cat("[INFO] Building WT bed from aa.seq.wt in complete_list_all_mers...\n")

    mers_file <- file.path(output_dir, "2023_0812_complete_list_all_mers.tsv")
    if (!file.exists(mers_file)) {
      warning("[WARN] complete_list_all_mers.tsv not found -- skipping WT bed")
    } else {
      mers_wt <- fread(mers_file, na.strings = c("", "NA"), quote = "",
                       select = c("junc.id", "enst.model", "aa.seq.wt"))
      mers_wt <- mers_wt[!is.na(aa.seq.wt) & nchar(aa.seq.wt) > 0]
      # Keep longest aa.seq.wt per junc.id
      mers_wt[, seq_len := nchar(aa.seq.wt)]
      setorder(mers_wt, junc.id, -seq_len)
      mers_wt <- mers_wt[, .SD[1], by = junc.id]

      # For each unique WT peptide x junction, find position in aa.seq.wt
      wt_unique <- unique(wt_t1[, .(peptide, junc.id, allele)])
      wt_joined <- merge(wt_unique, mers_wt, by = "junc.id",
                         all.x = TRUE, allow.cartesian = TRUE)
      wt_joined <- wt_joined[!is.na(aa.seq.wt)]

      # Find AA start position of peptide in wt sequence
      wt_joined[, aa_start := {
        pos <- regexpr(peptide, aa.seq.wt, fixed = TRUE)
        fifelse(pos > 0, as.integer(pos), NA_integer_)
      }]
      wt_joined <- wt_joined[!is.na(aa_start)]
      wt_joined[, aa_end := aa_start + nchar(peptide) - 1L]

      if (nrow(wt_joined) == 0) {
        warning("[WARN] No WT peptides found in aa.seq.wt -- skipping WT bed")
      } else {
        # Collapse alleles per transcript + position + peptide
        wt_bed <- wt_joined[, .(
          HLA_ALLELES = paste(sort(unique(allele)), collapse = ",")
        ), by = .(
          ENST_ID  = enst.model,
          AA_START = aa_start,
          AA_END   = aa_end,
          PEPTIDE  = peptide
        )]
        setorder(wt_bed, ENST_ID, AA_START)

        wt_bed_file <- paste0("wt_native_immunogenic_peptides_", current_date, ".bed")
        fwrite(wt_bed, file.path(output_dir, wt_bed_file),
               sep = "	", col.names = TRUE, quote = FALSE)
        cat(sprintf("[INFO] WT bed: %d rows | %d peptides | %d transcripts -> %s\n",
                    nrow(wt_bed), uniqueN(wt_bed$PEPTIDE),
                    uniqueN(wt_bed$ENST_ID), wt_bed_file))
      }
    }
  }
}

###########################################################################
#  Step 4: Build sample x neoantigen matrix (Tier 1 only)
###########################################################################

cat("\n[STEP 4] Building sample x neoantigen matrix...\n")

if (is.null(date_count)) {
  warning("[WARN] Count_Neojunctions file not found -- skipping matrix")
} else {
  count_file <- file.path(output_dir,
                           paste0("Count_Neojunctions_", date_count, ".tsv"))
  count_neo  <- fread(count_file, na.strings = c("", "NA"))
  cat(sprintf("[INFO] Count NJ matrix: %d rows\n", nrow(count_neo)))

  # Sample names from count matrix columns
  sample_names <- sub("\\.$", "", setdiff(colnames(count_neo), "junc.id"))
  cat(sprintf("[INFO] Samples in count matrix: %d\n", length(sample_names)))

  # Build neo_id for Tier 1
  alt_t1[, neo_id := paste(peptide, junc.id, allele, sep = "|")]

  # Melt count matrix to long
  count_long <- melt(count_neo, id.vars = "junc.id",
                     variable.name = "sample", value.name = "count")
  count_long[, sample := sub("\\.$", "", as.character(sample))]
  count_long <- count_long[count > 0]

  # Join Tier 1 neoantigens with sample counts via junction ID
  score_long <- merge(count_long,
                      alt_t1[, .(junc.id, neo_id, combined_score)],
                      by = "junc.id", allow.cartesian = TRUE)

  # Wide matrix: neo_id x sample
  score_wide <- dcast(score_long, neo_id ~ sample,
                      value.var = "combined_score",
                      fun.aggregate = max, fill = NA_real_)

  # Metadata columns
  meta_cols <- intersect(
    c("neo_id","peptide","junc.id","allele","concordance_tier",
      "combined_score","netmhcpan_EL_score","netmhcpan_EL_rank",
      "mhcflurry_affinity","mhcflurry_presentation_score",
      "symbol","type","fs"),
    colnames(alt_t1)
  )
  # Best row per neo_id for metadata
  alt_t1_meta <- alt_t1[, .SD[which.min(netmhcpan_EL_rank)], by = neo_id]

  neo_matrix <- merge(alt_t1_meta[, meta_cols, with = FALSE],
                      score_wide, by = "neo_id", all.x = TRUE)
  setDT(neo_matrix)

  # Add any missing sample columns
  missing_s <- setdiff(sample_names, colnames(neo_matrix))
  if (length(missing_s) > 0) neo_matrix[, (missing_s) := NA_real_]

  # Summary stats per row
  sample_cols_in_matrix <- intersect(sample_names, colnames(neo_matrix))
  neo_matrix[, samples_with_neo := rowSums(
    !is.na(.SD)), .SDcols = sample_cols_in_matrix]
  neo_matrix[, max_sample_score  := do.call(
    pmax, c(.SD, na.rm = TRUE)), .SDcols = sample_cols_in_matrix]

  cat(sprintf("[INFO] Matrix: %d neoantigens x %d samples\n",
              nrow(neo_matrix), length(sample_cols_in_matrix)))

  # Long format
  neo_long <- melt(neo_matrix,
                   id.vars       = meta_cols,
                   measure.vars  = sample_cols_in_matrix,
                   variable.name = "sample_id",
                   value.name    = "sample_score")
  neo_long <- neo_long[!is.na(sample_score)]

  # HLA summary
  hla_summary <- alt_t1[, .(
    n_tier1_peptides  = .N,
    avg_combined_score = mean(combined_score, na.rm = TRUE),
    max_combined_score = max(combined_score,  na.rm = TRUE),
    avg_EL_rank        = mean(netmhcpan_EL_rank, na.rm = TRUE),
    avg_mhc_affinity   = mean(mhcflurry_affinity, na.rm = TRUE),
    n_unique_junctions = uniqueN(junc.id)
  ), by = allele][order(-max_combined_score)]

  # Sample summary
  sample_summary <- data.table(
    sample_id         = sample_cols_in_matrix,
    n_tier1_neoantigens = colSums(!is.na(
      neo_matrix[, .SD, .SDcols = sample_cols_in_matrix])),
    max_score = sapply(sample_cols_in_matrix,
                       function(s) {
                         vals <- neo_matrix[[s]]
                         if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
                       }),
    avg_score = sapply(sample_cols_in_matrix,
                       function(s) {
                         vals <- neo_matrix[[s]]
                         if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
                       })
  )

  # Join patient purity if available
  patient_file <- file.path(output_dir,
                             "Patient_List_Post_TumorPurity_Filter_0.60.txt")
  if (file.exists(patient_file)) {
    patient_list <- fread(patient_file)
    sample_summary <- merge(sample_summary, patient_list,
                            by = "sample_id", all.x = TRUE)
  }

  # Write matrix outputs
  fwrite(neo_matrix,     file.path(output_dir,
    paste0("top_neoantigens_matrix_",     current_date, ".tsv")), sep="\t")
  fwrite(neo_long,       file.path(output_dir,
    paste0("top_neoantigens_long_",       current_date, ".tsv")), sep="\t")
  fwrite(hla_summary,    file.path(output_dir,
    paste0("hla_alleles_summary_",        current_date, ".tsv")), sep="\t")
  fwrite(sample_summary, file.path(output_dir,
    paste0("sample_neoantigen_summary_",  current_date, ".tsv")), sep="\t")

  cat(sprintf("[INFO] Matrix written: %d rows\n",     nrow(neo_matrix)))
  cat(sprintf("[INFO] Long written:   %d rows\n",     nrow(neo_long)))
  cat(sprintf("[INFO] HLA summary:    %d alleles\n",  nrow(hla_summary)))
  cat(sprintf("[INFO] Sample summary: %d samples\n",  nrow(sample_summary)))
}

###########################################################################
#  Step 5: Figures
###########################################################################

cat("\n[STEP 5] Generating figures...\n")

# Dynamic HLA colour palette
hla_alleles  <- sort(unique(alt_t1$allele))
n_hla        <- length(hla_alleles)
hla_cols     <- if (n_hla <= 12) {
  setNames(brewer.pal(max(3, n_hla), "Set3"), hla_alleles)
} else {
  setNames(colorRampPalette(brewer.pal(12, "Set3"))(n_hla), hla_alleles)
}

# Helper: save PDF
save_pdf <- function(p, filename, w = 10, h = 6) {
  f <- file.path(directory_figures, filename)
  ggsave(f, plot = p, width = w, height = h, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", f))
}

# --- Figure 1: Jitter -- frameshift vs in-frame ---
if ("fs" %in% colnames(alt_t1) && nrow(alt_t1) > 0) {
  df_fs <- alt_t1[, .(allele, combined_score, fs)]
  df_fs[, fs_label := fifelse(fs == "fs", "Frame-shift", "In-frame")]

  p_fs_jitter <- ggplot(df_fs,
                         aes(x = fs_label, y = combined_score, colour = allele)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_colour_manual(values = hla_cols) +
    theme_bw() +
    labs(x = "", y = "Combined score (NMP + MHC)",
         title = "Tier 1: Frame-shift vs In-frame neoantigens") +
    theme(axis.text.x  = element_text(size = 12, face = "bold"),
          legend.text  = element_text(size = 8)) +
    guides(colour = guide_legend(ncol = 2, title = "HLA Allele"))

  save_pdf(p_fs_jitter,
           paste0("figure_15d_tier1_fs_if_jitter_", current_date, ".pdf"))

  p_fs_box <- ggplot(df_fs,
                      aes(x = allele, y = combined_score, fill = fs_label)) +
    geom_boxplot(position = position_dodge()) +
    scale_fill_manual(values = c("Frame-shift" = "#006e90",
                                  "In-frame"    = "#f18f01")) +
    theme_bw() +
    labs(x = "HLA Allele", y = "Combined score", fill = "Type",
         title = "Tier 1: Combined score by allele and frameshift status") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

  save_pdf(p_fs_box,
           paste0("figure_15d_tier1_fs_if_boxplot_", current_date, ".pdf"),
           w = max(10, n_hla * 0.4))
}

# --- Figure 2: Jitter -- by NJ splice type ---
if ("type" %in% colnames(alt_t1) && nrow(alt_t1) > 0) {
  type_counts <- alt_t1[, .N, by = type]
  df_type     <- alt_t1[, .(allele, combined_score, type)]
  df_type     <- merge(df_type, type_counts, by = "type")
  df_type[, type_label := paste0(type, " (n=", N, ")")]

  p_type_jitter <- ggplot(df_type,
                           aes(x = type_label, y = combined_score,
                               colour = allele)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_colour_manual(values = hla_cols) +
    theme_bw() +
    labs(x = "", y = "Combined score",
         title = "Tier 1: Neoantigens by splice junction type") +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 9),
          legend.text  = element_text(size = 8)) +
    guides(colour = guide_legend(ncol = 2, title = "HLA Allele"))

  save_pdf(p_type_jitter,
           paste0("figure_15d_tier1_splice_type_jitter_", current_date, ".pdf"),
           w = 13)
}

# --- Figure 3: HLA allele bar chart ---
if (exists("hla_summary") && nrow(hla_summary) > 0) {
  p_hla <- ggplot(hla_summary,
                   aes(x = reorder(allele, -n_tier1_peptides),
                       y = n_tier1_peptides,
                       fill = allele)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n_tier1_peptides), vjust = -0.4, size = 3) +
    scale_fill_manual(values = hla_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw() +
    labs(x = "HLA Allele", y = "Tier 1 neoantigens",
         title = "Tier 1 neoantigens per HLA allele") +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  save_pdf(p_hla,
           paste0("figure_15d_tier1_hla_summary_", current_date, ".pdf"),
           w = max(10, n_hla * 0.35))
}

# --- Figure 4: Sample neoantigen counts ---
if (exists("sample_summary") && nrow(sample_summary) > 0) {
  fill_var <- if ("purity" %in% colnames(sample_summary)) "purity" else NULL

  # Build sample plot -- colour by purity if available, plain fill otherwise
  if (!is.null(fill_var)) {
    p_sample <- ggplot(sample_summary,
                        aes(x = reorder(sample_id, -n_tier1_neoantigens),
                            y = n_tier1_neoantigens, fill = purity)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis_c() +
      theme_bw() +
      labs(x = "Sample", y = "Tier 1 neoantigens",
           title = "Tier 1 neoantigens per sample") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
  } else {
    p_sample <- ggplot(sample_summary,
                        aes(x = reorder(sample_id, -n_tier1_neoantigens),
                            y = n_tier1_neoantigens)) +
      geom_bar(stat = "identity", fill = "#4d908e") +
      theme_bw() +
      labs(x = "Sample", y = "Tier 1 neoantigens",
           title = "Tier 1 neoantigens per sample") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
  }

  save_pdf(p_sample,
           paste0("figure_15d_tier1_sample_summary_", current_date, ".pdf"),
           w = 14)
}

###########################################################################
#  Summary
###########################################################################

cat("\n=== 15d Summary ===\n")
cat(sprintf("  ALT Tier 1 peptides:       %d\n", nrow(alt_t1)))
cat(sprintf("  ALT Tier 1 unique junctions: %d\n", uniqueN(alt_t1$junc.id)))
cat(sprintf("  ALT Tier 1 unique alleles:   %d\n", uniqueN(alt_t1$allele)))
if (exists("neo_matrix"))
  cat(sprintf("  Matrix dimensions:         %d neoantigens x %d samples\n",
              nrow(neo_matrix), length(sample_cols_in_matrix)))
cat(sprintf("  Figures written to: %s\n", directory_figures))
cat("\n[DONE] Step 15d complete.\n")