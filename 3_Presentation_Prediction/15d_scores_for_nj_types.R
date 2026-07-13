#!/usr/bin/env Rscript
# Step 15d: Generate Figure 5i - Top validated NJs and their corresponding neoantigens
# November 25, 2025 | Gaurav Raichand | The Institute of Cancer Research
# UPDATED: December 3, 2025 | Fixed for new MHCflurry format and dynamic HLA alleles

# Purpose: Generate visualizations for neoantigen scores and identify top HLA alleles/neoantigens
# Joins verified with perfect input overlap; substring matching for mapping metadata; no NAs/empty outputs
# Fixes: Deduplicate inputs to avoid many-to-many joins and OOM; memory-efficient data.table usage
# Patient list path hardcoded to "results/Patient_List_Post_TumorPurity_Filter_0.60.txt" based on your ls output


###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################


rm(list = ls(all.names = TRUE))


library(tidyverse)
library(ggsci)
library(data.table)
library(RColorBrewer)
library(gridExtra)
library(parallel)
library(doParallel)
library(foreach)
library(Cairo)
library(stringr)  # For substring matching


# Directories with fallbacks
directory_15 <- Sys.getenv("STEP15_OUTPUT_DIR", getwd())
directory_figures <- Sys.getenv("OUTPUT_DIR", getwd())
directory_step10 <- Sys.getenv("STEP10_OUTPUT_DIR", getwd())


current_date <- format(Sys.Date(), "%Y%m%d")


# Set working directory to match config.sh OUTPUT_DIR
setwd(Sys.getenv("OUTPUT_DIR", "results"))
cat("Working in:", getwd(), "\n")


# Define cache files
CACHE_DIR <- "cache"
dir.create(CACHE_DIR, showWarnings = FALSE)


CACHE_ALL_MAP <- file.path(CACHE_DIR, "df_all_map_cached.rds")
CACHE_NA_NJ_MAP <- file.path(CACHE_DIR, "df_na_nj_map_cached.rds")


###########################################################################
#  Helper Functions -------------------------------------------------------
###########################################################################


# Helper: pick most recent file matching a pattern
latest_file <- function(pattern, dir = NULL) {
  if (is.null(dir)) {
    dir <- "."
  }
  files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files match pattern ", pattern, " in ", dir)
  }
  files[which.max(file.info(files)$mtime)]
}


###########################################################################
#  Step 0.5: Load or Compute Main Data (COMPUTATIONALLY INTENSIVE) --------
###########################################################################


if (file.exists(CACHE_ALL_MAP)) {
  cat("Loading cached df_all_map...\n")
  df_all_map <- readRDS(CACHE_ALL_MAP)
} else {
  cat("Computing df_all_map (this will take time)...\n")
  
  # Load MHCflurry output files as data.table
  mf_patterns <- c("^mhcflurry_08mer_top_\\d{8}\\.tsv$",
                   "^mhcflurry_09mer_top_\\d{8}\\.tsv$",
                   "^mhcflurry_10mer_top_\\d{8}\\.tsv$",
                   "^mhcflurry_11mer_top_\\d{8}\\.tsv$")
  mf_files <- lapply(mf_patterns, latest_file)
  cat("MHCflurry TOP files:\n", paste(mf_files, collapse="\n"), "\n")
  df_all_map <- rbindlist(lapply(mf_files, fread), use.names = TRUE, fill = TRUE)
  if ("mhcflurry_affinity" %in% names(df_all_map))
    setnames(df_all_map, "mhcflurry_affinity", "mhcflurry_binding_affinity")

  ###########################################################################
  #  WT FILTER + BED FILES (added) -----------------------------------------
  ###########################################################################

  # 1. WT filter -- remove any ALT peptide+allele found in WT predictions.
  #    No score threshold: presence alone disqualifies (the sequence exists
  #    in the normal self repertoire regardless of how strongly it binds).
  wt_patterns <- c("^08mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
                   "^09mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
                   "^10mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
                   "^11mers_flank_mhcflurry_wt_[0-9_]+\\.csv$")
  wt_files <- tryCatch(lapply(wt_patterns, latest_file), error = function(e) NULL)
  if (!is.null(wt_files)) {
    cat("Loading WT predictions for filter...\n")
    wt_all <- rbindlist(lapply(wt_files, fread, na.strings = c("", "NA")),
                        use.names = TRUE, fill = TRUE)
    if ("mhcflurry_affinity" %in% names(wt_all))
      setnames(wt_all, "mhcflurry_affinity", "mhcflurry_binding_affinity")
    wt_set <- unique(wt_all[, .(peptide, allele)])
    rm(wt_all); gc()
    n_before <- nrow(df_all_map)
    df_all_map <- df_all_map[!wt_set, on = .(peptide, allele)]
    cat("WT filter: removed", n_before - nrow(df_all_map),
        "| retained", nrow(df_all_map), "tumor-specific rows\n")
    rm(wt_set)
  } else {
    cat("[WARNING] WT prediction files not found -- skipping WT filter\n")
  }

  # 2. Bed files -- written now while we have the data, before anything else.
  #    Format: ENST_ID / AA_START / AA_END / PEPTIDE / HLA_ALLELES / N_ALLELES
  #    Source: raw 14b predictions filtered at <=500nM -- NOT 14c top files.
  #    14c keeps only ONE best allele per peptide. Raw 14b has ALL alleles.
  #    A peptide like RVHYKGTGR can bind 11 alleles at <=500nM -- we want
  #    all of them comma-separated in HLA_ALLELES, not just the best one.
  AFFINITY_THRESHOLD <- 500  # nM

  coord_map_files <- c(
    "2023_0812_peptide_coordinate_map_08mers.tsv",
    "2023_0812_peptide_coordinate_map_09mers.tsv",
    "2023_0812_peptide_coordinate_map_10mers.tsv",
    "2023_0812_peptide_coordinate_map_11mers.tsv"
  )

  make_bed <- function(raw_dt, coord_map, label) {
    # raw_dt: peptide, n_flank, c_flank, allele -- all alleles at <=500nM
    joined <- merge(raw_dt, coord_map,
                    by.x = c("peptide","n_flank","c_flank"),
                    by.y = c("n_mer","n_flank","c_flank"),
                    all.x = TRUE, allow.cartesian = TRUE)
    # Collapse all qualifying alleles per (transcript, position, peptide)
    bed <- joined[!is.na(enst.model), .(
      HLA_ALLELES = paste(sort(unique(allele)), collapse = ","),
      N_ALLELES   = uniqueN(allele)
    ), by = .(ENST_ID  = enst.model,
              AA_START = aa_start,
              AA_END   = aa_end,
              PEPTIDE  = peptide)]
    setorder(bed, ENST_ID, AA_START)
    cat(label, "bed:", nrow(bed), "rows |",
        uniqueN(bed$PEPTIDE), "unique peptides |",
        uniqueN(bed$ENST_ID), "transcripts |",
        round(mean(bed$N_ALLELES), 2), "alleles/peptide avg\n")
    bed
  }

  load_raw_14b <- function(patterns, threshold, filter_peps = NULL) {
    files <- tryCatch(lapply(patterns, latest_file), error = function(e) NULL)
    if (is.null(files)) return(NULL)
    dt <- rbindlist(lapply(files, fread, na.strings = c("","NA")),
                    use.names = TRUE, fill = TRUE)
    if ("mhcflurry_affinity" %in% names(dt))
      setnames(dt, "mhcflurry_affinity", "mhcflurry_binding_affinity")
    dt <- dt[!is.na(mhcflurry_binding_affinity) &
               mhcflurry_binding_affinity <= threshold]
    if (!is.null(filter_peps)) dt <- dt[peptide %in% filter_peps]
    dt
  }

  if (all(file.exists(coord_map_files))) {
    coord_map <- rbindlist(lapply(coord_map_files, fread, na.strings = c("", "NA")),
                           use.names = TRUE, fill = TRUE)
    setnames(coord_map, c("ctex_up","ctex_dn"), c("n_flank","c_flank"))

    # ALT bed -- raw 14b ALT predictions, tumor-specific peptides only
    tumor_peps <- unique(df_all_map$peptide)
    alt_raw <- load_raw_14b(
      c("^08mers_flank_mhcflurry_[0-9_]+\\.csv$",
        "^09mers_flank_mhcflurry_[0-9_]+\\.csv$",
        "^10mers_flank_mhcflurry_[0-9_]+\\.csv$",
        "^11mers_flank_mhcflurry_[0-9_]+\\.csv$"),
      AFFINITY_THRESHOLD, filter_peps = tumor_peps
    )
    if (!is.null(alt_raw)) {
      alt_bed <- make_bed(alt_raw[, .(peptide, n_flank, c_flank, allele)],
                          coord_map, "ALT filtered")
      fwrite(alt_bed, paste0("immunopeptidome_alt_filtered_", current_date, ".bed"),
             sep = "\t", col.names = TRUE, quote = FALSE)
      cat("Wrote immunopeptidome_alt_filtered_", current_date, ".bed\n", sep = "")
      rm(alt_raw, alt_bed)
    }
    rm(tumor_peps)

    # WT bed -- raw 14b WT predictions, all alleles <=500nM
    wt_nmer_files <- c(
      "2023_0802_all_iterations_wt_list_08mers.tsv",
      "2023_0802_all_iterations_wt_list_09mers.tsv",
      "2023_0802_all_iterations_wt_list_10mers.tsv",
      "2023_0802_all_iterations_wt_list_11mers.tsv"
    )
    wt_raw <- load_raw_14b(
      c("^08mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
        "^09mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
        "^10mers_flank_mhcflurry_wt_[0-9_]+\\.csv$",
        "^11mers_flank_mhcflurry_wt_[0-9_]+\\.csv$"),
      AFFINITY_THRESHOLD, filter_peps = NULL
    )
    if (!is.null(wt_raw) && all(file.exists(wt_nmer_files))) {
      wt_nmer_dt <- rbindlist(lapply(wt_nmer_files, fread, na.strings = c("","NA")),
                              use.names = TRUE, fill = TRUE)
      wt_coord <- wt_nmer_dt[, .(n_mer, n_flank, c_flank, enst.model, aa_start, aa_end)]
      wt_coord[, n_flank := stringr::str_pad(ifelse(is.na(n_flank),"",n_flank),30,"left", "-")]
      wt_coord[, c_flank := stringr::str_pad(ifelse(is.na(c_flank),"",c_flank),30,"right","-")]
      rm(wt_nmer_dt); gc()
      wt_bed <- make_bed(wt_raw[, .(peptide, n_flank, c_flank, allele)],
                         wt_coord, "WT native")
      fwrite(wt_bed, paste0("immunopeptidome_wt_native_", current_date, ".bed"),
             sep = "\t", col.names = TRUE, quote = FALSE)
      cat("Wrote immunopeptidome_wt_native_", current_date, ".bed\n", sep = "")
      rm(wt_raw, wt_coord, wt_bed); gc()
    } else {
      cat("[WARNING] WT raw 14b or n-mer files not found -- skipping WT native bed\n")
    }
    rm(coord_map); gc()
  } else {
    cat("[WARNING] Coordinate map files not found -- skipping bed files\n")
  }

  ###########################################################################
  #  END OF NEW ADDITIONS -- original script continues below ---------------
  ###########################################################################

  # Load MHCflurry input files as data.table and deduplicate
  input_patterns <- c("^08mer_mhcflurry_input_.*\\.csv$",
                      "^09mer_mhcflurry_input_.*\\.csv$",
                      "^10mer_mhcflurry_input_.*\\.csv$",
                      "^11mer_mhcflurry_input_.*\\.csv$")
  input_files <- lapply(input_patterns, latest_file)
  cat("Input files:\n", paste(input_files, collapse="\n"), "\n")
  df_input <- rbindlist(lapply(input_files, fread), use.names = TRUE, fill = TRUE)
  df_input <- unique(df_input, by = c("peptide", "allele"))
  
  # Join with inputs
  setkey(df_all_map, peptide, allele)
  setkey(df_input, peptide, allele)
  df_all_map <- df_all_map[df_input, nomatch = NA]
  
  # Fast junction mapping via coordinate map (replaces slow substring match)
  mapping_file <- latest_file("^.*complete_list_all_mers\\.tsv$")
  cat("Mapping file:", mapping_file, "\n")
  df_mapping <- fread(mapping_file)

  if ("junc.id" %in% names(df_all_map)) {
    # Already have junc.id from coordinate map -- just look up type and fs
    map_lookup <- unique(df_mapping[, .(
      junc.id,
      type,
      fs = fifelse(grepl("shift|fs", aa.change) | (ln.diff %% 3 != 0), "fs", "in-frame")
    )])
    df_all_map <- merge(df_all_map, map_lookup, by = "junc.id", all.x = TRUE)
    cat("Junction type/fs mapped via junc.id join\n")
  } else {
    # Build junc.id lookup from coordinate map first
    coord_map_files2 <- c(
      "2023_0812_peptide_coordinate_map_08mers.tsv",
      "2023_0812_peptide_coordinate_map_09mers.tsv",
      "2023_0812_peptide_coordinate_map_10mers.tsv",
      "2023_0812_peptide_coordinate_map_11mers.tsv"
    )
    if (all(file.exists(coord_map_files2))) {
      cm2 <- rbindlist(lapply(coord_map_files2, fread, na.strings = c("","NA")),
                       use.names = TRUE, fill = TRUE)
      setnames(cm2, c("ctex_up","ctex_dn"), c("n_flank","c_flank"))
      # One junc.id per peptide+flank combination
      cm2 <- cm2[, .(junc.id = junc.id[1]), by = .(n_mer, n_flank, c_flank)]
      df_all_map <- merge(df_all_map, cm2,
                          by.x = c("peptide","n_flank","c_flank"),
                          by.y = c("n_mer","n_flank","c_flank"),
                          all.x = TRUE)
      map_lookup <- unique(df_mapping[, .(
        junc.id,
        type,
        fs = fifelse(grepl("shift|fs", aa.change) | (ln.diff %% 3 != 0), "fs", "in-frame")
      )])
      df_all_map <- merge(df_all_map, map_lookup, by = "junc.id", all.x = TRUE)
      cat("Attached junc.id to", sum(!is.na(df_all_map$junc.id)),
          "of", nrow(df_all_map), "rows via coordinate map\n")
      rm(cm2); gc()
    } else {
      cat("[WARNING] Coordinate map files not found -- junc.id/type/fs will be Unknown\n")
      df_all_map[, c("junc.id","type","fs") := .("Unknown","OTHERS","in-frame")]
    }
  }
  
  # FIX: Calculate percentile from binding affinity (since column doesn't exist in new files)
  cat("Calculating binding affinity percentiles...\n")
  df_all_map[, mhcflurry_affinity_percentile := rank(mhcflurry_binding_affinity) / .N * 100]
  
  # FIX: Normalize binding affinity for score calculation (lower affinity = better)
  df_all_map[, binding_norm := 1 - (mhcflurry_binding_affinity / max(mhcflurry_binding_affinity, na.rm = TRUE))]
  df_all_map[, score_average := (binding_norm + mhcflurry_presentation_score) / 2]
  
  # FIX: Use calculated percentile for classification
  df_all_map[, shared := fifelse(mhcflurry_affinity_percentile <= 10, "Top 10%tile in MF", "Other")]
  
  df_all_map[, hla_allele := allele]
  df_all_map[, type := fifelse(type == "Unknown", "OTHERS", type)]
  df_all_map[, fs := fifelse(fs == "Unknown", "in-frame", fs)]
  
  # Save to cache
  saveRDS(df_all_map, CACHE_ALL_MAP)
  cat("Saved df_all_map to cache\n")
}


if (file.exists(CACHE_NA_NJ_MAP)) {
  cat("Loading cached df_na_nj_map...\n")
  df_na_nj_map <- readRDS(CACHE_NA_NJ_MAP)
} else {
  cat("Computing df_na_nj_map...\n")
  
  # Filter for top
  df_na_nj_map <- df_all_map[shared == "Top 10%tile in MF"]
  if (nrow(df_na_nj_map) == 0) stop("No top peptides. Adjust criteria.")
  
  # Load neojunction details (use latest_file for flexibility)
  psr_neo <- fread(latest_file("^PSR_Neojunctions_\\d{8}\\.tsv$"))
  count_neo_detail <- fread(latest_file("^Count_Table_Retained_and_Passed_Junctions_\\d{8}\\.tsv$"))
  
  # Join with data.table
  setkey(df_na_nj_map, junc.id)
  setkey(psr_neo, junc.id)
  setkey(count_neo_detail, junc.id)
  df_na_nj_map <- df_na_nj_map[psr_neo][count_neo_detail]
  
  # Save to cache
  saveRDS(df_na_nj_map, CACHE_NA_NJ_MAP)
  cat("Saved df_na_nj_map to cache\n")
}


# Load patient list (this is quick)
patient_list_file <- "Patient_List_Post_TumorPurity_Filter_0.60.txt"
if (!file.exists(patient_list_file)) stop("Patient list not found.")
patient_list <- fread(patient_list_file)


###########################################################################
#  Step 1. Plot FS NJ's neoantigen scores for the top neoantigens ---------
###########################################################################


df_fs <- df_na_nj_map[fs == "fs", .(hla_allele, score_average)]
df_if <- df_na_nj_map[fs == "in-frame", .(hla_allele, score_average)]


df_combined <- rbind(
  df_fs[, type := paste0("Frame-shift (n=", .N, ")")],
  df_if[, type := paste0("In-frame (n=", .N, ")")]
)


# FIX: Generate dynamic color palette for HLA alleles
hla_alleles <- unique(df_combined$hla_allele)
n_hla <- length(hla_alleles)
cat("Found", n_hla, "HLA alleles:", paste(hla_alleles, collapse=", "), "\n")

# Create color palette (using Set3 for distinct colors)
if (n_hla <= 12) {
  hla_colors <- brewer.pal(max(3, n_hla), "Set3")
} else {
  hla_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_hla)
}
names(hla_colors) <- hla_alleles


p_jitter_fs_if <- ggplot(df_combined, aes(x = type, y = score_average, color = hla_allele)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.5, "cm")) +
  labs(y = "Average presentation score") +
  guides(color = guide_legend(ncol = 2, title = "HLA Allele"))


setwd(directory_figures)
ggsave(paste0("figure_5i_fs_if_jitter_", current_date, ".pdf"), plot = p_jitter_fs_if, width = 10, height = 6)
ggsave(paste0("figure_5i_fs_if_jitter_", current_date, ".png"), plot = p_jitter_fs_if, width = 10, height = 6)


df_all <- rbind(
  df_fs[, type := "Frame-shift"],
  df_if[, type := "In-frame"]
)


p_box_fs_if <- ggplot(df_all, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("Frame-shift" = "#006e90", "In-frame" = "#f18f01")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "bottom")


ggsave(paste0("figure_5i_fs_if_boxplot_", current_date, ".pdf"), plot = p_box_fs_if, 
       width = min(50, max(8, n_hla * 0.5)), height = 6, limitsize = FALSE)
ggsave(paste0("figure_5i_fs_if_boxplot_", current_date, ".png"), plot = p_box_fs_if, 
       width = min(50, max(8, n_hla * 0.5)), height = 6, limitsize = FALSE)

###########################################################################
#  Step 2. Plot FS NJ's neoantigen scores for ALL neoantigens ------------
###########################################################################


hist <- df_all_map[, .(hla_allele, score_average, fs)]
hist[, fs := fifelse(fs == "fs", "Frame-shift", "In-frame")]
hist[, score_log2 := log2(score_average + 0.001)]


p_density_fs_if <- ggplot(hist, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  scale_fill_manual(values = c("Frame-shift" = "#006e90", "In-frame" = "#f18f01")) +
  labs(x = "log2(Average presentation score)", y = "Density", fill = "Type") +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(paste0("figure_5i_fs_if_density_all_", current_date, ".pdf"), plot = p_density_fs_if, width = 8, height = 5)
ggsave(paste0("figure_5i_fs_if_density_all_", current_date, ".png"), plot = p_density_fs_if, width = 8, height = 5)


###########################################################################
#  Step 3. Plot TYPE NJ's neoantigen scores for the top neoantigens -------
###########################################################################


df_a3loss <- df_na_nj_map[type == "A3.loss", .(hla_allele, score_average)][, type := paste0("A3 loss (n=", .N, ")")]
df_a3gain <- df_na_nj_map[type == "A3.gain", .(hla_allele, score_average)][, type := paste0("A3 gain (n=", .N, ")")]
df_a5loss <- df_na_nj_map[type == "A5.loss", .(hla_allele, score_average)][, type := paste0("A5 loss (n=", .N, ")")]
df_a5gain <- df_na_nj_map[type == "A5.gain", .(hla_allele, score_average)][, type := paste0("A5 gain (n=", .N, ")")]
df_juncin <- df_na_nj_map[type == "JUNC.WITHIN.EXON", .(hla_allele, score_average)][, type := paste0("JWE (n=", .N, ")")]
df_juncex <- df_na_nj_map[type == "JUNC.WITHIN.INTRON", .(hla_allele, score_average)][, type := paste0("JWI (n=", .N, ")")]
df_exskip <- df_na_nj_map[type == "ES", .(hla_allele, score_average)][, type := paste0("ES (n=", .N, ")")]
df_others <- df_na_nj_map[type == "OTHERS", .(hla_allele, score_average)][, type := paste0("OTHERS (n=", .N, ")")]


df_combined_types <- rbind(df_a3loss, df_a3gain, df_a5loss, df_a5gain, df_juncin, df_juncex, df_exskip, df_others, fill = TRUE)


p_jitter_types <- ggplot(df_combined_types, aes(x = type, y = score_average, color = hla_allele)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold")) +
  labs(y = "Average presentation score") +
  guides(color = guide_legend(ncol = 2, title = "HLA Allele"))


ggsave(paste0("figure_5i_splice_types_jitter_", current_date, ".pdf"), plot = p_jitter_types, width = 12, height = 6)
ggsave(paste0("figure_5i_splice_types_jitter_", current_date, ".png"), plot = p_jitter_types, width = 12, height = 6)


df_all_types <- rbind(
  df_a3loss[, type := "A3.loss"],
  df_a3gain[, type := "A3.gain"],
  df_a5loss[, type := "A5.loss"],
  df_a5gain[, type := "A5.gain"],
  df_juncin[, type := "JUNC.WITHIN.EXON"],
  df_juncex[, type := "JUNC.WITHIN.INTRON"],
  df_exskip[, type := "ES"],
  df_others[, type := "OTHERS"],
  fill = TRUE
)


p_box_types <- ggplot(df_all_types, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                               "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                               "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                               "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, title = "Junction Type"))

width_box <- max(10, n_hla * 0.6)
ggsave(paste0("figure_5i_splice_types_boxplot_", current_date, ".pdf"), plot = p_box_types, 
       width = width_box, height = 6, limitsize = FALSE)
ggsave(paste0("figure_5i_splice_types_boxplot_", current_date, ".png"), plot = p_box_types, 
       width = width_box, height = 6, limitsize = FALSE)
       
###########################################################################
#  Step 4. Plot TYPE NJ's neoantigen scores for ALL neoantigens -----------
###########################################################################

hist_all <- df_all_map[, .(hla_allele, score_average, type)]
hist_all[, score_log2 := log2(score_average + 0.001)]

# Get all unique HLA alleles
hla_alleles_all <- unique(hist_all$hla_allele)
n_hla_all <- length(hla_alleles_all)
cat("Found", n_hla_all, "HLA alleles\n")

# STRATEGY: Select top N alleles by neoantigen count
TOP_N_ALLELES <- 15

# Count neoantigens per allele and select top N
allele_counts <- hist_all[, .N, by = hla_allele][order(-N)]
top_alleles <- allele_counts[1:min(TOP_N_ALLELES, nrow(allele_counts)), hla_allele]

cat("Selected top", length(top_alleles), "HLA alleles by neoantigen count\n")

# Filter for top alleles
hist_top <- hist_all[hla_allele %in% top_alleles]

# Create individual density plots for each top allele
plot_list <- lapply(top_alleles, function(allele) {
  hist_subset <- hist_top[hla_allele == allele]
  n_neo <- nrow(hist_subset)
  
  ggplot(hist_subset, aes(x = score_log2, fill = type)) +
    geom_density(alpha = 0.5) +
    labs(x = paste0("log2(Average presentation score) (", allele, ", n=", n_neo, ")"), 
         y = "Density",
         title = allele) +
    scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                                 "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                                 "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                                 "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 9),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 2, title = "Junction Type"))
})

# FIX: Use pdf() device instead of ggsave()
setwd(directory_figures)

pdf(paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".pdf"), 
    width = 10, height = 6)

p_arranged <- grid.arrange(grobs = plot_list, ncol = 1, heights = rep(1, length(plot_list)))
print(p_arranged)

dev.off()

cat("✅ Saved top", TOP_N_ALLELES, "alleles to PDF\n")

CairoPNG(paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".png"), 
         width = 10, height = 6, units = "in", res = 150)

p_arranged <- grid.arrange(grobs = plot_list, ncol = 1, heights = rep(1, length(plot_list)))
print(p_arranged)

dev.off()

cat("✅ Saved to:", paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".pdf\n"))
cat("✅ Saved to:", paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".png\n"))

###########################################################################
#  OPTIONAL: Save all 105 alleles as paginated PDF (supplementary) --------
###########################################################################

cat("Generating supplementary archive (all", n_hla_all, "alleles)...\n")

# Use pdf() device for paginated output - NO SIZE LIMIT!
pdf(paste0("figure_5i_splice_types_density_ALL_archive_", current_date, ".pdf"), 
    width = 10, height = 6)

for(i in seq_along(hla_alleles_all)) {
  allele <- hla_alleles_all[i]
  hist_subset <- hist_all[hla_allele == allele]
  n_neo <- nrow(hist_subset)
  
  p <- ggplot(hist_subset, aes(x = score_log2, fill = type)) +
    geom_density(alpha = 0.5) +
    labs(title = paste0(allele, " (n=", n_neo, ")"),
         x = "log2(Average presentation score)", 
         y = "Density") +
    scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                                 "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                                 "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                                 "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 9),
          legend.key.size = unit(0.3, "cm")) +
    guides(fill = guide_legend(ncol = 2, title = "Junction Type"))
  
  print(p)
  
  if(i %% 10 == 0) cat("Processed", i, "/", n_hla_all, "alleles\n")
}

dev.off()

cat("✅ Saved all", n_hla_all, "alleles to:", 
    paste0("figure_5i_splice_types_density_ALL_archive_", current_date, ".pdf\n"))

###########################################################################
#  Summary ----------------------------------------------------------------
###########################################################################

cat("\n=== FIGURE 5i GENERATION COMPLETE ===\n")
cat("Main figure (top", TOP_N_ALLELES, "alleles):\n")
cat("  PDF:", paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".pdf\n"))
cat("  PNG:", paste0("figure_5i_splice_types_density_top", TOP_N_ALLELES, "_", current_date, ".png\n"))
cat("\nSupplementary archive (all", n_hla_all, "alleles):\n")
cat("  PDF:", paste0("figure_5i_splice_types_density_ALL_archive_", current_date, ".pdf\n"))
cat("  Pages:", n_hla_all, "\n")

###########################################################################
#  Step 5. Identify Top HLA Alleles and Neoantigens (MATRIX FORMAT) -------
###########################################################################


cat("Using cached df_all_map for Step 5...\n")


# Select top 10% by score_average
top_neo <- df_all_map[order(-score_average)][1:round(0.1 * .N)]
cat("Selected", nrow(top_neo), "top neoantigens\n")


# Load Count_Neojunctions
count_neo_file <- latest_file("^Count_Neojunctions_\\d{8}\\.tsv$")
count_neo <- fread(count_neo_file)
cat("Loaded Count_Neojunctions from:", count_neo_file, "\n")


# Get sample names (remove trailing dots)
sample_names <- names(count_neo)[-1]  # Remove junc.id column
sample_names <- sub("\\.$", "", sample_names)
cat("Found", length(sample_names), "samples in Count_Neojunctions\n")


# Create a unique neoantigen identifier
top_neo[, neo_id := paste(peptide, junc.id, hla_allele, sep = "|")]


# Initialize empty matrix: rows = neoantigens, columns = samples
neo_matrix <- as.data.table(matrix(NA,
                                   nrow = nrow(top_neo),
                                   ncol = length(sample_names),
                                   dimnames = list(NULL, sample_names)))


# Add neoantigen metadata as first columns
neo_matrix[, c("neo_id", "peptide", "junc.id", "hla_allele", "score_average", "type", "fs", "shared", "mhcflurry_affinity_percentile", "mhcflurry_binding_affinity", "mhcflurry_presentation_score") :=
             list(top_neo$neo_id, top_neo$peptide, top_neo$junc.id, top_neo$hla_allele,
                  top_neo$score_average, top_neo$type, top_neo$fs, top_neo$shared,
                  top_neo$mhcflurry_affinity_percentile, top_neo$mhcflurry_binding_affinity, top_neo$mhcflurry_presentation_score)]


# Reorder columns to have metadata first, then samples
setcolorder(neo_matrix, c("neo_id", "peptide", "junc.id", "hla_allele", "score_average", 
                          "mhcflurry_affinity_percentile", "mhcflurry_binding_affinity", "mhcflurry_presentation_score",
                          "type", "fs", "shared", sample_names))


# Fast per-sample matrix using dcast (replaces slow nested for loop)
cat("Mapping scores to samples (fast)...\n")
count_long <- melt(as.data.table(count_neo), id.vars = "junc.id",
                   variable.name = "sample", value.name = "count")
count_long[, sample := sub("\\.$", "", as.character(sample))]
count_long <- count_long[count > 0]

score_long <- merge(count_long,
                    top_neo[, .(junc.id, neo_id, score_average)],
                    by = "junc.id", allow.cartesian = TRUE)
score_wide <- dcast(score_long, neo_id ~ sample,
                    value.var = "score_average",
                    fun.aggregate = max, fill = NA_real_)

neo_matrix <- merge(top_neo[, .(neo_id, peptide, junc.id, hla_allele, score_average,
                                  mhcflurry_affinity_percentile, mhcflurry_binding_affinity,
                                  mhcflurry_presentation_score, type, fs, shared)],
                    score_wide, by = "neo_id", all.x = TRUE)
setDT(neo_matrix)
missing_samples <- setdiff(sample_names, colnames(neo_matrix))
if (length(missing_samples) > 0) neo_matrix[, (missing_samples) := NA_real_]
setcolorder(neo_matrix, c("neo_id","peptide","junc.id","hla_allele","score_average",
                           "mhcflurry_affinity_percentile","mhcflurry_binding_affinity",
                           "mhcflurry_presentation_score","type","fs","shared",
                           intersect(sample_names, colnames(neo_matrix))))
cat("Score matrix built:", nrow(neo_matrix), "neoantigens x",
    length(intersect(sample_names, colnames(neo_matrix))), "samples\n")


# Add summary statistics
cat("Adding summary statistics...\n")
neo_matrix[, samples_with_neo := rowSums(!is.na(.SD)), .SDcols = sample_names]
neo_matrix[, max_sample_score := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = sample_names]


# Create sample-focused summary
sample_summary <- data.table(
  sample_id = sample_names,
  total_neoantigens = colSums(!is.na(neo_matrix[, .SD, .SDcols = sample_names])),
  max_score = sapply(sample_names, function(s) max(neo_matrix[[s]], na.rm = TRUE)),
  avg_score = sapply(sample_names, function(s) mean(neo_matrix[[s]], na.rm = TRUE))
)


# Join with patient list for purity
patient_list <- fread(patient_list_file)
setkey(sample_summary, sample_id)
setkey(patient_list, sample_id)
sample_summary <- sample_summary[patient_list, nomatch = 0]


# Create HLA summary (dynamic for all alleles)
hla_summary <- top_neo[, .(
  count_neoantigens = .N,
  avg_score = mean(score_average),
  max_score = max(score_average),
  avg_percentile = mean(mhcflurry_affinity_percentile),
  avg_binding_affinity = mean(mhcflurry_binding_affinity),
  avg_presentation_score = mean(mhcflurry_presentation_score),
  unique_junctions = uniqueN(junc.id),
  fs_count = sum(fs == "fs"),
  in_frame_count = sum(fs == "in-frame")
), by = hla_allele][order(-max_score)]


# Create dynamic color palette for HLA summary plot
if (nrow(hla_summary) <= 12) {
  hla_summary_colors <- brewer.pal(max(3, nrow(hla_summary)), "Set3")
} else {
  hla_summary_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(hla_summary))
}
names(hla_summary_colors) <- hla_summary$hla_allele


# Write outputs
cat("Writing output files...\n")


# 1. Main neoantigen matrix (wide format)
fwrite(neo_matrix, paste0("top_neoantigens_matrix_", current_date, ".tsv"), sep = "\t")


# 2. Sample summary
fwrite(sample_summary, paste0("sample_neoantigen_summary_", current_date, ".tsv"), sep = "\t")


# 3. HLA summary
fwrite(hla_summary, paste0("hla_alleles_summary_", current_date, ".tsv"), sep = "\t")


# 4. Long format for specific analyses
neo_long <- melt(neo_matrix,
                 id.vars = c("neo_id", "peptide", "junc.id", "hla_allele", "score_average", 
                            "mhcflurry_affinity_percentile", "mhcflurry_binding_affinity", "mhcflurry_presentation_score",
                            "type", "fs", "shared"),
                 variable.name = "sample_id",
                 value.name = "sample_score",
                 measure.vars = sample_names)
neo_long <- neo_long[!is.na(sample_score)]
fwrite(neo_long, paste0("top_neoantigens_long_", current_date, ".tsv"), sep = "\t")


# Create visualization for HLA summary
p_top_hla <- ggplot(hla_summary, aes(x = reorder(hla_allele, -max_score), y = max_score, fill = hla_allele)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = hla_summary_colors) +
  theme_bw() +
  labs(x = "HLA Allele", y = "Highest Score", title = "HLA Alleles by Highest Neoantigen Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


pdf(paste0("figure_top_hla_distribution_", current_date, ".pdf"), 
    width = 12, height = 6)
print(p_top_hla)
dev.off()

CairoPNG(paste0("figure_top_hla_distribution_", current_date, ".png"), width = 12, height = 6, units = "in", res = 300)
print(p_top_hla)
dev.off()

# Sample-level visualization
p_sample_neo <- ggplot(sample_summary, aes(x = reorder(sample_id, -total_neoantigens), y = total_neoantigens, fill = purity)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Sample", y = "Number of Top Neoantigens", title = "Top Neoantigens per Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))


# p_sample_neo
pdf(paste0("figure_sample_neoantigens_", current_date, ".pdf"), 
    width = 12, height = 6)
print(p_sample_neo)
dev.off()

# Create p_hla_summary plot (add this before the saving code)
p_hla_summary <- ggplot(hla_summary, aes(x = reorder(hla_allele, -count_neoantigens), 
                                         y = count_neoantigens, 
                                         fill = hla_allele)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count_neoantigens), vjust = -0.5, size = 3) +
  scale_fill_manual(values = hla_summary_colors) +
  theme_bw() +
  labs(x = "HLA Allele", y = "Number of Top Neoantigens", 
       title = "HLA Alleles by Number of Top Neoantigens") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Then save it
pdf(paste0("figure_hla_summary_", current_date, ".pdf"), 
    width = 12, height = 6)
print(p_hla_summary)
dev.off()

# Also save as PNG if you want
CairoPNG(paste0("figure_hla_summary_", current_date, ".png"), 
         width = 12, height = 6, units = "in", res = 300)
print(p_hla_summary)
dev.off()


# Additional: Create a summary table of all HLA alleles found
hla_overview <- data.table(
  HLA_Allele = hla_alleles_all,
  Total_Neoantigens = sapply(hla_alleles_all, function(x) sum(df_all_map$hla_allele == x)),
  Top_10pct_Neoantigens = sapply(hla_alleles_all, function(x) sum(df_na_nj_map$hla_allele == x)),
  Mean_Score = sapply(hla_alleles_all, function(x) mean(df_all_map[hla_allele == x, score_average]))
)[order(-Total_Neoantigens)]

fwrite(hla_overview, paste0("hla_alleles_overview_", current_date, ".tsv"), sep = "\t")

cat("HLA alleles overview:\n")
print(hla_overview)


cat("Step 5 completed successfully!\n")
cat("Generated files:\n")
cat("1. top_neoantigens_matrix_", current_date, ".tsv - Main matrix format\n")
cat("2. sample_neoantigen_summary_", current_date, ".tsv - Sample-level summary\n")
cat("3. hla_alleles_summary_", current_date, ".tsv - HLA-level summary\n")
cat("4. top_neoantigens_long_", current_date, ".tsv - Long format for analysis\n")
cat("5. hla_alleles_overview_", current_date, ".tsv - HLA alleles overview\n")


print("All tasks completed successfully.")
