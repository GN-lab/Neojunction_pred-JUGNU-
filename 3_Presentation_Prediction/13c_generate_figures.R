#!/usr/bin/env Rscript
# Step 13c: Generate figures for NetMHCpan predictions
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Visualise the NetMHCpan 4.2 strong-binder candidates selected in
#          Step 13b. Produces histograms (EL rank distribution by allele and
#          by n-mer length) and pie charts (allele / length breakdown) for
#          the top 10 percentile of binders.
#
# Input:   netmhcpan_XXmer_selected_alleles_YYYYMMDD.tsv  (from Step 13b)
# Output:  PDF figures in STEP13_FIGURES_DIR
#
# X-axis:  netmhcpan_EL_rank  (lower = better presentation; inverted so
#          strong binders appear on the RIGHT of histograms)
#
# Note:    All unique peptides only (distinct on allele + peptide).

###########################################################################
#  Step 0: Load packages and config ---------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(data.table)
library(ggplot2)

# Directories from environment
# config.sh sets all STEPXX_OUTPUT_DIR to OUTPUT_DIR (results/).
# STEP13_FIGURES_DIR is not exported by config.sh -- fall back gracefully.
output_dir        <- Sys.getenv("OUTPUT_DIR")
directory_13      <- Sys.getenv("STEP13_OUTPUT_DIR")
directory_figures <- Sys.getenv("STEP13_FIGURES_DIR")

if (nchar(output_dir) == 0)   stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_13) == 0) directory_13 <- output_dir
if (nchar(directory_figures) == 0)
  directory_figures <- file.path(output_dir, "figures", "step13")

cat("[INFO] Reading 13b files from:", directory_13,     "\n")
cat("[INFO] Writing figures to:    ", directory_figures, "\n")

dir.create(directory_figures, showWarnings = FALSE, recursive = TRUE)

# Detect out_date from the most recently written 13b selected-allele file.
# Avoids mismatch when 13c runs on a different day from 13b.
sb_files <- list.files(directory_13, pattern = "^netmhcpan_09mer_selected_alleles_[0-9]{8}\\.tsv$")
if (length(sb_files) == 0) {
  stop("[ERROR] No netmhcpan_09mer_selected_alleles_YYYYMMDD.tsv found in: ", directory_13,
       "\n  Has Step 13b finished running?")
}
sb_mtimes <- file.info(file.path(directory_13, sb_files))$mtime
out_date  <- sub("^netmhcpan_09mer_selected_alleles_([0-9]{8})\\.tsv$", "\\1",
                 sb_files[which.max(sb_mtimes)])
cat("[INFO] Detected 13b out_date:", out_date, "\n")

###########################################################################
#  Step 1: Load selected-allele files from 13b ----------------------------
###########################################################################

load_nmer <- function(len) {
  f <- file.path(directory_13, paste0("netmhcpan_", len, "mer_selected_alleles_", out_date, ".tsv"))
  if (!file.exists(f)) {
    warning("[WARN] File not found: ", f)
    return(data.table())
  }
  dt <- fread(f, na.strings = c("", "NA"))
  # Keep unique peptide x allele combinations
  dt <- unique(dt, by = c("allele", "peptide"))
  dt[, len_label := paste0(len, "-mer")]
  cat(sprintf("[INFO] %smer: %d unique peptides loaded\n", len, nrow(dt)))
  dt
}

nmer_08 <- load_nmer("08")
nmer_09 <- load_nmer("09")
nmer_10 <- load_nmer("10")
nmer_11 <- load_nmer("11")

nmer_all <- rbindlist(list(nmer_08, nmer_09, nmer_10, nmer_11), fill = TRUE)

cat(sprintf("[INFO] Total unique peptides across all lengths: %d\n", nrow(nmer_all)))

###########################################################################
#  Colour palettes --------------------------------------------------------
###########################################################################

# Allele colours -- up to 10 distinct HLA alleles
# Add more hex codes if your cohort has >10 unique alleles
ALLELE_COLOURS <- c(
  "#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

# Length colours (8-11mer)
LEN_COLOURS <- c(
  "08-mer" = "#e6c229",
  "09-mer" = "#f17105",
  "10-mer" = "#d11149",
  "11-mer" = "#6610f2"
)

###########################################################################
#  Helper: safe colour mapping (handles variable number of alleles) -------
###########################################################################

allele_colour_map <- function(alleles) {
  u <- sort(unique(alleles))
  cols <- ALLELE_COLOURS[seq_along(u)]
  setNames(cols, u)
}

###########################################################################
#  Step 2: Full-distribution histograms (all SB peptides) ----------------
###########################################################################
# EL rank is inverted on x-axis: lower rank = stronger binder.
# We flip so strongest binders appear on the right (more intuitive).

setwd(directory_figures)

plot_hist <- function(dt, title, filename, fill_var, colour_map, legend_title) {
  if (nrow(dt) == 0) { warning("[WARN] No data for: ", title); return(invisible(NULL)) }

  p <- ggplot(dt, aes(x = -log10(netmhcpan_EL_rank + 0.001), fill = .data[[fill_var]])) +
    geom_histogram(bins = 50, colour = "white", linewidth = 0.1) +
    scale_fill_manual(values = colour_map) +
    theme_minimal() +
    theme(
      plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
      text             = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    xlab("-log10(NetMHCpan EL Rank)  [higher = stronger binder]") +
    ylab("Count") +
    ggtitle(title) +
    labs(fill = legend_title)

  ggsave(filename, plot = p, width = 8, height = 5, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", filename))
}

# Per-length histograms coloured by allele
for (len in c("08", "09", "10", "11")) {
  dt_i <- get(paste0("nmer_", len))
  if (nrow(dt_i) == 0) next
  cmap  <- allele_colour_map(dt_i$allele)
  n     <- nrow(dt_i)
  title <- paste0(len, "-mers | SB candidates (n=", n, ")")
  fname <- paste0("histogram_all_", len, "mer_n", n, "_", out_date, ".pdf")
  plot_hist(dt_i, title, fname, "allele", cmap, "HLA allele")
}

# All n-mers combined, coloured by length
if (nrow(nmer_all) > 0) {
  n     <- nrow(nmer_all)
  title <- paste0("All n-mers | SB candidates (n=", n, ")")
  fname <- paste0("histogram_all_nmers_n", n, "_", out_date, ".pdf")
  plot_hist(nmer_all, title, fname, "len_label", LEN_COLOURS, "n-mer length")
}

###########################################################################
#  Step 3: Top 10 percentile (strongest binders) -------------------------
###########################################################################

top10pct <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  # Sort ascending by EL rank (lower = better)
  setorder(dt, netmhcpan_EL_rank)
  n_top <- max(1L, as.integer(nrow(dt) * 0.1))
  dt[1:n_top]
}

top_08  <- top10pct(nmer_08)
top_09  <- top10pct(nmer_09)
top_10  <- top10pct(nmer_10)
top_11  <- top10pct(nmer_11)
top_all <- top10pct(nmer_all)

cat(sprintf("[INFO] Top 10pct counts: 08=%d 09=%d 10=%d 11=%d all=%d\n",
            nrow(top_08), nrow(top_09), nrow(top_10), nrow(top_11), nrow(top_all)))

# Top 10% histograms -- per length coloured by allele
for (len in c("08", "09", "10", "11")) {
  dt_i <- get(paste0("top_", len))
  if (nrow(dt_i) == 0) next
  cmap  <- allele_colour_map(dt_i$allele)
  n     <- nrow(dt_i)
  title <- paste0(len, "-mers | Top 10 Percentile (n=", n, ")")
  fname <- paste0("histogram_top10pct_", len, "mer_", out_date, ".pdf")
  plot_hist(dt_i, title, fname, "allele", cmap, "HLA allele")
}

# Top 10% all n-mers, coloured by length
if (nrow(top_all) > 0) {
  n     <- nrow(top_all)
  title <- paste0("All n-mers | Top 10 Percentile (n=", n, ")")
  fname <- paste0("histogram_top10pct_all_nmers_len_", out_date, ".pdf")
  plot_hist(top_all, title, fname, "len_label", LEN_COLOURS, "n-mer length")
}

# Top 10% all n-mers, coloured by allele
if (nrow(top_all) > 0) {
  cmap  <- allele_colour_map(top_all$allele)
  n     <- nrow(top_all)
  title <- paste0("All n-mers | Top 10 Percentile (n=", n, ")")
  fname <- paste0("histogram_top10pct_all_nmers_hla_", out_date, ".pdf")
  plot_hist(top_all, title, fname, "allele", cmap, "HLA allele")
}

###########################################################################
#  Step 4: Pie charts (top 10 percentile) ---------------------------------
###########################################################################

plot_pie <- function(dt, fill_var, title, filename, colour_map, legend_title) {
  if (nrow(dt) == 0) { warning("[WARN] No data for pie: ", title); return(invisible(NULL)) }

  counts <- as.data.frame(table(dt[[fill_var]]))
  colnames(counts) <- c("label", "count")

  p <- ggplot(counts, aes(x = "", y = count, fill = label)) +
    geom_bar(stat = "identity", width = 1, colour = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = colour_map) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      text       = element_text(size = 16)
    ) +
    labs(fill = legend_title)

  ggsave(filename, plot = p, width = 8, height = 5, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", filename))
}

# Per-length pie: allele distribution
for (len in c("08", "09", "10", "11")) {
  dt_i <- get(paste0("top_", len))
  if (nrow(dt_i) == 0) next
  cmap  <- allele_colour_map(dt_i$allele)
  title <- paste0(len, "-mers | Top 10 Percentile | Allele breakdown")
  fname <- paste0("pie_top10pct_", len, "mer_allele_", out_date, ".pdf")
  plot_pie(dt_i, "allele", title, fname, cmap, "HLA allele")
}

# All n-mers: allele distribution
if (nrow(top_all) > 0) {
  cmap  <- allele_colour_map(top_all$allele)
  title <- "All n-mers | Top 10 Percentile | Allele breakdown"
  fname <- paste0("pie_top10pct_all_nmers_allele_", out_date, ".pdf")
  plot_pie(top_all, "allele", title, fname, cmap, "HLA allele")
}

# All n-mers: length distribution
if (nrow(top_all) > 0) {
  title <- "All n-mers | Top 10 Percentile | Length breakdown"
  fname <- paste0("pie_top10pct_all_nmers_len_", out_date, ".pdf")
  plot_pie(top_all, "len_label", title, fname, LEN_COLOURS, "n-mer length")
}

###########################################################################
#  Summary ----------------------------------------------------------------
###########################################################################

cat("\n=== 13c Summary ===\n")
cat(sprintf("  Figures written to: %s\n", directory_figures))
n_pdf <- length(list.files(directory_figures, pattern = "\\.pdf$"))
cat(sprintf("  Total PDF files: %d\n", n_pdf))
cat("\n[DONE] Step 13c complete.\n")