#!/usr/bin/env Rscript
# Step 14d: Make histogram 
# Updated November 25, 2025 | Gaurav Raichand | The Institute of Cancer Research
# Purpose: Visualize MHCflurry 2.0 presentation scores (histograms + top 10% pies)

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggplot2)
library(data.table)

#  Load Directories -------------------------------------------------------
directory_14      <- Sys.getenv("STEP14_OUTPUT_DIR")  # where 14c outputs live
directory_figures <- Sys.getenv("OUTPUT_DIR")         # where figures will be written

# Helper: pick most recent file matching a pattern in current working directory
latest_file <- function(pattern) {
  files <- list.files(pattern = pattern)
  if (length(files) == 0) {
    stop("No files match pattern ", pattern, " in ", getwd())
  }
  files[which.max(file.info(files)$mtime)]
}

###########################################################################
#  Load MHCflurry Top Files (cohort-wide best alleles) --------------------
###########################################################################

setwd(directory_14)

# Use the *_top_YYYYMMDD.tsv files from Step 14c
f08 <- latest_file("^mhcflurry_08mer_top_\\d{8}\\.tsv$")
f09 <- latest_file("^mhcflurry_09mer_top_\\d{8}\\.tsv$")
f10 <- latest_file("^mhcflurry_10mer_top_\\d{8}\\.tsv$")
f11 <- latest_file("^mhcflurry_11mer_top_\\d{8}\\.tsv$")

message("Using input files:\n",
        "  ", f08, "\n",
        "  ", f09, "\n",
        "  ", f10, "\n",
        "  ", f11, "\n")

nmer_08 <- fread(f08, na = c("", "NA"))
nmer_09 <- fread(f09, na = c("", "NA"))
nmer_10 <- fread(f10, na = c("", "NA"))
nmer_11 <- fread(f11, na = c("", "NA"))

# Use today's date for figure filenames (decoupled from input dates)
current_date <- format(Sys.Date(), "%Y%m%d")

###########################################################################
#  Step 1: Remove duplicate n-mers (ADD THIS VERSION) --------------------
###########################################################################

cat("nmer_08:", nrow(nmer_08), "rows\n")
nmer_08_unique <- distinct(nmer_08[1:min(500000, nrow(nmer_08))], allele, peptide, n_flank, c_flank, .keep_all = TRUE)
cat("  -> unique:", nrow(nmer_08_unique), "\n")

cat("nmer_09:", nrow(nmer_09), "rows\n")
nmer_09_unique <- distinct(nmer_09[1:min(500000, nrow(nmer_09))], allele, peptide, n_flank, c_flank, .keep_all = TRUE)
cat("  -> unique:", nrow(nmer_09_unique), "\n")

cat("nmer_10:", nrow(nmer_10), "rows\n")
nmer_10_unique <- distinct(nmer_10[1:min(500000, nrow(nmer_10))], allele, peptide, n_flank, c_flank, .keep_all = TRUE)
cat("  -> unique:", nrow(nmer_10_unique), "\n")

cat("nmer_11:", nrow(nmer_11), "rows\n")
nmer_11_unique <- distinct(nmer_11[1:min(500000, nrow(nmer_11))], allele, peptide, n_flank, c_flank, .keep_all = TRUE)
cat("  -> unique:", nrow(nmer_11_unique), "\n")

nmer_all_unique <- bind_rows(nmer_08_unique, nmer_09_unique, nmer_10_unique, nmer_11_unique)
cat("All unique:", nrow(nmer_all_unique), "\n")

###########################################################################
#  Step 2: Stacked histograms for all n-mers ------------------------------
###########################################################################

setwd(directory_figures)

for (i in 1:5) {
  if (i == 1) {plot_i <- nmer_08_unique; n_count <- nrow(plot_i); title_i <- paste0("8-mers (n=", n_count, ")");  filename_i <- paste0("histogram_mhcflurry_all_08mer_n",  n_count, "_", current_date, ".pdf")}
  if (i == 2) {plot_i <- nmer_09_unique; n_count <- nrow(plot_i); title_i <- paste0("9-mers (n=", n_count, ")");  filename_i <- paste0("histogram_mhcflurry_all_09mer_n",  n_count, "_", current_date, ".pdf")}
  if (i == 3) {plot_i <- nmer_10_unique; n_count <- nrow(plot_i); title_i <- paste0("10-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_10mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 4) {plot_i <- nmer_11_unique; n_count <- nrow(plot_i); title_i <- paste0("11-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_11mer_n", n_count, "_", current_date, ".pdf")}
  if (i == 5) {plot_i <- nmer_all_unique; n_count <- nrow(plot_i); title_i <- paste0("All n-mers (n=", n_count, ")"); filename_i <- paste0("histogram_mhcflurry_all_all_nmers_n", n_count, "_", current_date, ".pdf")}
  
  p <- ggplot(plot_i, aes(x = mhcflurry_presentation_score, fill = allele)) +
    geom_histogram(bins = 50) + 
    theme_minimal() + 
    # LET GGPLOT AUTO-ASSIGN COLORS (handles 100+ alleles)
    xlab("MHCflurry 2.0 presentation score") +
    ylab("Count") +
    ggtitle(title_i) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      text       = element_text(size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 14)
    ) +
    labs(fill = "HLA allele")
  
  ggsave(filename_i, p, limitsize = FALSE, width = 12, height = 6)  # wider for more alleles
}

###########################################################################
#  Step 3: Top 10 percentile ----------------------------------------------
###########################################################################

# 1. Order by presentation score and take top 10%
nmer_08_ordered <- nmer_08_unique[order(nmer_08_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_09_ordered <- nmer_09_unique[order(nmer_09_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_10_ordered <- nmer_10_unique[order(nmer_10_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_11_ordered <- nmer_11_unique[order(nmer_11_unique$mhcflurry_presentation_score, decreasing = TRUE), ]
nmer_all_ordered <- nmer_all_unique[order(nmer_all_unique$mhcflurry_presentation_score, decreasing = TRUE), ]

for (i in 1:5) {
  if (i == 1) nmer_i <- nmer_08_ordered
  if (i == 2) nmer_i <- nmer_09_ordered
  if (i == 3) nmer_i <- nmer_10_ordered
  if (i == 4) nmer_i <- nmer_11_ordered
  if (i == 5) nmer_i <- nmer_all_ordered
  
  percentile_10 <- as.integer(nrow(nmer_i) / 10)
  
  if (i == 1) nmer_10percentile_08  <- nmer_i[1:percentile_10, ]
  if (i == 2) nmer_10percentile_09  <- nmer_i[1:percentile_10, ]
  if (i == 3) nmer_10percentile_10  <- nmer_i[1:percentile_10, ]
  if (i == 4) nmer_10percentile_11  <- nmer_i[1:percentile_10, ]
  if (i == 5) nmer_10percentile_all <- nmer_i[1:percentile_10, ]
}

# 2. Pie charts: allele distribution in top 10% (by n-mer length)
pie_08 <- as.data.frame(t(table(nmer_10percentile_08$allele)))[, c(2, 3)]
pie_09 <- as.data.frame(t(table(nmer_10percentile_09$allele)))[, c(2, 3)]
pie_10 <- as.data.frame(t(table(nmer_10percentile_10$allele)))[, c(2, 3)]
pie_11 <- as.data.frame(t(table(nmer_10percentile_11$allele)))[, c(2, 3)]

# All n-mers: derive length + allele tables
nmer_10percentile_all_edit <- nmer_10percentile_all %>%
  rowwise() %>%
  mutate(len = nchar(peptide)) %>%
  ungroup()

pie_all_nmer <- as.data.frame(t(table(nmer_10percentile_all_edit$len)))[, c(2, 3)]
pie_all_hla  <- as.data.frame(t(table(nmer_10percentile_all_edit$allele)))[, c(2, 3)]
cat("✅ Top 10% length calculation done\n")

# Pie charts
for (i in 1:6) {
  if (i == 1) {pie_i <- pie_08;       title_i <- "8-mers (Top 10%)";  filename_i <- paste0("pie_chart_mhcflurry_08mer_top10percentile_", current_date, ".pdf")}
  if (i == 2) {pie_i <- pie_09;       title_i <- "9-mers (Top 10%)";  filename_i <- paste0("pie_chart_mhcflurry_09mer_top10percentile_", current_date, ".pdf")}
  if (i == 3) {pie_i <- pie_10;       title_i <- "10-mers (Top 10%)"; filename_i <- paste0("pie_chart_mhcflurry_10mer_top10percentile_", current_date, ".pdf")}
  if (i == 4) {pie_i <- pie_11;       title_i <- "11-mers (Top 10%)"; filename_i <- paste0("pie_chart_mhcflurry_11mer_top10percentile_", current_date, ".pdf")}
  if (i == 5) {pie_i <- pie_all_hla;  title_i <- "All n-mers (Top 10%)"; filename_i <- paste0("pie_chart_mhcflurry_all_nmers_top10percentile_hla_distribution_", current_date, ".pdf")}
  if (i == 6) {pie_i <- pie_all_nmer; title_i <- "All n-mers (Top 10%)"; filename_i <- paste0("pie_chart_mhcflurry_all_nmers_top10percentile_len_distribution_", current_date, ".pdf")}
  
  if (i < 6) {
    p <- ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      ggtitle(title_i) +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text       = element_text(size = 20)
      ) +
      labs(fill = "HLA allele")
    
    ggsave(filename_i, p, limitsize = FALSE, width = 8, height = 5)
  }
  
  if (i == 6) {
    p <- ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      scale_fill_manual(values = c("#e6c229", "#f17105", "#d11149", "#6610f2")) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      ggtitle(title_i) +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text       = element_text(size = 20)
      ) +
      labs(fill = "n-mer length")
    
    ggsave(filename_i, p, limitsize = FALSE, width = 8, height = 5)
  }
}

###########################################################################
#  Step 4: Histogram of top 10% (SIMPLE VERSION) -------------------------
###########################################################################

cat("Making top 10% plots...\n")

# Limit top 10% to max 50k rows each for fast plotting
nmer_10percentile_08 <- head(nmer_08_ordered, min(50000, nrow(nmer_08_ordered)/10))
nmer_10percentile_09 <- head(nmer_09_ordered, min(50000, nrow(nmer_09_ordered)/10))
nmer_10percentile_10 <- head(nmer_10_ordered, min(50000, nrow(nmer_10_ordered)/10))
nmer_10percentile_11 <- head(nmer_11_ordered, min(50000, nrow(nmer_11_ordered)/10))

for (i in 1:4) {
  if (i == 1) {plot_i <- nmer_10percentile_08; title_i <- "8-mers (Top 10%)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_08mer_", current_date, ".pdf")}
  if (i == 2) {plot_i <- nmer_10percentile_09; title_i <- "9-mers (Top 10%)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_09mer_", current_date, ".pdf")}
  if (i == 3) {plot_i <- nmer_10percentile_10; title_i <- "10-mers (Top 10%)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_10mer_", current_date, ".pdf")}
  if (i == 4) {plot_i <- nmer_10percentile_11; title_i <- "11-mers (Top 10%)"; filename_i <- paste0("histogram_mhcflurry_top10percentile_11mer_", current_date, ".pdf")}
  
  p <- ggplot(plot_i, aes(x = mhcflurry_presentation_score, fill = allele)) +
    geom_histogram(bins = 30) + 
    theme_minimal() + 
    xlab("MHCflurry presentation score") + ylab("Count") + ggtitle(title_i) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(filename_i, p, width = 10, height = 6)
  cat("Saved:", filename_i, "\n")
}

# Combined plots
ggplot(nmer_10percentile_all_edit, aes(x = mhcflurry_presentation_score, fill = allele)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~len) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(paste0("top10percentile_overview_", current_date, ".pdf"), width = 12, height = 6)

cat("✅ All plots complete!\n")
