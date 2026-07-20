#!/usr/bin/env Rscript
# Step 15b: Plot Cross Analysis (NetMHCpan vs MHCFlurry) -- Concordance Tiers
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Visualise the concordance tier results from Step 15a.
#   Plot 1: Scatter -- NetMHCpan EL score vs MHCflurry presentation score
#            coloured by concordance tier (Tier1/2/3)
#   Plot 2: Scatter -- NetMHCpan EL rank vs MHCflurry affinity (nM)
#            with threshold lines at rank=0.5%, aff=500nM
#   Plot 3: Bar chart -- peptide counts per concordance tier
#
# Input (Step 15a):
#   alt_concordance_all_YYYYMMDD.tsv   -- clean ALT all tiers
#   wt_concordance_all_YYYYMMDD.tsv    -- WT native all tiers
#   cross_alg_all_nmers_YYYYMMDD.tsv   -- full join (for downstream)

###########################################################################
#  Step 0: Packages and config
###########################################################################

rm(list = ls(all.names = TRUE))
library(data.table)
library(ggplot2)

output_dir        <- Sys.getenv("OUTPUT_DIR")
directory_15      <- Sys.getenv("STEP15_OUTPUT_DIR")
directory_figures <- Sys.getenv("STEP15_FIGURES_DIR")

if (nchar(output_dir)   == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_15) == 0) directory_15 <- output_dir
if (nchar(directory_figures) == 0)
  directory_figures <- file.path(output_dir, "figures", "step15")

dir.create(directory_figures, showWarnings = FALSE, recursive = TRUE)

# Thresholds (must match 15a)
NMP_SB_RANK <- 0.5
NMP_WB_RANK <- 2.0
MHC_AFF_NM  <- 500

# Concordance tier colours
TIER_COLOURS <- c(
  "Tier1_HighConfidence"          = "#2dc653",  # green
  "Tier2_MediumConfidence"        = "#f4a100",  # orange
  "Tier3_Discordant_NMPstrong"    = "#e63946",  # red
  "Tier3_Discordant_MHCstrong"    = "#e07a5f"   # salmon
)

TIER_LABELS <- c(
  "Tier1_HighConfidence"          = "Tier 1: Both strong (SB + <500nM)",
  "Tier2_MediumConfidence"        = "Tier 2: WB + MHCflurry <500nM",
  "Tier3_Discordant_NMPstrong"    = "Tier 3: NMP strong, MHC weak",
  "Tier3_Discordant_MHCstrong"    = "Tier 3: MHC strong, NMP weak"
)

###########################################################################
#  Step 1: Detect date from 15a output files
###########################################################################

alt_scan <- list.files(directory_15,
                        pattern = "^alt_concordance_all_[0-9]{8}\\.tsv$")
if (length(alt_scan) == 0)
  stop("[ERROR] No alt_concordance_all_YYYYMMDD.tsv in: ", directory_15,
       "\n  Has Step 15a finished?")
current_date <- sub("^alt_concordance_all_([0-9]{8})\\.tsv$", "\\1",
                    alt_scan[which.max(
                      file.info(file.path(directory_15, alt_scan))$mtime)])

cat("[INFO] Detected 15a output date:", current_date, "\n")
cat("[INFO] Reading files from:      ", directory_15, "\n")
cat("[INFO] Writing figures to:      ", directory_figures, "\n")

###########################################################################
#  Step 2: Load 15a outputs
###########################################################################

alt_all <- fread(file.path(directory_15,
                            paste0("alt_concordance_all_", current_date, ".tsv")),
                 na.strings = c("", "NA"), quote = "")
wt_all  <- fread(file.path(directory_15,
                            paste0("wt_concordance_all_", current_date, ".tsv")),
                 na.strings = c("", "NA"), quote = "")

cat(sprintf("[INFO] ALT all tiers rows: %d\n", nrow(alt_all)))
cat(sprintf("[INFO] WT  all tiers rows: %d\n", nrow(wt_all)))
cat("\nALT concordance tier breakdown:\n")
print(alt_all[, .N, by = concordance_tier][order(-N)])
cat("\nWT concordance tier breakdown:\n")
print(wt_all[, .N, by = concordance_tier][order(-N)])

###########################################################################
#  Helper: generate the three plots for a given dataset (ALT or WT)
###########################################################################

make_plots <- function(dt, label, prefix) {

  dt[, netmhcpan_EL_score          := as.numeric(netmhcpan_EL_score)]
  dt[, netmhcpan_EL_rank           := as.numeric(netmhcpan_EL_rank)]
  dt[, mhcflurry_presentation_score := as.numeric(mhcflurry_presentation_score)]
  dt[, mhcflurry_affinity          := as.numeric(mhcflurry_affinity)]

  # Only plot rows where both scores are available
  dt_plot <- dt[!is.na(netmhcpan_EL_score) & !is.na(mhcflurry_presentation_score)]

  # Tier counts for titles
  tier_n <- dt[, .N, by = concordance_tier]

  # ------------------------------------------------------------------
  # Plot 1: EL score vs MHCflurry presentation score
  # ------------------------------------------------------------------
  n1 <- nrow(dt_plot[concordance_tier == "Tier1_HighConfidence"])
  n2 <- nrow(dt_plot[concordance_tier == "Tier2_MediumConfidence"])
  n3 <- nrow(dt_plot[grepl("Tier3", concordance_tier)])

  p1 <- ggplot(dt_plot,
               aes(x = netmhcpan_EL_score,
                   y = mhcflurry_presentation_score,
                   colour = concordance_tier)) +
    geom_point(size = 1.8, alpha = 0.6) +
    scale_colour_manual(values = TIER_COLOURS, labels = TIER_LABELS,
                        na.value = "grey80") +
    theme_minimal() +
    theme(plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
          text            = element_text(size = 12),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(ncol = 1, title = "Concordance tier",
                                 override.aes = list(size = 3, alpha = 1))) +
    labs(x     = "NetMHCpan 4.2 EL Score  [0-1, higher = stronger]",
         y     = "MHCflurry 2.0 Presentation Score  [0-1, higher = stronger]",
         title = paste0(label, ": EL Score vs Presentation Score\n",
                        "Tier1=", n1, "  Tier2=", n2, "  Tier3=", n3))

  f1 <- file.path(directory_figures,
                  paste0(prefix, "_scatter_ELscore_vs_pres_", current_date, ".pdf"))
  ggsave(f1, plot = p1, width = 8, height = 7, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", f1))

  # ------------------------------------------------------------------
  # Plot 2: EL rank vs MHCflurry affinity (nM) with threshold lines
  # Log scale on both axes -- lower rank + lower affinity = better
  # ------------------------------------------------------------------
  dt_plot2 <- dt[!is.na(netmhcpan_EL_rank) & !is.na(mhcflurry_affinity) &
                   mhcflurry_affinity > 0 & netmhcpan_EL_rank > 0]

  p2 <- ggplot(dt_plot2,
               aes(x = netmhcpan_EL_rank,
                   y = mhcflurry_affinity,
                   colour = concordance_tier)) +
    geom_point(size = 1.8, alpha = 0.6) +
    geom_vline(xintercept = NMP_SB_RANK, linetype = "dashed",
               colour = "#2dc653", linewidth = 0.7) +
    geom_vline(xintercept = NMP_WB_RANK, linetype = "dotted",
               colour = "#f4a100", linewidth = 0.7) +
    geom_hline(yintercept = MHC_AFF_NM,  linetype = "dashed",
               colour = "#2dc653", linewidth = 0.7) +
    annotate("text", x = NMP_SB_RANK, y = max(dt_plot2$mhcflurry_affinity, na.rm=TRUE),
             label = "SB 0.5%", hjust = -0.1, size = 3.5, colour = "#2dc653") +
    annotate("text", x = NMP_WB_RANK, y = max(dt_plot2$mhcflurry_affinity, na.rm=TRUE),
             label = "WB 2%",   hjust = -0.1, size = 3.5, colour = "#f4a100") +
    annotate("text", x = max(dt_plot2$netmhcpan_EL_rank, na.rm=TRUE),
             y = MHC_AFF_NM, label = "500nM", vjust = -0.5, hjust = 1,
             size = 3.5, colour = "#2dc653") +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_manual(values = TIER_COLOURS, labels = TIER_LABELS,
                        na.value = "grey80") +
    theme_minimal() +
    theme(plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
          text             = element_text(size = 12),
          legend.position  = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(ncol = 1, title = "Concordance tier",
                                 override.aes = list(size = 3, alpha = 1))) +
    labs(x     = "NetMHCpan 4.2 EL Rank % [log scale, lower = stronger]",
         y     = "MHCflurry 2.0 Affinity nM [log scale, lower = stronger]",
         title = paste0(label, ": EL Rank vs Affinity (nM)\n",
                        "Dashed lines = Tier 1 thresholds"))

  f2 <- file.path(directory_figures,
                  paste0(prefix, "_scatter_ELrank_vs_affinity_", current_date, ".pdf"))
  ggsave(f2, plot = p2, width = 8, height = 7, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", f2))

  # ------------------------------------------------------------------
  # Plot 3: Bar chart -- peptide counts per tier
  # ------------------------------------------------------------------
  tier_counts <- dt[, .(n_rows     = .N,
                         n_peptides = uniqueN(peptide)),
                     by = concordance_tier]
  tier_counts <- tier_counts[concordance_tier != "Excluded"]
  tier_counts[, tier_label := TIER_LABELS[concordance_tier]]
  tier_counts[is.na(tier_label), tier_label := concordance_tier]

  p3 <- ggplot(tier_counts,
               aes(x = reorder(tier_label, -n_peptides),
                   y = n_peptides,
                   fill = concordance_tier)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = scales::comma(n_peptides)),
              vjust = -0.4, size = 4) +
    scale_fill_manual(values = TIER_COLOURS, guide = "none") +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.15))) +
    theme_minimal() +
    theme(plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
          text             = element_text(size = 11),
          axis.text.x      = element_text(angle = 20, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank()) +
    labs(x     = "",
         y     = "Unique peptides",
         title = paste0(label, ": Peptides per concordance tier"))

  f3 <- file.path(directory_figures,
                  paste0(prefix, "_barplot_tier_counts_", current_date, ".pdf"))
  ggsave(f3, plot = p3, width = 8, height = 5, limitsize = FALSE)
  cat(sprintf("[INFO] Saved: %s\n", f3))
}

###########################################################################
#  Step 3: Generate plots for ALT (tumour-specific) and WT (native)
###########################################################################

cat("\n[STEP 3] Generating ALT plots...\n")
make_plots(alt_all, "ALT tumour-specific", "alt")

cat("\n[STEP 4] Generating WT native plots...\n")
make_plots(wt_all, "WT native immunopeptidome", "wt")

###########################################################################
#  Step 4: Write cross-analysis summary
###########################################################################

setwd(directory_15)
out_file <- paste0("cross_analysis_summary_nmp_mf_", current_date, ".tsv")
fwrite(alt_all, out_file, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("\n[INFO] Cross-analysis summary written: %s\n", out_file))

###########################################################################
#  Summary
###########################################################################

cat("\n=== 15b Summary ===\n")
cat("\nALT (tumour-specific):\n")
print(alt_all[concordance_tier != "Excluded",
              .(n_rows = .N, n_peptides = uniqueN(peptide)),
              by = concordance_tier][order(concordance_tier)])
cat("\nWT (native reference):\n")
print(wt_all[concordance_tier != "Excluded",
             .(n_rows = .N, n_peptides = uniqueN(peptide)),
             by = concordance_tier][order(concordance_tier)])
cat(sprintf("\nFigures written to: %s\n", directory_figures))
cat("\n[DONE] Step 15b complete.\n")