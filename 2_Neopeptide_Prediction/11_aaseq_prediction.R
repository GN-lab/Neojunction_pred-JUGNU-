#!/usr/bin/env Rscript
# Title: "Step 11: Amino-Acid Sequence Prediction"
# September 26, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: Predict the "template" or "isoform" of each candidate junction derived from Step 10
#          This is based on the following two information/assumptions:
#              1. Genomic coordinate
#              2. Premise that one or two ends of each junction will overlap with that of the canonical junction
#                 In the case of Exon Skipping, both edges will end at either two of the canonical junctions.
#                 In the case of A3 and A5, either of the edges is expected to be located at the canonical junction.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

start.time <- proc.time()  # Start timing -- moved after rm() so it actually survives (previously wiped immediately, which is why the later "Does start.time exist?" fallback check was needed at all)

# Accumulator for the exon -> CDS/AA coordinate map (NEW). One row per
# (junc.id, enst.model, exon_idx) for the ALT exon set, giving each exon's
# contribution in CDS-relative nt and AA-relative coordinates. This is what
# Step 12 joins against to attach peptide coordinates and flag n-mers that
# span an exon-exon boundary. Placed AFTER rm() above -- putting it before,
# like start.time originally was, meant it got wiped immediately and stayed
# undefined for the entire run.
EXON_COORD_MAP <- tibble::tibble()

# Suppress package startup messages and load explicitly
suppressPackageStartupMessages({
  library(BiocManager)
  library(tidyverse)
  library(XML)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)  # Updated for GRCh38
  library(GenomicRanges)
  library(AnnotationHub)
  library(Rsamtools)
  library(data.table)  # For faster loading
  library(BSgenome.Hsapiens.UCSC.hg38)  # Reference genome for getSeq
  library(Biostrings)  # For DNAString and translate
})

# Reference genomes - use both BSgenome and FASTA file
dna_bs <- BSgenome.Hsapiens.UCSC.hg38  # For main chromosomes

# Establish Directories ---------------------------------------------------

directory_00 <- Sys.getenv("EXTERNAL_PATH")
directory_03 <- Sys.getenv("STEP03_OUTPUT_DIR")
directory_04 <- Sys.getenv("STEP04_OUTPUT_DIR")
directory_10 <- Sys.getenv("STEP10_OUTPUT_DIR")
directory_11 <- Sys.getenv("OUTPUT_DIR")

# If environment variables are not set, use default paths
if (directory_00 == "") directory_00 <- "0_Input_Files/"
if (directory_03 == "") directory_03 <- "results/"
if (directory_04 == "") directory_04 <- "results/"
if (directory_10 == "") directory_10 <- "results/"
if (directory_11 == "") directory_11 <- "results/"

# Also load the FASTA file for alt contigs -- LAZY AND OPTIONAL (was a hard
# stop() before, even though this is only a fallback for sequences on
# non-primary-chromosome contigs). Protein-coding neojunctions almost never
# land on alt/patch scaffolds, so most runs never need this file at all. We
# only attempt to open it the first time get_sequence_smart() genuinely
# needs an alt-contig lookup, and if it's missing at that point we return NA
# (with a warning) instead of halting a multi-hour run over a file that may
# never actually be used.
fasta_file <- file.path(directory_00, "Homo_sapiens.GRCh38.dna.toplevel.fa")
dna_fa <- NULL
fasta_file_checked <- FALSE

get_alt_contig_fasta <- function() {
  if (is.null(dna_fa) && !fasta_file_checked) {
    fasta_file_checked <<- TRUE
    if (file.exists(fasta_file)) {
      dna_fa <<- FaFile(fasta_file)
      cat("[INFO] Alt-contig FASTA loaded:", fasta_file, "\n")
    } else {
      cat("[WARNING] Alt-contig FASTA not found at:", fasta_file,
          "-- any sequence on a non-primary-chromosome contig will return NA. ",
          "If this warning never appears again with an actual coordinate attached, ",
          "you never needed this file.\n")
    }
  }
  dna_fa
}

# Helper function to get sequence from appropriate source
get_sequence_smart <- function(seqname, start, end, strand) {
  if (length(seqname) == 0 || all(is.na(seqname))) {
    return(NA_character_)
  }
  seqname <- seqname[1]  # scalar per row
  start <- start[1]
  end <- end[1] 
  strand <- strand[1]

  main_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")
  
  if (seqname %in% main_chromosomes) {
    res <- tryCatch({
      gr <- GRanges(seqnames = seqname, ranges = IRanges(start, end), strand = strand)
      as.vector(getSeq(dna_bs, gr))
    }, error = function(e) NA_character_)
    if (!is.na(res[1])) return(res)
  }
  
  # Alt contigs → FASTA (strips 'chr' prefix)
  seqname_clean <- gsub("^chr", "", seqname)
  fa <- get_alt_contig_fasta()
  if (is.null(fa)) {
    cat("[WARNING] Needed alt-contig sequence at", seqname, ":", start, "-", end,
        "but no FASTA is available -- returning NA for this exon.\n")
    return(NA_character_)
  }
  tryCatch({
    gr <- GRanges(seqnames = seqname_clean, ranges = IRanges(start, end), strand = strand)
    open(fa)
    seq_val <- as.vector(getSeq(fa, gr))
    close(fa)
    seq_val
  }, error = function(e) {
    close(fa)
    NA_character_
  })
}

#  Load files: Auto-detect and Load Latest Input Files -----------------------

# Helper to pick the most recent file matching a pattern in the current working directory
latest_file <- function(pattern) {
  files <- list.files(pattern = pattern)
  if (length(files) == 0) stop("No files match pattern ", pattern, " in ", getwd())
  files[which.max(file.info(files)$mtime)]
}

# Get current date for output files (universal, not tied to input files)
current_date <- format(Sys.Date(), "%Y%m%d")

# From Step 10: PSR Tables for Neojunctions
setwd(directory_10)
psr_pattern <- "^PSR_Neojunctions_[0-9]{8}\\.tsv$"
psr_file    <- latest_file(psr_pattern)
dataframe_psr.neo <- data.table::fread(psr_file, na = c("", "NA"))

# Convert to tibble so dplyr verbs apply
dataframe_psr.neo <- tibble::as_tibble(dataframe_psr.neo)

# Extract needed columns - since your PSR file only has junction IDs, we need to match to genes
neo_ids <- dataframe_psr.neo %>% 
  dplyr::select(junc.id)

# From Step 03: GTF Table of Junctions Filtered for Protein-Coding & TPM
setwd(directory_03)
tpm_thresh  <- as.numeric(Sys.getenv("MIN_TPM", unset = "10"))
gtf_pattern <- paste0("^GTF_ProteinCoding_Filter", tpm_thresh, "_[0-9]{8}\\.tsv$")
gtf_file    <- latest_file(gtf_pattern)
dataframe_gtf <- data.table::fread(gtf_file, colClasses = list(character = "chr"))
dataframe_gtf <- tibble::as_tibble(dataframe_gtf)

# From Step 04: Annotated SJ Junction List Filtered by Step 03
setwd(directory_04)
sj_pattern  <- "^SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_[0-9]{8}\\.tsv$"
sj_file     <- latest_file(sj_pattern)
dataframe_sj <- data.table::fread(sj_file, colClasses = list(character = "chr"))
dataframe_sj <- tibble::as_tibble(dataframe_sj)

# From Step 00: Original Annotated Junctions List
setwd(directory_00)

# Load the raw file (5 columns: seqname, start, end, strand, motif)
sj_raw <- data.table::fread(
  "sjdbList.fromGTF.out.tab",
  header       = FALSE,
  sep          = "\t",
  quote        = "",
  showProgress = FALSE
)

# Assign column names based on actual structure
data.table::setnames(
  sj_raw,
  c("seqname", "int.start", "int.end", "strand", "motif")
)

# Convert to tibble, prefix seqname with "chr", and select needed columns
dataframe_sj.orig <- sj_raw %>%
  tibble::as_tibble() %>%
  dplyr::mutate(chr = paste0("chr", seqname)) %>%  # e.g., "1" becomes "chr1"
  dplyr::select(chr, int.start, int.end, strand)

# Confirmation Prints
cat("Loaded PSR file: ", basename(psr_file), "\n")
cat("Loaded GTF file: ", basename(gtf_file), "\n")
cat("Loaded SJ file:  ", basename(sj_file), "\n\n")

# Quick dimension checks
cat("Dimensions:\n")
cat("neo_ids:", dim(neo_ids), "\n")
cat("GTF:", dim(dataframe_gtf), "\n")
cat("SJ:", dim(dataframe_sj), "\n")
cat("SJ orig:", dim(dataframe_sj.orig), "\n")

###########################################################################
#  Step 1: Edit the Loaded Dataframes -------------------------------------
###########################################################################

# Since your PSR file only contains junction IDs without gene information,
# we need to recreate the gene matching logic from the original script

# Add informational columns for the neojunction IDs (obtained from Step 10)
neo_ids.2 <- neo_ids %>% 
  dplyr::mutate(chr = gsub("chr", "", sapply(strsplit(junc.id, ":"), "[[", 1))) %>% 
  dplyr::mutate(strand = sapply(strsplit(junc.id, ":"), "[[", 2)) %>% 
  dplyr::mutate(int.start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3))) + 1) %>% 
  dplyr::mutate(int.end = as.numeric(gsub(".*-", "", sapply(strsplit(junc.id, ":"), "[[", 3))))

# Optional smoke-test toggle (NEW): set TEST_N_JUNCTIONS before running to
# limit this run to just the first N junctions, e.g. `export
# TEST_N_JUNCTIONS=5`. Unset (the default) means a normal full run --
# this has zero effect otherwise. Everything downstream derives from
# neo_ids.2, so subsetting it here cascades correctly through the rest of
# the script without touching any loop logic.
test_n <- Sys.getenv("TEST_N_JUNCTIONS", unset = "")
if (test_n != "") {
  test_n <- as.integer(test_n)
  cat("[TEST MODE] TEST_N_JUNCTIONS is set -- limiting to the first", test_n,
      "of", nrow(neo_ids.2), "junctions.\n")
  neo_ids.2 <- neo_ids.2 %>% dplyr::slice(1:min(test_n, nrow(neo_ids.2)))
}

# Edit the column names of the complete junction list (obtained from Step 00)
dataframe_sj.orig <- dataframe_sj.orig %>% 
  dplyr::mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
  dplyr::select(junc.id, chr, strand, int.start, int.end)

# Match junctions to genes like in the original script
neo_ids.3 <- tibble()  # Initialize as empty tibble
for(i in 1:nrow(neo_ids.2)) {
  if(i %% 100 == 0) print(i / nrow(neo_ids.2) * 100)  # Progress Bar
  
  # 1. Slice out row(i) from the list of neoantigens (Step 10)
  neo_ids.i <- neo_ids.2 %>% 
    dplyr::slice(i)
  
  JUNC.ID <- neo_ids.i %>% dplyr::pull(junc.id)
  CHR <- neo_ids.i %>% dplyr::pull(chr)
  STRAND <- neo_ids.i %>% dplyr::pull(strand)
  INT.START <- neo_ids.i %>% dplyr::pull(int.start)
  INT.END <- neo_ids.i %>% dplyr::pull(int.end)
  
  # 2. Identify whether the neojunction has a canonical junction from the annotated SJ's
  sj.hit <- dataframe_sj %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>% 
    dplyr::filter(int.start == INT.START | int.end == INT.END) %>% 
    dplyr::pull(junc.id)
  
  # 3. In the case where the neojunction does not hit any canonical junction, 
  # see if they fall in any of the non-protein coding transcripts
  sj.orig.hit <- dataframe_sj.orig %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>% 
    dplyr::filter(int.start == INT.START | int.end == INT.END) %>% 
    dplyr::pull(junc.id)
  
  # 4. Based on #2/3, label the neojunctions
  sj.hit <- if(length(sj.hit) == 0 & length(sj.orig.hit) == 0){
    "no.hit"
  } else if(length(sj.hit) == 0 & length(sj.orig.hit) > 0){ 
    "out.of.scope"
  } else {                                                
    paste(sj.hit, collapse = ";")
  }
  
  # 5. Using the GTF file, identify the gene by which the neojunction's ends lie between
  gtf.i <- dataframe_gtf %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>% 
    dplyr::filter(start < INT.START & INT.END < end)
  
  # 6. Add columns for the Canonical ID, Gene Symbol, ENSG, and ENST 
  neo_ids.i <- neo_ids.i %>% 
    dplyr::mutate(canonical = sj.hit) %>%
    dplyr::mutate(symbol = ifelse(nrow(gtf.i) > 0, 
                                 paste(gtf.i %>% dplyr::pull(symbol) %>% unique(), collapse = ";"), 
                                 "")) %>%
    dplyr::mutate(ensg = ifelse(nrow(gtf.i) > 0, 
                               paste(gtf.i %>% dplyr::pull(ensg) %>% unique(), collapse = ";"), 
                               "")) %>%
    dplyr::mutate(enst = ifelse(nrow(gtf.i) > 0, 
                               paste(gtf.i %>% dplyr::pull(enst) %>% unique(), collapse = ";"), 
                               ""))
  
  # 7. Concatenate rows together for every iteration
  if(i == 1){
    neo_ids.3 <- neo_ids.i
  } else {
    neo_ids.3 <- neo_ids.3 %>% 
      bind_rows(neo_ids.i)
  }
}

# Now use the matched data for amino acid prediction
neo_ids_final <- neo_ids.3 %>% 
  dplyr::select(junc.id, symbol, ensg, enst) %>% 
  dplyr::filter(symbol != "" & !is.na(symbol))  # Filter out junctions without gene annotations

cat("Junctions with gene annotations:", nrow(neo_ids_final), "out of", nrow(neo_ids.3), "\n")

if(nrow(neo_ids_final) == 0) {
  stop("No junctions found with gene annotations. Cannot proceed with amino acid prediction.")
}

###########################################################################
#  Step 2: Amino Acid Sequence Prediction (ORIGINAL PIPELINE) ------------
###########################################################################

neo_ids.4 <- tibble()  # Initialize as empty tibble
for(i in 1:nrow(neo_ids_final)) {
  print(i / nrow(neo_ids_final) * 100)  # Progress Bar
  
  # 1. Slice out row(i) from the list of neoantigens (Step 10)
  neo_ids.i <- neo_ids_final %>% 
    dplyr::slice(i)
  
  JUNC.ID <- neo_ids.i %>% dplyr::pull(junc.id)
  ENSG <- neo_ids.i %>% dplyr::pull(ensg) %>% strsplit(";") %>% unlist()
  ENSTs <- neo_ids.i %>% dplyr::pull(enst) %>% strsplit(";") %>% unlist()
  CHR <- gsub("chr", "", sapply(strsplit(JUNC.ID, ":"), "[[", 1))
  STRAND <- sapply(strsplit(JUNC.ID, ":"), "[[", 2)
  INT.START <- as.numeric(gsub("-.*", "", sapply(strsplit(JUNC.ID, ":"), "[[", 3))) + 1
  INT.END <- as.numeric(gsub(".*-", "", sapply(strsplit(JUNC.ID, ":"), "[[", 3)))
  
  # load exonic info of all the protein-coding transcripts (isoforms) from the ENSG
  suppressWarnings({
    exons <- exons(
      EnsDb.Hsapiens.v86,  # Updated for GRCh38
      columns = c("tx_id", "exon_idx"), 
      filter = list(GeneIdFilter(ENSG), TxBiotypeFilter("protein_coding"))
    ) %>% 
      as_tibble() %>% 
      dplyr::mutate(seqnames = paste0("chr", as.character(seqnames)),  # Convert Rle to character
                    strand = as.character(strand)) %>%  # Convert Rle to character
      dplyr::arrange(tx_id, exon_idx) %>% 
      dplyr::select(-exon_id, -gene_id, -tx_biotype)
    
    # Get sequences using our smart function
    exons <- exons %>% 
      dplyr::filter(!is.na(seqnames), !is.na(start), !is.na(end), !is.na(strand)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(seq = get_sequence_smart(seqnames, start, end, strand)) %>%
      dplyr::ungroup()
    
    # Remove any exons where we couldn't get sequence
    exons <- exons %>% 
      dplyr::filter(!is.na(seq))
  })
  
  # Skip if no exons left
  if (nrow(exons) == 0) {
    cat("Skipping gene", ENSG, "- could not get sequences for any exons\n")
    next
  }
  
  # check if 1 or 2 of the edge(s) of the junction will match the canonical junc --------
  
  ENSTs.test <- exons %>% # tier 1 # retain ENSTs that fulfill both (highly-expressed / 1 or 2 edge(s) match)
    dplyr::filter(tx_id %in% ENSTs & (end == INT.START - 1 | INT.END + 1 == start)) %>% 
    dplyr::pull(tx_id) %>% 
    unique()
  
  if (length(ENSTs.test) == 0) { # tier 2 # retain ENSTs that fulfill either
    ENSTs.test <- exons %>% 
      dplyr::filter(tx_id %in% ENSTs | (end == INT.START - 1 | INT.END + 1 == start)) %>% 
      dplyr::pull(tx_id) %>% 
      unique()
  }
  
  if (length(ENSTs.test) == 0) { # tier 3 # others
    ENSTs.test <- exons %>% 
      dplyr::pull(tx_id) %>% 
      unique()
  }
  
  if (length(ENSTs.test) == 0) { # break in case of zero-hit
    print("check required")
    next
  }
  
  # genomic range-based filter: further validation: exclude the isoforms of which genomic range don't include JUNC.ID
  ENSTs.test <- transcripts(
    EnsDb.Hsapiens.v86, 
    columns = c("tx_id"), 
    filter = list(TxIdFilter(ENSTs.test))
  ) %>% 
    as_tibble() %>% 
    dplyr::filter(start < INT.START & INT.END < end) %>% 
    dplyr::pull(tx_id)
  
  # transcripts to be excluded (manually added) -----------------------------
  
  ENSTs.exclude <- c(
    "ENST00000400991", # "no protein" in ensembl 
    "ENST00000368242", # "no protein" in ensembl
    "ENST00000404436", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000429422", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000423049", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000316851", # the beginning part of aaseq is a bit different from that in Ensembl. 
    "ENST00000084795", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000467825"  # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
  )
  
  ENSTs.test <- ENSTs.test[!ENSTs.test %in% ENSTs.exclude]
  
  # loop for each tx --------------------------------------------------------
  
  if (length(ENSTs.test) == 0) {
    next
  } else {
    for (j in 1:length(ENSTs.test)) {
      ENST.j <- ENSTs.test[j]
      
      exons.j <- exons %>% 
        dplyr::filter(tx_id == ENST.j)
      
      # prep for AS type --------------------------------------------------------------------
      
      exons.j <- exons.j %>% 
        dplyr::arrange(start)
      
      if (nrow(exons.j) == 1) {
        exons.introns.j <- exons.j %>% 
          dplyr::mutate(lab = "exon") %>% 
          dplyr::mutate(idx = paste0("E", exon_idx)) %>% 
          dplyr::select(start, end, lab, id = exon_idx, idx)
      } else {
        introns.j <- exons.j %>% 
          dplyr::select(start, end, exon_idx) %>% 
          dplyr::slice(1:(nrow(exons.j) - 1)) %>% 
          dplyr::mutate(int.start = end + 1) %>% 
          dplyr::mutate(int.end = exons.j$start[2:nrow(exons.j)] - 1) %>% 
          dplyr::mutate(int.idx = if (STRAND == "+") 1:(nrow(exons.j) - 1) else (nrow(exons.j) - 1):1) %>% 
          dplyr::select(int.start, int.end, int.idx)
        
        exons.introns.j <- bind_rows(
          exons.j %>% 
            dplyr::mutate(lab = "exon") %>% 
            dplyr::mutate(idx = paste0("E", exon_idx)) %>% 
            dplyr::select(start, end, lab, id = exon_idx, idx),
          introns.j %>% 
            dplyr::mutate(lab = "intron") %>% 
            dplyr::mutate(idx = paste0("I", int.idx)) %>% 
            dplyr::select(start = int.start, end = int.end, lab, id = int.idx, idx)
        )
      }
      
      # AS type -----------------------------------------------------------------
      TYPE <- "UNDEFINED"
      NOTE <- "UNDEFINED"
      # discriminate 9 different patterns
      
      # 1) ES; exon-skipping
      if (nrow(exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(end + 1 == INT.START | INT.END == start - 1)) == 2) {
        if (nrow(exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(INT.START < start & end < INT.END)) > 1) {
          TYPE <- "ES"
          IDX <- exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(end == INT.START - 1 | INT.END + 1 == start) %>% dplyr::pull(idx)
          NOTE <- paste0("lt.", IDX[1], "-rt.", IDX[2])
        }
      } else if (nrow(exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(INT.END + 1 == start)) == 1) {  # LEFT
        RT.SIDE <- exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(INT.END + 1 == start)
        LT.SIDE <- exons.introns.j %>% dplyr::filter(start <= INT.START & INT.START <= end)
        
        if (LT.SIDE$lab == "intron" & (if (STRAND == "+") LT.SIDE$id == RT.SIDE$id - 1 else LT.SIDE$id == RT.SIDE$id)) {  # 2) LEFT - GAIN
          TYPE <- if (STRAND == "+") "A5.gain" else "A3.gain"
        } else {  # 3) LEFT - LOSS
          TYPE <- if (STRAND == "+") "A5.loss" else "A3.loss"
        }
        NOTE <- paste0("lt.", tolower(LT.SIDE$idx), "-rt.", RT.SIDE$idx)
      } else if (nrow(exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(end == INT.START - 1)) == 1) {  # RIGHT
        LT.SIDE <- exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter(end == INT.START - 1)
        RT.SIDE <- exons.introns.j %>% dplyr::filter(start <= INT.END & INT.END <= end)
        
        if (RT.SIDE$lab == "intron" & (if (STRAND == "+") LT.SIDE$id == RT.SIDE$id else LT.SIDE$id == RT.SIDE$id + 1)) {  # 4) RIGHT - GAIN
          TYPE <- if (STRAND == "+") "A3.gain" else "A5.gain"
        } else {  # 5) RIGHT - LOSS
          TYPE <- if (STRAND == "+") "A3.loss" else "A5.loss"
        }
        NOTE <- paste0("lt.", LT.SIDE$idx, "-rt.", tolower(RT.SIDE$idx))
      } else if (nrow(exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter((start < INT.START & INT.START < end) & (start < INT.END & INT.END < end))) == 1) {  # 6) WITHIN SINGLE EXONIC REGION
        TYPE <- "JUNC.WITHIN.EXON"
        NOTE <- paste0("within.", tolower((exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter((start < INT.START & INT.START < end) & (start < INT.END & INT.END < end)) %>% dplyr::pull(idx))))
      } else if (nrow(exons.introns.j %>% dplyr::filter(lab == "intron") %>% dplyr::filter((start < INT.START & INT.START < end) & (start < INT.END & INT.END < end))) == 1) {  # 7) WITHIN SINGLE INTRONIC REGION
        TYPE <- "JUNC.WITHIN.INTRON"
        NOTE <- paste0("within.", tolower((exons.introns.j %>% dplyr::filter(lab == "intron") %>% dplyr::filter((start < INT.START & INT.START < end) & (start < INT.END & INT.END < end)) %>% dplyr::pull(idx))))
      } else {  # 8) THE OTHERS
        LT.SIDE <- exons.introns.j %>% dplyr::filter(start <= INT.START & INT.START <= end)
        RT.SIDE <- exons.introns.j %>% dplyr::filter(start <= INT.END & INT.END <= end)
        TYPE <- "OTHERS"
        NOTE <- paste0("lt.", tolower(LT.SIDE$idx), "-rt.", tolower(RT.SIDE$idx))
      }
      
      # nucleotide sequence of "altered" tx -------------------------------------      
      
      # extract "involved" exons
      
      LT.SIDE <- exons.introns.j %>% 
        dplyr::arrange(start) %>% 
        dplyr::mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(start <= INT.START & INT.START <= end)
      
      RT.SIDE <- exons.introns.j %>% 
        dplyr::arrange(start) %>% 
        dplyr::mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(start <= INT.END & INT.END <= end)
      
      EXONS.NOT.INVOLVED <- exons.introns.j %>% 
        dplyr::arrange(start) %>% 
        dplyr::mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(nr < (LT.SIDE %>% dplyr::pull(nr)) | (RT.SIDE %>% dplyr::pull(nr)) < nr) %>% 
        dplyr::filter(lab == "exon")
      
      # FIX: Make sure EXONS.NOT.INVOLVED has the correct column name for joining
      if(nrow(EXONS.NOT.INVOLVED) > 0) {
        EXONS.NOT.INVOLVED <- EXONS.NOT.INVOLVED %>% 
          dplyr::rename(exon_idx = id)  # Rename 'id' to 'exon_idx' for proper joining
      }
      
      # edit "involved" exons
      
      for (k in 1:2) { # 1: LT; 2: RT
        NEW.EXON <- list(LT.SIDE, RT.SIDE)[[k]] %>% 
          dplyr::mutate(end = if (k == 1) INT.START - 1 else end) %>% 
          dplyr::mutate(start = if (k == 2) INT.END + 1 else start) %>% 
          dplyr::mutate(idx = id + 0.5)
        
        if (!is.null(NEW.EXON$start) && !is.null(NEW.EXON$end) && NEW.EXON$start > NEW.EXON$end) {
          NEW.EXON <- NULL
        } else {
          # Use our smart function to get sequence
          seq_val <- get_sequence_smart(paste0("chr", CHR), NEW.EXON$start, NEW.EXON$end, STRAND)
          
          NEW.EXON <- GRanges(
            seqnames = paste0("chr", CHR),
            strand = STRAND, 
            ranges = IRanges(NEW.EXON$start, NEW.EXON$end)
          ) %>% 
            as_tibble() %>% 
            dplyr::mutate(tx_id = ENST.j) %>% 
            dplyr::mutate(exon_idx = NEW.EXON$idx) %>% 
            dplyr::mutate(seqnames = as.character(seqnames)) %>% 
            dplyr::mutate(seq = seq_val)
        }
        if (k == 1) {
          LT.NEW.EXON <- NEW.EXON
        } else {
          RT.NEW.EXON <- NEW.EXON
        }
      }
      
      # FIX: Handle the case where EXONS.NOT.INVOLVED might be empty
      if(nrow(EXONS.NOT.INVOLVED) > 0) {
        exons.new <- exons.j %>% 
          dplyr::semi_join(EXONS.NOT.INVOLVED, by = "exon_idx")
      } else {
        exons.new <- exons.j[0,]  # Empty dataframe with same structure
      }
      
      # Add the new exons if they exist
      if(!is.null(LT.NEW.EXON)) {
        exons.new <- exons.new %>% bind_rows(LT.NEW.EXON)
      }
      if(!is.null(RT.NEW.EXON)) {
        exons.new <- exons.new %>% bind_rows(RT.NEW.EXON)
      }
      
      if (STRAND == "+") {
        exons.new <- exons.new %>% dplyr::arrange(start)
      } else {
        exons.new <- exons.new %>% dplyr::arrange(desc(start))
      }
      
      exons.j <- exons.j %>% dplyr::arrange(exon_idx)
      
      # tx.seq ------------------------------------------------------------------
      
      # transcript dna seq - wildtype and altered
      tx.seq.wt <- paste(exons.j$seq, collapse = "")
      tx.seq.alt <- if(nrow(exons.new) > 0) {
        paste(exons.new$seq, collapse = "")
      } else {
        ""  # Handle case where no exons remain
      }
      
      # rm. 5'UTR / 3'UTR -----------------------------------------------------------
      
      # remove "untranslated regions" from the above dna seq
      suppressWarnings({
        utr5 <- fiveUTRsByTranscript(EnsDb.Hsapiens.v86, columns = c("tx_id", "exon_idx"), filter = list(TxIdFilter(ENST.j)))
        ln.utr5 <- if (length(utr5) == 0) 0 else sum(utr5 %>% as_tibble() %>% dplyr::pull(width))
        
        utr3 <- threeUTRsByTranscript(EnsDb.Hsapiens.v86, columns = c("tx_id", "exon_idx"), filter = list(TxIdFilter(ENST.j)))
        ln.utr3 <- if (length(utr3) == 0) 0 else sum(utr3 %>% as_tibble() %>% dplyr::pull(width))
      })
      
      cds.seq.wt <- substr(tx.seq.wt, (ln.utr5 + 1), (nchar(tx.seq.wt) - ln.utr3))
      cds.seq.alt <- if(nchar(tx.seq.alt) > 0) {
        substr(tx.seq.alt, (ln.utr5 + 1), (nchar(tx.seq.alt) - ln.utr3))
      } else {
        ""
      }
      
      # exon -> CDS/AA coordinate map + validation (NEW) -------------------------
      # For each exon (in the SAME order used to build tx.seq via paste(collapse)),
      # compute how much of that exon's own sequence falls inside the CDS (i.e.
      # outside the 5'/3' UTR), in both CDS-relative nt coordinates and AA
      # coordinates. Then reconstruct the CDS purely from this per-exon map and
      # check it's byte-identical to cds.seq.wt/alt built above by simple
      # substr(). A mismatch means the coordinate math is wrong -- almost always
      # right at an exon boundary -- and is flagged rather than trusted silently.
      
      add_cds_coords <- function(exon_df, ln.utr5, ln.utr3, tx_total_len) {
        widths <- nchar(exon_df$seq)
        widths[is.na(widths)] <- 0
        tx_end   <- cumsum(widths)
        tx_start <- tx_end - widths + 1
        exon_df$tx_start <- tx_start
        exon_df$tx_end   <- tx_end
        
        cds_lo <- ln.utr5 + 1
        cds_hi <- tx_total_len - ln.utr3
        
        ov_start <- pmax(tx_start, cds_lo)
        ov_end   <- pmin(tx_end,   cds_hi)
        has_cds  <- ov_start <= ov_end
        
        exon_df$cds_start <- ifelse(has_cds, ov_start - ln.utr5, NA_real_)
        exon_df$cds_end   <- ifelse(has_cds, ov_end   - ln.utr5, NA_real_)
        
        local_start <- ifelse(has_cds, ov_start - tx_start + 1, NA_real_)
        local_end   <- ifelse(has_cds, ov_end   - tx_start + 1, NA_real_)
        exon_df$cds_seq <- ifelse(has_cds, substr(exon_df$seq, local_start, local_end), NA_character_)
        
        # AA coordinates: 1-based, inclusive. A codon split across two exons
        # gets ceiling()'d on both ends, so it's correctly attributed to BOTH
        # exons -- this is deliberate, not a bug: that codon really did come
        # from both.
        exon_df$aa_start <- ifelse(has_cds, ceiling(exon_df$cds_start / 3), NA_real_)
        exon_df$aa_end   <- ifelse(has_cds, ceiling(exon_df$cds_end   / 3), NA_real_)
        
        exon_df
      }
      
      COORD_CHECK <- tryCatch({
        exons.j.map <- add_cds_coords(exons.j, ln.utr5, ln.utr3, nchar(tx.seq.wt))
        exons.new.map <- if (nrow(exons.new) > 0) {
          add_cds_coords(exons.new, ln.utr5, ln.utr3, nchar(tx.seq.alt))
        } else {
          exons.new[0, ]
        }
        
        reconstructed_cds_wt <- paste(exons.j.map$cds_seq[!is.na(exons.j.map$cds_seq)], collapse = "")
        reconstructed_cds_alt <- if (nrow(exons.new.map) > 0) {
          paste(exons.new.map$cds_seq[!is.na(exons.new.map$cds_seq)], collapse = "")
        } else {
          ""
        }
        
        check_wt  <- identical(reconstructed_cds_wt,  cds.seq.wt)
        check_alt <- identical(reconstructed_cds_alt, cds.seq.alt)
        
        if (!check_wt || !check_alt) {
          cat("[COORD MAP MISMATCH] junc.id =", JUNC.ID, "enst.model =", ENST.j,
              "| wt_ok =", check_wt, "| alt_ok =", check_alt, "\n")
        }
        
        # Persist the ALT exon->coordinate map -- only rows that actually
        # contribute to the CDS, only if the reconstruction check passed for
        # ALT (a failed map is worse than no map -- don't propagate bad
        # coordinates downstream to Step 12).
        if (check_alt && nrow(exons.new.map) > 0) {
          exon_rows <- exons.new.map %>%
            dplyr::filter(!is.na(cds_start)) %>%
            dplyr::mutate(junc.id = JUNC.ID, enst.model = ENST.j) %>%
            dplyr::select(junc.id, enst.model, exon_idx, seqnames, start, end, strand,
                          cds_start, cds_end, aa_start, aa_end)
          EXON_COORD_MAP <<- dplyr::bind_rows(EXON_COORD_MAP, exon_rows)
        }
        
        if (check_wt && check_alt) "pass" else "MISMATCH"
      }, error = function(e) {
        cat("[COORD MAP ERROR] junc.id =", JUNC.ID, "enst.model =", ENST.j,
            "|", conditionMessage(e), "\n")
        "ERROR"
      })
      
      # translate ------------------------------------------------------------------
      
      aa.seq.wt <- if(nchar(cds.seq.wt) > 0) {
        suppressWarnings(as.character(translate(DNAString(cds.seq.wt))))
      } else {
        ""
      }
      
      aa.seq.alt <- if(nchar(cds.seq.alt) > 0) {
        suppressWarnings(as.character(translate(DNAString(cds.seq.alt))))
      } else {
        ""
      }
      
      # remove "*(stop-codon)" in the end ---------------------------------------
      
      aa.seq.wt <- if (nchar(aa.seq.wt) > 0 && str_sub(aa.seq.wt, -1, -1) == "*") str_sub(aa.seq.wt, 1, -2) else aa.seq.wt
      aa.seq.alt <- if (nchar(aa.seq.alt) > 0 && str_sub(aa.seq.alt, -1, -1) == "*") str_sub(aa.seq.alt, 1, -2) else aa.seq.alt
      
      # i = 224, SELENOW gene ---------------------------------------------------------------
      
      # "Selenocystein (Sec or U)" ref: https://en.wikipedia.org/wiki/Selenocysteine"
      aa.seq.wt <- gsub("MALAVRVVYCGA\\*GYKSKYLQLKKKLE", "MALAVRVVYCGAUGYKSKYLQLKKKLE", aa.seq.wt)
      aa.seq.alt <- gsub("MALAVRVVYCGA\\*GYKSKYLQLKKKLE", "MALAVRVVYCGAUGYKSKYLQLKKKLE", aa.seq.alt)
      
      # evaluate some elements --------------------------------------------------------------------
      
      # frame-shift or not 
      FS <- if (abs(nchar(cds.seq.wt) - nchar(cds.seq.alt)) %% 3 == 0) "in-frame" else "fs"
      # stop-codon creating or not
      SC <- if (grepl("\\*", aa.seq.alt)) "sc" else "no.sc"
      
      # check -------------------------------------------------------------------
      
      CHECK <- if (nchar(aa.seq.wt) > 0 && aa.seq.wt != proteins(EnsDb.Hsapiens.v86, columns = c("protein_sequence"), filter = TxIdFilter(ENST.j)) %>% as_tibble() %>% dplyr::pull(protein_sequence)) {
        "not.match"
      } else {
        "pass"
      }
      
      # check 
      print(paste0("i=", i, "; j=", j))
      print(JUNC.ID)
      print(ENST.j)
      print(c(TYPE, NOTE, FS, SC, CHECK))
      
      print("wt")
      print(exons.j %>% dplyr::select(start, end, strand, exon_idx) %>% anti_join(exons.new, by = "exon_idx"))
      print("alt")
      print(exons.new %>% dplyr::select(start, end, strand, exon_idx) %>% anti_join(exons.j, by = "exon_idx"))
      
      print("wt")
      print(aa.seq.wt)
      print("alt")
      print(aa.seq.alt)
      
      # prep for output ------------------------------------------------------------
      
      neo_ids.4 <- bind_rows(neo_ids.4, 
                             tibble(
                               junc.id = JUNC.ID, 
                               enst.model = ENST.j, 
                               aa.change = if (nchar(aa.seq.wt) > nchar(aa.seq.alt)) "loss" else "gain", 
                               aa.seq.wt = aa.seq.wt, 
                               aa.seq.alt = aa.seq.alt, 
                               cds.seq.wt = cds.seq.wt,
                               cds.seq.alt = cds.seq.alt,
                               ln.wt = nchar(aa.seq.wt), 
                               ln.alt = nchar(aa.seq.alt), 
                               ln.diff = ln.alt - ln.wt, 
                               fs = FS, 
                               sc = SC, 
                               type = TYPE, 
                               note = NOTE, 
                               check = CHECK,
                               coord_check = COORD_CHECK
                             ))
    }
  }
}

# Debug: Check if start.time exists
cat("Debug: Does start.time exist before RUNTIME?", exists("start.time"), "\n")
if (!exists("start.time")) {
  cat("Debug: Re-defining start.time as fallback\n")
  start.time <- proc.time()  # Fallback
}

RUNTIME <- tryCatch(proc.time() - start.time, error = function(e) {
  cat("Warning: Runtime calculation failed (", e$message, "). Skipping and proceeding to write file.\n")
  return(NULL)
})  # Calculate runtime before the final checks

# check -------------------------------------------------------------------

dim(neo_ids.4) %>% print()
distinct(neo_ids.4, junc.id, .keep_all = TRUE) %>% print()

# edit --------------------------------------------------------------------

neo_ids.5 <- neo_ids_final %>% 
  dplyr::select(junc.id, symbol, ensg, enst) %>% 
  right_join(neo_ids.4, by = "junc.id")

setwd(directory_11)
filename_out <- paste0("Res_AA_Prediction_Confirmed_", current_date, ".tsv")
write_tsv(neo_ids.5, filename_out, na = "NA", col_names = TRUE)

# Write the exon -> CDS/AA coordinate map (NEW). Step 12 reads this to attach
# peptide coordinates and flag n-mers that span an exon-exon boundary.
filename_coord_map <- paste0("Exon_CDS_AA_Map_", current_date, ".tsv")
write_tsv(EXON_COORD_MAP, filename_coord_map, na = "NA", col_names = TRUE)

n_mismatch <- sum(neo_ids.5$coord_check == "MISMATCH", na.rm = TRUE)
n_error    <- sum(neo_ids.5$coord_check == "ERROR", na.rm = TRUE)
cat("\nCoordinate map validation: ", sum(neo_ids.5$coord_check == "pass", na.rm = TRUE),
    " pass, ", n_mismatch, " MISMATCH, ", n_error, " ERROR (out of ", nrow(neo_ids.5), ")\n", sep = "")
if (n_mismatch > 0) {
  cat("[WARNING] ", n_mismatch, " junction(s) failed CDS reconstruction from the exon map -- ",
      "their rows are excluded from ", filename_coord_map, " and should not be trusted for peptide coordinates. ",
      "Filter neo_ids.5 (or the output tsv) on coord_check == 'MISMATCH' to inspect which junc.id/enst.model these are.\n", sep = "")
}

cat("Script completed successfully!\n")
cat("Output file:", filename_out, "\n")
if (!is.null(RUNTIME)) cat("Total runtime:", RUNTIME[3], "seconds\n")
cat("Processed", nrow(neo_ids.4), "amino acid predictions\n")
cat("Unique junctions:", nrow(distinct(neo_ids.4, junc.id)), "\n")

# Print summary of results
if(nrow(neo_ids.5) > 0) {
  cat("\nSummary of results:\n")
  cat("- Total predictions:", nrow(neo_ids.5), "\n")
  cat("- Unique junctions:", nrow(distinct(neo_ids.5, junc.id)), "\n")
  cat("- Frame-shift events:", sum(neo_ids.5$fs == "fs"), "\n")
  cat("- Stop-codon events:", sum(neo_ids.5$sc == "sc"), "\n")
  cat("- Most common types:\n")
  print(table(neo_ids.5$type) %>% sort(decreasing = TRUE) %>% head(10))
}
