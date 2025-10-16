# Neojunction_prediction

Neojunction prediction from 307 transcriptomic dataset acquired from Hartwig (Breast sample)

<img width="4170" height="7770" alt="neojunction_flowchart_clean" src="https://github.com/user-attachments/assets/5ab05221-1258-40fa-bcc4-8da275d0d7f4" />

***This flowchart explains the process in detail:
1. Select high-purity tumors: Apply a tumor purity threshold of ≥0.60 to obtain 40 high-purity tumors, which are used for all downstream counts and PSR calculations. (Positive Sample Rate (PSR): fraction of tumor samples with a given junction supported by ≥10 reads; for 40 samples, PSR = number of positive samples ÷ 40 (threshold used: PSR ≥0.1 for recurrence).
2. Extract normal spliced junctions: Extract those spliced junctions which are considered as normal through the GTF file (150k+ transcripts, 20k gene).
3. Filter for highly expressed transcripts: Filter out those with highly expressed transcripts (e.g., only those transcripts which have >=10 TPM) on the basis of their read counts. Hence, making sure that we are choosing only those junctions where the most transcriptions happened, so downstream calls occur in active genes (2,828 transcripts, 2,206 genes).
4. Derive annotated splice junctions: Derive annotated (normal) splice junctions by connecting consecutive exons from those expressed transcripts; this creates the reference blacklist of normal junctions (≈41k total, 2,737 unique, 2.2k genes) to subtract later.
5. Process tumor splice-junctions: Load nfcore-rna seq generated SJ.out.tab files and subtract the annotated set; retain non-annotated candidates that have ≥10 supporting reads in at least one sample and fall in protein-coding regions (16,608 candidates; ≈1.6% novel rate).
6. Apply PSR threshold for recurrence: Keep junctions present in ≥0.1 PSR of tumor samples (≥4/40 with ≥10 reads), which reduces noise and focuses on shared tumor splicing events (3,368 retained junctions).
7. Map to genes/transcripts: Map each retained junction to genes/transcripts (overlap with GTF features) ≈3k unique overlaps across ≈300 genes.
8. Apply reliability filtering: Apply reliability (judge) filtering that requires sufficient read count, splice-site depth, and junction frequency; only junctions that pass all criteria are retained (1,541 survivors across ≈200 genes).
9. Confirm tumor specificity: Confirm tumor specificity by requiring negligible signal in normal tissues (GTEx PSR ≤0.01; in your run all survivors had 0), yielding a tumor-only set (still 1,541). GTEx is the NIH Genotype‑Tissue Expression program that profiles gene expression and splicing across many non‑diseased human tissues from hundreds of donors to link genetic variation to tissue‑specific RNA patterns.
10. Compile final neojunctions: Compile final neojunctions with per-sample counts, PSR, depth/frequency, and a composite “judge” pass flag—ready for neoantigen prediction and immunogenicity ranking (1,541 total).


**STAR aligner** 

STAR produces SJ.out.tab: Contains every splice junction STAR detected, with per-junction read counts (unique + multi-mapping reads) and annotation flags.

- srun nextflow run nf-core/rnaseq \
 
  --aligner star_salmon \
  --skip_quantification \
  --skip_qc \
   --star_align_args "--outSAMtype BAM SortedByCoordinate \
                     --quantMode TranscriptomeSAM \
                     --outSAMunmapped Within \
                     --outSJfilterOverhangMin 15 20 20 20 \
                     --alignSJoverhangMin 8 \
                     --alignIntronMin 20 \
                     --alignIntronMax 1000000" \
  -profile singularity \
  -c /data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/config_files/icr_alma.config \
  -r 3.18.0 \
  -resume

--outSJfilterOverhangMin 15 20 20 20
Controls the minimum overhang (length of mapped sequence flanking a splice junction) for different splice junction types to be reported. This increases stringency on junction detection:
  -  15 for canonical GT/AG junctions
  -  20 for other types (e.g., GC/AG, AT/AC, non-canonical)
This helps reduce false positives and ensures only well-supported junctions are output.

--alignSJoverhangMin 8
Sets the minimum overhang length required on each side of a splice junction to consider it a valid splice junction during alignment. A value of 8 bases means STAR must see at least 8 matching bases flanking the junction.

--alignIntronMin 20
Specifies the minimum allowed intron length (20 bases) for splice junctions. Very short introns are usually sequencing or alignment artifacts, so this filters those out.

--alignIntronMax 1000000
Sets the maximum allowed intron length (1,000,000 bases). This limits extremely long introns that could be biologically implausible or mapping artifacts.

<img width="2400" height="1800" alt="figure_5i_fs_if_boxplot_20250928" src="https://github.com/user-attachments/assets/afada8f3-e0b4-4ad9-8843-6276ad1510eb" />

<img width="3000" height="1800" alt="figure_5i_splice_types_jitter_20250928" src="https://github.com/user-attachments/assets/cdb3c880-395a-4753-9251-12ddf68b5f1f" />

<img width="2400" height="1500" alt="figure_5i_fs_if_density_all_20250928" src="https://github.com/user-attachments/assets/1b15fbf6-babc-4fb4-bcd7-caaca6ac60d3" />

<img width="2400" height="3000" alt="figure_5i_splice_types_density_all_20250928" src="https://github.com/user-attachments/assets/73be694c-fd42-4bdf-89d6-a7c3983c3f9f" />

<img width="2400" height="1800" alt="figure_5i_fs_if_jitter_20250928" src="https://github.com/user-attachments/assets/40a758b0-f2d4-464c-b5f4-51e530a5b98b" />

<img width="3600" height="1800" alt="figure_sample_neoantigens_20250929" src="https://github.com/user-attachments/assets/0271d17b-6499-4232-8e17-b50e91b3a1b2" />

