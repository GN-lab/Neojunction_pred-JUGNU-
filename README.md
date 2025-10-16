# Neojunction_prediction

Neojunction prediction from 307 transcriptomic dataset acquired from Hartwig (Breast sample)

______________________________________________________________________________________________________
Phase 1: Neojunction prediction
______________________________________________________________________________________________________

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

______________________________________________________________________________________________________
Phase 2: Sequence Processing Pipeline
______________________________________________________________________________________________________

<img width="3840" height="2760" alt="phase2_flowchart" src="https://github.com/user-attachments/assets/abef3a8b-758e-4231-afbe-dca1c4b94b68" />

***This flowchart explains the process in detail:
1. Translate neojunctions into altered amino acid sequences: For each retained neojunction from Phase 1, translate the aberrant spliced transcripts into amino acid (AA) sequences using a three-frame translation approach (forward, +1, and +2 frames) to account for potential novel open reading frames (ORFs); this captures all possible protein-level alterations arising from the non-canonical splicing events (yielding 2,237 altered AA sequences from 1,513 neojunctions across 666 unique genes). Translation begins at the first in-frame start codon downstream of the junction or the junction itself if no codon is available, ensuring comprehensive coverage of potential neoantigens while discarding frames without viable ORFs.
2. Extract candidate peptides (N-mers) from altered sequences: Apply a sliding window to generate all possible 8- to 11-mer peptides from the translated AA sequences, as these lengths are optimal for binding to HLA class I molecules (major histocompatibility complex class I proteins on cell surfaces that present intracellular peptides to T cells); focus on peptides spanning the neojunction to maximize immunogenicity (producing 580,455 peptides, of which 79,114 are unique), with each peptide tagged to its originating neojunction, gene, and flanking residues for subsequent HLA binding predictions. This step prioritizes HLA class I-restricted epitopes, which are crucial for cytotoxic T-cell recognition in cancer immunotherapy.

______________________________________________________________________________________________________
Phase 3: Neoantigen Prediction Pipeline
______________________________________________________________________________________________________

<img width="3840" height="3300" alt="phase3_flowchart" src="https://github.com/user-attachments/assets/b446eb8c-a7f6-41fc-8f52-106f94ef81e8" />

***This flowchart explains the process in detail:
1. Prepare input for MHCflurry binding predictions: Construct a comprehensive input DataFrame by combining the generated 8-11mer peptides with their 3-residue N-terminal and C-terminal flanks (contextual sequences that influence MHC processing and presentation) and the patient's HLA alleles (human leukocyte antigen types, polymorphic proteins that dictate peptide presentation; here, using up to 5 alleles per patient based on prior selection criteria); this results in 580,455 rows ready for batch prediction, ensuring each peptide is evaluated against the individual's HLA repertoire to identify personalized neoantigens (79,114 unique peptides prepared). The flanks are essential as MHCflurry's model incorporates proteasomal cleavage and TAP transport signals for realistic immunogenicity scoring.
2. Run MHCflurry for affinity, processing, and presentation predictions: Execute the MHCflurry algorithm—a machine learning tool trained on mass spectrometry and binding assay data—to predict each peptide-HLA pair's binding affinity (IC50 value in nanomolar, where <500 nM indicates strong binders), proteasomal processing efficiency (cleavage likelihood at flanks), and overall presentation score (composite metric integrating binding, processing, and stability); process across all patient HLAs in parallel (yielding 342,344 total predictions for the cohort, covering 79,114 unique peptides against 5 common alleles). MHCflurry outperforms traditional netMHC tools by leveraging pan-allele models, providing robust scores for rare HLAs and reducing false positives in neoantigen discovery.
3. Construct neoantigen presence/absence matrix: Aggregate the MHCflurry outputs into a binary matrix indicating neoantigen presence (peptides passing a presentation score threshold, e.g., top percentile or <500 nM binding with favorable processing) across 307 patient samples, linking back to the 945 unique neojunctions; include binding scores, HLA matches, and per-sample counts to enable downstream analyses like tumor specificity and immunogenicity ranking (42,157 putative neoantigens identified). This matrix facilitates cohort-level insights, such as neoantigen burden per patient and shared epitopes, while flagging high-confidence candidates for validation in T-cell assays or vaccine design.

______________________________________________________________________________________________________
SJ.out.tab file generating prameters used :-
______________________________________________________________________________________________________

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

______________________________________________________________________________________________________
Results :-
______________________________________________________________________________________________________

<img width="2400" height="1800" alt="figure_5i_fs_if_boxplot_20250928" src="https://github.com/user-attachments/assets/afada8f3-e0b4-4ad9-8843-6276ad1510eb" />

<img width="3000" height="1800" alt="figure_5i_splice_types_jitter_20250928" src="https://github.com/user-attachments/assets/cdb3c880-395a-4753-9251-12ddf68b5f1f" />

<img width="2400" height="1500" alt="figure_5i_fs_if_density_all_20250928" src="https://github.com/user-attachments/assets/1b15fbf6-babc-4fb4-bcd7-caaca6ac60d3" />

<img width="2400" height="3000" alt="figure_5i_splice_types_density_all_20250928" src="https://github.com/user-attachments/assets/73be694c-fd42-4bdf-89d6-a7c3983c3f9f" />

<img width="2400" height="1800" alt="figure_5i_fs_if_jitter_20250928" src="https://github.com/user-attachments/assets/40a758b0-f2d4-464c-b5f4-51e530a5b98b" />

<img width="3600" height="1800" alt="figure_sample_neoantigens_20250929" src="https://github.com/user-attachments/assets/0271d17b-6499-4232-8e17-b50e91b3a1b2" />

