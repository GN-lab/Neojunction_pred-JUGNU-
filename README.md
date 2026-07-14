# JUGNU (Junction Guided Neoantigen Uncoverer)

<img width="1536" height="1024" alt="Designer-4" src="https://github.com/user-attachments/assets/56a2ac52-43d4-4b80-9c89-b05997d9ea4b" />

---

## Overview

JUGNU (Junction-guided Neoantigen Uncoverer) is a computational pipeline for identifying tumour-specific aberrant splicing events and predicting splice-derived neoantigens from bulk RNA-sequencing data.

The framework integrates splice-junction discovery, intron retention detection, transcript reconstruction, peptide generation, HLA typing, and neoantigen prediction into a unified workflow.

---

## Workflow

The pipeline consists of three major stages:

### Aberrant Splicing Discovery

RNA-seq reads are processed using STAR and IRFinder.

Two complementary event types are detected:

- Novel splice junctions (SJ.out.tab)
- Intron retention events (IRFinder)

Candidate events are filtered against:

- Reference GTF annotation
- GTEx normal tissue splicing database
- Read support thresholds
- Positive Sample Rate (PSR)

Only high-confidence tumour-specific events are retained.

---

### Transcript Reconstruction & Translation

Retained events are reconstructed into alternative transcripts.

Supported event classes include:

- Exon Skipping (ES)
- Alternative 3' Splice Site (A3SS)
- Alternative 5' Splice Site (A5SS)
- Exon Gain/Loss
- Intron Retention (IR)

Alternative transcripts are translated using a multi-frame strategy to recover potential novel open reading frames.

---

### Neoantigen Prediction

Altered peptide sequences are generated and evaluated against patient-specific HLA alleles.

Predictions are performed using:

- MHCFlurry 2.0

Reported outputs include:

- Binding affinity
- Antigen processing
- Presentation probability

High-confidence peptide-HLA pairs are prioritised as candidate neoantigens.

---

## Inputs

Required inputs:

```text
samples.txt
FASTQ files
Reference genome
GTF annotation
GTEx junction database
Tumour metadata
HLA typing information
```

---

## Outputs

JUGNU generates:

```text
Novel junction calls
Intron retention events
Alternative transcripts
Protein sequences
Candidate peptides
Neoantigen predictions
Sample-level neoantigen matrices
```

---

## Software Dependencies

```text
Nextflow
STAR
Salmon
IRFinder-S
Python >= 3.10
MHCFlurry 2.0
R >= 4.3
```

---

## STAR Parameters

Splice junction discovery is performed using STAR with the following settings:

```bash
--star_align_args "--outSAMtype BAM Unsorted \
                     --twopassMode Basic \
                     --outFilterMultimapScoreRange 1 \
                     --outFilterMultimapNmax 20 \
                     --outFilterMismatchNmax 10 \
                     --alignIntronMax 500000 \
                     --alignMatesGapMax 1000000 \
                     --sjdbScore 2 \
                     --alignSJDBoverhangMin 1 \
                     --genomeLoad NoSharedMemory \
                     --limitBAMsortRAM 80000000000 \
                     --readFilesCommand gunzip -c \
                     --outFilterMatchNminOverLread 0.33 \
                     --outFilterScoreMinOverLread 0.33 \
                     --sjdbOverhang 100 \
                     --outSAMstrandField intronMotif \
                     --outSAMattributes NH HI NM MD AS XS \
                     --limitSjdbInsertNsj 2000000 \
                     --outSAMunmapped Within \
                     --outSAMheaderHD @HD VN:1.4 \
                     --outSAMmultNmax 1" \
```

These parameters increase splice-junction confidence by enforcing minimum read support and biologically plausible intron lengths.

---

## Applications

JUGNU can be applied to:

- Cancer neoantigen discovery
- Alternative splicing analysis
- Intron retention studies
- Tumour-specific transcript discovery
- Immunotherapy target identification
- Precision oncology research

---

## Contact

Gaurav Nagar 

Institute of Cancer Research (ICR), London

For bug reports, feature requests, or collaborations, please open a GitHub issue.
