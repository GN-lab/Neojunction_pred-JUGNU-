#!/bin/bash

# ------------------------- Core Configuration --------------------------

# Working directory (SSNIP root) - update if needed
export WORKDIR=""

# Input and output directories
export INPUT_DIR="${WORKDIR}/0_Input_Files"  # Purity metadata, GTF, TPM files, etc.
export OUTPUT_DIR="${WORKDIR}/results"  # Where all outputs will be saved
mkdir -p "${OUTPUT_DIR}"

# Set step-specific output directories to match OUTPUT_DIR
for i in $(seq -w 1 15); do
  export STEP${i}_OUTPUT_DIR="${OUTPUT_DIR}"
done

# Annotation and reference paths
export GTF_FILE="${INPUT_DIR}/genome.gtf"  # GTF for protein-coding extraction (Step 02)
export PURITY_FILE="${INPUT_DIR}/mapped_sample_purity.txt"  # Tumor purity metadata (Step 01)
export TPM_FILE="${INPUT_DIR}/transcript_tpm_matrix.tsv"  # Transcript-level TPM matrix (Step 03; use gene_tpm_matrix.tsv if gene-level needed)
export STAR_SJ_DIR="../output_hartwig/star_salmon/log"  # STAR SJ.out.tab files (Step 04) - fixed double slash
export GTEX_FILE="${INPUT_DIR}"  # GTEx metadata (Step 08)
export GTEX_ANNOTATION_PATH="${INPUT_DIR}"  # GTEx GTF (Step 09)
export META_FILES_PATH="${INPUT_DIR}"  # Sample TPM matrix (Step 10)
export EXTERNAL_PATH="${INPUT_DIR}"  # sjdbList.fromGTF.out.tab files (Step 11)
export PATH_TO_FASTA_UNIPROT="${INPUT_DIR}" # UP000005640_9606.fasta file (Step 12)
export HLATHENA_OUTPUT_PATH="${OUTPUT_DIR}" # results of step 12
export OPTITYPE_OUTPUT_DIR="../nextflow_HLAtyping/output_hartwig/optitype"

# Thresholds (from SSNIP README)
export MIN_PURITY=0.60  # Tumor purity threshold (Step 01)
export MIN_TPM=10  # TPM expression threshold (Step 03)
export MIN_READ_COUNT=10  # Minimum reads for junction filtering (Step 05)

# Additional filters (customizable, based on SSNIP logic)
export MIN_JUNC_OVERHANG=8  # Minimum junction overhang
export MIN_INTRON_SIZE=20  # Minimum intron size
export MAX_INTRON_SIZE=1000000  # Maximum intron size
export MIN_TOTAL_READS=20  # Cohort-wide depth
export MIN_FREQUENCY=0.01  # Minimum spliced frequency
export MIN_SAMPLES_PCT=0.10  # Minimum samples percentage for "public" junctions
export MAX_GTEX_FREQUENCY=0.01  # Maximum GTEx frequency

# R environment
export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/graichand/R_libs"  # Your R library path

# Script directories
export NEOJUNCTION_DIR="${WORKDIR}/1_Neojunction_Calling"
export NEOPEPTIDE_DIR="${WORKDIR}/2_Neopeptide_Prediction"
export PRESENTATION_DIR="${WORKDIR}/3_Presentation_Prediction"
# Path to the master pipeline driver itself (needed so child SLURM jobs can
# call back into it for validation, since SLURM stages a copy elsewhere)
export MASTER_SCRIPT_PATH="${WORKDIR}/master.sh"

# Log and checkpoint dirs
mkdir -p "${WORKDIR}/logs" "${WORKDIR}/.checkpoints"

# Other parameters
export THREADS=24  # Adjust based on SLURM cpus-per-task (from your example)

echo "[CONFIG] SSNIP config loaded from ${BASH_SOURCE[0]}"
