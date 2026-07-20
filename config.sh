#!/bin/bash

# ------------------------- Core Configuration --------------------------
# Use the main dir for all conda/mamba caches and configs
export HOME="/data/rds/DMP/UCEC/EVOLIMMU/graichand/fake_home"
export CONDA_PKGS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_cache"
export CONDA_ENVS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs"
export CONDARC="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.condarc"
export CONDA_CONFIG_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_config"
export TMPDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/tmp"

mkdir -p "$HOME" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS" "$CONDA_CONFIG_DIR" "$TMPDIR"

source "/data/scratch/DMP/UCEC/EVOLIMMU/graichand/miniconda3/etc/profile.d/conda.sh"

module load Mamba/23.1.0-0
module load CUDA/12.1.1 

# Working directory (SSNIP root) - update if needed
export WORKDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/Neojuction_pred/SSNIP"

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
export STAR_SJ_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/nextflow_rnaseq/output_hartwig/star_salmon/log"  # STAR SJ.out.tab files (Step 04) - fixed double slash
export GTEX_FILE="${INPUT_DIR}"  # GTEx metadata (Step 08)
export GTEX_ANNOTATION_PATH="${INPUT_DIR}"  # GTEx GTF (Step 09)
export META_FILES_PATH="${INPUT_DIR}"  # Sample TPM matrix (Step 10)
export EXTERNAL_PATH="${INPUT_DIR}"  # sjdbList.fromGTF.out.tab files (Step 11)
export PATH_TO_FASTA_UNIPROT="${INPUT_DIR}" # UP000005640_9606.fasta file (Step 12)
export HLATHENA_OUTPUT_PATH="${OUTPUT_DIR}" # results of step 12
export OPTITYPE_OUTPUT_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/nextflow_HLAtyping/output_hartwig/optitype"

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

# python3 and env

conda activate /data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs/neojunction_viz/

# NetMHCPan4.2
export NETMHCPAN_BIN="/data/rds/DMP/UCEC/EVOLIMMU/graichand/netMHCpan-4.2/netMHCpan"
export NETMHCPAN_INSTALL_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/netMHCpan-4.2"

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
export MIN_JUNC_OVERHANG=8  # Minimum junction overhang
export MIN_INTRON_SIZE=20  # Minimum intron size
export MAX_INTRON_SIZE=1000000  # Maximum intron size
export MIN_TOTAL_READS=20  # Cohort-wide depth
export MIN_FREQUENCY=0.01  # Minimum spliced frequency
export MIN_SAMPLES_PCT=0.10  # Minimum samples percentage for "public" junctions
export MAX_GTEX_FREQUENCY=0.01  # Maximum GTEx frequency

# R environment
export R_LIBS_USER="/data/R_libs"  # Your R library path

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

echo "[CONFIG] JUGNU config loaded from ${BASH_SOURCE[0]}"
