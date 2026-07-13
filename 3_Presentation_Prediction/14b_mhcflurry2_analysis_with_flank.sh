#!/bin/bash
# SLURM directives (submit with sbatch this_script.slurm)
#SBATCH --job-name=mhcflurry_gpu
#SBATCH --partition=gpu  # or 'gpu' if preferred
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # Match partition limits
#SBATCH --mem-per-cpu=13600
#SBATCH --gpus-per-node=1  # Request 1 GPU
#SBATCH --time=7-00:00:00  # 7 days max
#SBATCH --output=mhcflurry_%j.out
#SBATCH --error=mhcflurry_%j.err

set -euo pipefail
# Title: Step 14: MHCFlurry 2.0 Analysis (GPU-Enabled)
# November 24, 2025 | Gaurav Raichand | The Institute of Cancer Research

# Purpose: MHCFlurry 2.0 allows us to assess protein processing likelihoods for the peptides we 
# predicted in the previous steps. This algorithm offers us two experimental predictors,
# both of which are trained on mass spec-identified MHC-ligands:
#    1. Antigen Processing - models MHC allele-independent effects such as proteosomal cleavage
#    2. Presentation - integrates processing predictions with binding affinity predictions to 
#       give a composite "presentation score".

# Use the main dir for all conda/mamba caches and configs
export HOME="/data/rds/DMP/UCEC/EVOLIMMU/graichand/fake_home"
export CONDA_PKGS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_cache"
export CONDA_ENVS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs"
export CONDARC="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.condarc"
export CONDA_CONFIG_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_config"
export TMPDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/tmp"

mkdir -p "$HOME" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS" "$CONDA_CONFIG_DIR" "$TMPDIR"

source "/data/scratch/DMP/UCEC/EVOLIMMU/csalas/miniconda3/etc/profile.d/conda.sh"

module load Mamba/23.1.0-0
module load CUDA/12.1.1  # Matches your available module; bridges with Conda cudatoolkit=11.7.0

################################################################
# 0. (DISABLED) Install MHCFlurry environment and models
#    Kept for future re-install; currently commented via here-doc
################################################################

: <<'INSTALL_MHCFLURRY'
# Optional: clean any previous env to avoid PATH/TensorFlow conflicts (non-fatal if not present)
conda deactivate || true
mamba env remove -n mhcflurry-env -y || true

# Configure channels with strict priority for reliable Bioconda resolution on HPC
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nvidia  # For any needed CUDA fallbacks
conda config --set channel_priority strict

# Switch to classic solver to resolve libmamba strict priority warnings
conda config --set solver classic

# Disable mamba's low-speed download timeout (common on proxies/slow HPC networks)
export MAMBA_NO_LOW_SPEED_LIMIT=1

# Create base env with MHCflurry (no explicit CUDA here—handled by pip below)
mamba create -n mhcflurry-env -c conda-forge -c bioconda -y python=3.10 mhcflurry=2.1.5
conda activate mhcflurry-env

# Install TensorFlow with GPU support via pip (auto-resolves CUDA/cuDNN)
pip install --upgrade pip
pip install tensorflow[and-cuda] nvidia-cudnn-cu12  # Explicit cuDNN for CUDA 12.x

# Enable TensorFlow GPU memory growth (prevents full allocation) and test
python - <<'EOF'
import tensorflow as tf
print("TensorFlow version:", tf.__version__)
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        logical_gpus = tf.config.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
        print(e)
else:
    print("No GPUs detected—check setup.")
EOF

# Download datasets and trained models (run once—fetches ~1GB of data)
mhcflurry-downloads fetch

################################################################
# 1. Download MHCFlurry models (specific to presentation if needed)
################################################################

# models_class1_presentation is usually what we need as the presentation predictor includes
# a peptide/MHC I binding affinity predictor and an antigen processing predictor
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads path models_class1_presentation  # This gets the path to the downloaded files
INSTALL_MHCFLURRY

################################################################
# 0b. Runtime: activate existing env and check OUTPUT_DIR
################################################################

conda deactivate || true
conda activate mhcflurry-env

# Require OUTPUT_DIR (and thus STEP14_OUTPUT_DIR via config.sh) to be set
: "${OUTPUT_DIR:?OUTPUT_DIR is not set. Please source config.sh before sbatch.}"

# Use OUTPUT_DIR for this step (equivalent to STEP14_OUTPUT_DIR from config.sh)
STEP14_OUTPUT_DIR="${OUTPUT_DIR}"

################################################################
# 2. Generate Predictions
################################################################

# Hardcode the date matching your input CSVs (from R script)
# If your R script writes *_optitype_samples.csv, adjust RUN_DATE or filenames accordingly.
RUN_DATE="2026_1124"

# Run MHCflurry predictions for each input CSV (no --alleles since CSV has 'allele' column)

# 8 mers
mhcflurry-predict \
  "${STEP14_OUTPUT_DIR}/08mer_mhcflurry_input_${RUN_DATE}.csv" \
  --out "${STEP14_OUTPUT_DIR}/08mers_flank_mhcflurry_${RUN_DATE}.csv"

# 9 mers
mhcflurry-predict \
  "${STEP14_OUTPUT_DIR}/09mer_mhcflurry_input_${RUN_DATE}.csv" \
  --out "${STEP14_OUTPUT_DIR}/09mers_flank_mhcflurry_${RUN_DATE}.csv"

# 10 mers
mhcflurry-predict \
  "${STEP14_OUTPUT_DIR}/10mer_mhcflurry_input_${RUN_DATE}.csv" \
  --out "${STEP14_OUTPUT_DIR}/10mers_flank_mhcflurry_${RUN_DATE}.csv"

# 11 mers
mhcflurry-predict \
  "${STEP14_OUTPUT_DIR}/11mer_mhcflurry_input_${RUN_DATE}.csv" \
  --out "${STEP14_OUTPUT_DIR}/11mers_flank_mhcflurry_${RUN_DATE}.csv"

echo "MHCflurry predictions complete. Outputs saved in ${STEP14_OUTPUT_DIR}."

################################################################
# 3. Generate Predictions - WT/self side (NEW)
#    Feeds Step 15e, which removes any ALT neoantigen (peptide + allele)
#    that also arises from these WT/self peptides.
#
#    DEDUP BEFORE PREDICT: the raw WT input files have ~196 million rows
#    but only ~371,000 unique (peptide, allele, n_flank, c_flank)
#    combinations -- 99.8% duplication. MHCflurry gives identical results
#    for identical inputs, so predicting duplicates wastes ~99.8% of GPU
#    time. We deduplicate first (seconds), predict on the unique set
#    (~2h total vs 20h+ per length), then the full 196M-row file is not
#    needed again -- Step 15e only needs the unique scored peptide+allele
#    pairs to do its filtering, so no re-expansion step is required.
################################################################

dedup_and_predict() {
  local IN_FILE="$1"
  local OUT_FILE="$2"
  local DEDUP_FILE="${IN_FILE%.csv}_dedup.csv"

  echo "Deduplicating ${IN_FILE}..."
  awk -F',' '
    NR==1 { print; next }
    !seen[$0]++
  ' "$IN_FILE" > "$DEDUP_FILE"

  local before after
  before=$(wc -l < "$IN_FILE")
  after=$(wc -l < "$DEDUP_FILE")
  echo "  Rows: ${before} -> ${after} after dedup (removed $((before - after)) duplicates)"

  echo "Running MHCflurry on deduplicated file..."
  mhcflurry-predict "$DEDUP_FILE" --out "$OUT_FILE"

  rm -f "$DEDUP_FILE"
  echo "Done: $OUT_FILE"
}

# 8 mers WT
dedup_and_predict \
  "${STEP14_OUTPUT_DIR}/08mer_mhcflurry_input_wt_${RUN_DATE}.csv" \
  "${STEP14_OUTPUT_DIR}/08mers_flank_mhcflurry_wt_${RUN_DATE}.csv"

# 9 mers WT
dedup_and_predict \
  "${STEP14_OUTPUT_DIR}/09mer_mhcflurry_input_wt_${RUN_DATE}.csv" \
  "${STEP14_OUTPUT_DIR}/09mers_flank_mhcflurry_wt_${RUN_DATE}.csv"

# 10 mers WT
dedup_and_predict \
  "${STEP14_OUTPUT_DIR}/10mer_mhcflurry_input_wt_${RUN_DATE}.csv" \
  "${STEP14_OUTPUT_DIR}/10mers_flank_mhcflurry_wt_${RUN_DATE}.csv"

# 11 mers WT
dedup_and_predict \
  "${STEP14_OUTPUT_DIR}/11mer_mhcflurry_input_wt_${RUN_DATE}.csv" \
  "${STEP14_OUTPUT_DIR}/11mers_flank_mhcflurry_wt_${RUN_DATE}.csv"

echo "MHCflurry WT/self predictions complete. Outputs saved in ${STEP14_OUTPUT_DIR}."
