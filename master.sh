#!/bin/bash
#SBATCH --job-name=JUGNU_driver
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8042
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/master_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "**ERROR** rc=$rc at line $LINENO cmd: $BASH_COMMAND" >&2; exit $rc' ERR

WORKDIR="$(pwd)"

# SCRIPT_PATH must be the DURABLE path to this file, not derived from $0.
# When you `sbatch master.sh`, SLURM copies this file into a temporary
# per-job spool dir on the compute node (that's what
# /tmp/slurmd/job<id>/slurm_script is) and runs THAT copy. $0 inside the
# running script resolves to that ephemeral spool path, not your real file
# -- so every child job's callback into "$SCRIPT_PATH __validate_step__ ..."
# needs a path that survives after the (few-second) submitting job has
# already finished and its spool dir is gone.
#
# Rather than hardcoding a filename here (which breaks the moment you
# rename/move this file, as just happened), MASTER_SCRIPT_PATH is read from
# config.sh -- set it there once, e.g.:
#   export MASTER_SCRIPT_PATH="/data/rds/DMP/UCEC/EVOLIMMU/graichand/Neojuction_pred/SSNIP/master.sh"
# and you never need to touch this file again if you rename/move it.
if [[ -z "${MASTER_SCRIPT_PATH:-}" ]]; then
  # config.sh may not be sourced yet this early -- source it now just to
  # pick up MASTER_SCRIPT_PATH, then again later in normal flow (harmless).
  [[ -f "${WORKDIR}/config.sh" ]] && source "${WORKDIR}/config.sh"
fi
SCRIPT_PATH="${MASTER_SCRIPT_PATH:-${WORKDIR}/master.sh}"
if [[ ! -f "$SCRIPT_PATH" ]]; then
  echo "[ERROR] Expected to find this script at $SCRIPT_PATH but it's not there." >&2
  echo "[ERROR] Add 'export MASTER_SCRIPT_PATH=/full/path/to/this/script.sh' to config.sh." >&2
  exit 1
fi

###############################################################################
#  SELF-DISPATCH: this same file is what runs *inside* every SLURM step's
#  --wrap command to check the step's output before the checkpoint is
#  touched. Keeps everything in one script instead of a separate validator
#  file. Normal `sbatch master.sh` submissions never hit
#  these branches -- only `bash "$SCRIPT_PATH" __validate_*__ ...` calls do.
###############################################################################

_v_single() {
  # Generic check: newest file matching a pattern exists, is >50 bytes, and
  # (for .tsv/.csv/.txt) has more than just a header row.
  local step_name="$1" search_dir="$2" pattern="$3" min_lines="${4:-2}"

  if [[ "$search_dir" == "SKIP" || "$pattern" == "SKIP" ]]; then
    echo "[VALIDATE SKIP] $step_name: no pattern configured -- fill in the phase table below." >&2
    return 0
  fi
  if [[ ! -d "$search_dir" ]]; then
    echo "[VALIDATE FAIL] $step_name: directory does not exist: $search_dir" >&2
    return 1
  fi

  local match
  match=$(find "$search_dir" -maxdepth 1 -type f -regextype posix-extended -regex ".*/${pattern}" 2>/dev/null | sort | tail -n 1)
  if [[ -z "$match" ]]; then
    echo "[VALIDATE FAIL] $step_name: no file matching '$pattern' in $search_dir" >&2
    return 1
  fi

  local size
  size=$(stat -c%s "$match" 2>/dev/null || echo 0)
  if [[ "$size" -lt 50 ]]; then
    echo "[VALIDATE FAIL] $step_name: $match is essentially empty (${size} bytes)" >&2
    return 1
  fi

  case "$match" in
    *.tsv|*.csv|*.txt)
      local nlines
      nlines=$(wc -l < "$match")
      if [[ "$nlines" -lt "$min_lines" ]]; then
        echo "[VALIDATE FAIL] $step_name: $match has only $nlines line(s) (need >= $min_lines) -- looks header-only/empty" >&2
        return 1
      fi
      ;;
  esac

  echo "[VALIDATE OK] $step_name: $match (${size} bytes)"
  return 0
}

_v_dir_nonempty() {
  # Directory-based check (for 05b's IRFinder reference build: many files, no single output)
  local step_name="$1" dir="$2" min_files="${3:-5}"
  if [[ ! -d "$dir" ]]; then
    echo "[VALIDATE FAIL] $step_name: directory missing: $dir" >&2
    return 1
  fi
  local n
  n=$(find "$dir" -type f | wc -l)
  if [[ "$n" -lt "$min_files" ]]; then
    echo "[VALIDATE FAIL] $step_name: only $n file(s) in $dir, expected >= $min_files" >&2
    return 1
  fi
  echo "[VALIDATE OK] $step_name: $n files in $dir"
  return 0
}

_v_irfinder_bulk() {
  # Per-sample check (for 05c, which tolerates individual sample failures by
  # design -- job exit code alone can't tell you if most samples actually ran)
  local step_name="$1" irfinder_dir="$2" samples_file="$3" min_ratio="${4:-0.9}"
  if [[ ! -f "$samples_file" ]]; then
    echo "[VALIDATE FAIL] $step_name: samples file missing: $samples_file" >&2
    return 1
  fi
  local total=0 ok=0 sample f
  while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue
    total=$((total+1))
    f="${irfinder_dir}/${sample}/IRFinder-IR-nondir.txt"
    [[ -s "$f" ]] && ok=$((ok+1))
  done < "$samples_file"
  if [[ "$total" -eq 0 ]]; then
    echo "[VALIDATE FAIL] $step_name: samples file was empty" >&2
    return 1
  fi
  local ratio pass
  ratio=$(awk -v ok="$ok" -v total="$total" 'BEGIN{printf "%.3f", ok/total}')
  pass=$(awk -v r="$ratio" -v m="$min_ratio" 'BEGIN{print (r>=m)?1:0}')
  if [[ "$pass" -ne 1 ]]; then
    echo "[VALIDATE FAIL] $step_name: only $ok/$total samples have non-empty IRFinder output (ratio=$ratio, need >=$min_ratio)" >&2
    return 1
  fi
  echo "[VALIDATE OK] $step_name: $ok/$total samples have IRFinder output (ratio=$ratio)"
  return 0
}

_v_phase() {
  local phase="$1"
  local fail=0

  check()      { _v_single      "$1" "$2" "$3" "${4:-2}"     || fail=1; }
  check_dir()  { _v_dir_nonempty "$1" "$2" "${3:-5}"          || fail=1; }
  check_irf()  { _v_irfinder_bulk "$1" "$2" "$3" "${4:-0.9}"  || fail=1; }

  if [[ "$phase" == "phase1" ]]; then
    check     "01_tumor_purity"                    "${OUTPUT_DIR:-SKIP}" 'Patient_List_Post_TumorPurity_Filter_[0-9]+\.[0-9]{2}\.txt'
    check     "02_protein_coding_genes"             "${OUTPUT_DIR:-SKIP}" 'GTF_Protein_Coding_Genes_[0-9]{8}\.txt'
    check     "03_tpm_filter_10"                    "${OUTPUT_DIR:-SKIP}" 'GTF_ProteinCoding_Filter[0-9]+_[0-9]{8}\.tsv'
    check     "04_extract_annot_sj_to_analyze"      "${OUTPUT_DIR:-SKIP}" 'SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_[0-9]{8}\.tsv'
    check     "05a_extract_sj_nonannot_hartwig"     "${OUTPUT_DIR:-SKIP}" 'SJ_List_NonAnnotated_Candidates_Protein_Coding_[0-9]{8}\.tsv'
    check_dir "05b_build_irfinder_reference"        "${WORKDIR}/0_Input_Files/irfinder_reference_grch38" 5
    check_irf "05c_run_irfinder_hartwig_samples"    "${WORKDIR}/results/irfinder_hartwig_results" "${WORKDIR}/samples.txt" 0.9
    check     "05d_aggregate_irfinder_results"      "${OUTPUT_DIR:-SKIP}" 'IR_Candidates_PSR10_High_Confidence\.tsv'
    check     "05e_extract_sj_nonannot_psr10"       "${OUTPUT_DIR:-SKIP}" 'SJ_PSR_NonAnnotated_Candidates_Protein_Coding_[0-9]{8}\.tsv'
    check     "05f_merge_SJ_IR_junctions"           "${OUTPUT_DIR:-SKIP}" 'SJ_Novel_Potential_Neojunctions\.tsv'
    check     "06_prep_overlap_table"               "${OUTPUT_DIR:-SKIP}" 'SJ_Overlap_Table_[0-9]{8}\.tsv'
    check     "07_count_depth_freq_judge"           "${OUTPUT_DIR:-SKIP}" 'PSR_Table_[0-9]{8}\.tsv'
    check     "08_overlap_table_gtex"               "${OUTPUT_DIR:-SKIP}" 'GTEx_SJ_Overlap_Table_[0-9]{8}\.tsv'
    check     "09_count_depth_freq_judge_psr_gtex"  "${OUTPUT_DIR:-SKIP}" 'PSR_Retained_Passed_GTEx_[0-9]{8}\.tsv'
    check     "10_extract_neojunctions"             "${OUTPUT_DIR:-SKIP}" 'PSR_Neojunctions_[0-9]{8}\.tsv'
  elif [[ "$phase" == "phase2" ]]; then
    check "11_aaseq_prediction"               "${STEP11_OUTPUT_DIR:-SKIP}" 'Res_AA_Prediction_Confirmed_[0-9]{8}\.tsv'
    check "12_nmer_generation"                "${OUTPUT_DIR:-SKIP}"        '2023_0812_complete_list_all_mers\.tsv'
    check "12_nmer_generation_wt"             "${OUTPUT_DIR:-SKIP}"        '2023_0812_hlathenalist_msic_08mers_wt\.tsv'
    check "14a_mhcflurry2_input_alt"          "${STEP14_OUTPUT_DIR:-SKIP}" '08mer_mhcflurry_input_[0-9_]+\.csv'
    check "14a_mhcflurry2_input_wt"           "${STEP14_OUTPUT_DIR:-SKIP}" '08mer_mhcflurry_input_wt_[0-9_]+\.csv'
    check "14b_mhcflurry2_analysis_alt"       "${STEP14_OUTPUT_DIR:-SKIP}" '08mers_flank_mhcflurry_[0-9_]+\.csv'
    check "14b_mhcflurry2_analysis_wt"        "${STEP14_OUTPUT_DIR:-SKIP}" '08mers_flank_mhcflurry_wt_[0-9_]+\.csv'
    check "14c_select_top_alleles"            "${STEP14_OUTPUT_DIR:-SKIP}" 'mhcflurry_08mer_top_[0-9]{8}\.tsv'
    check "15d_scores_for_nj_types"           "${STEP15_OUTPUT_DIR:-SKIP}" 'top_neoantigens_matrix_[0-9]{8}\.tsv'
    check "15e_filter_alt_neoantigens_vs_wt"  "${STEP15_OUTPUT_DIR:-SKIP}" 'tumor_specific_neoantigens_wt_filtered_[0-9]{8}\.tsv'
  else
    echo "[VALIDATE FAIL] unknown phase '$phase'" >&2
    return 1
  fi

  if [[ "$fail" -eq 1 ]]; then
    echo "[PHASE VALIDATE FAIL] $phase: one or more steps failed the data-presence check -- see [VALIDATE FAIL] lines above." >&2
    return 1
  fi
  echo "[PHASE VALIDATE OK] $phase"
  return 0
}

if [[ "${1:-}" == "__validate_step__" ]]; then
  shift
  source "${WORKDIR}/config.sh"
  _v_single "$@"
  exit $?
elif [[ "${1:-}" == "__validate_dir__" ]]; then
  shift
  source "${WORKDIR}/config.sh"
  _v_dir_nonempty "$@"
  exit $?
elif [[ "${1:-}" == "__validate_irfinder__" ]]; then
  shift
  source "${WORKDIR}/config.sh"
  _v_irfinder_bulk "$@"
  exit $?
elif [[ "${1:-}" == "__validate_phase__" ]]; then
  shift
  source "${WORKDIR}/config.sh"
  _v_phase "$@"
  exit $?
fi

###############################################################################
#  NORMAL MODE: submit the pipeline, phase-gated, dependency-chained.
###############################################################################

source "${WORKDIR}/config.sh"

if [[ ! -s "${WORKDIR}/samples.txt" ]]; then
  echo "[INFO] Generating samples.txt from STAR SJ.out.tab files..."
  find "${STAR_SJ_DIR}" -type f -name "*.SJ.out.tab" | \
    awk -F'/' '{fname=$NF; sub(/\.SJ\.out\.tab$/, "", fname); print fname}' | \
    sort -u > "${WORKDIR}/samples.txt"
  echo "[INFO] Created samples.txt with $(wc -l < "${WORKDIR}/samples.txt") samples."
fi

export SAMPLES_TXT="${WORKDIR}/samples.txt"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$THREADS"

module load R
module load Mamba/23.1.0-0
source /data/scratch/DMP/UCEC/EVOLIMMU/graichand/miniconda3/etc/profile.d/conda.sh
conda activate /data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs/neojunction_viz/

echo "=== Submitting JUGNU Pipeline (dependency-driven, 2 phases, self-validating) ==="

for f in "${GTF_FILE}" "${PURITY_FILE}" "${TPM_FILE}" "${STAR_SJ_DIR}" "${GTEX_FILE}"; do
  [[ -e "$f" ]] || { echo "[ERROR] Missing required file/directory: $f" >&2; exit 1; }
done

mkdir -p "${WORKDIR}/.checkpoints" "${WORKDIR}/logs"

step_done() { [[ -f "${WORKDIR}/.checkpoints/$1.done" ]]; }

jobPrev=""

###############################################################################
#  Durable job scripts, not --wrap
#
#  --wrap tells sbatch to stage your command as a script under the compute
#  node's LOCAL /tmp (e.g. /tmp/slurmd/job<id>/slurm_script). On shared HPC
#  that path can get swept by a /tmp cleanup daemon or lost to a node hiccup
#  mid-job, which is what caused:
#      bash: /tmp/slurmd/job22018213/slurm_script: No such file or directory
#  even though the R script itself had already finished successfully.
#
#  write_job_script() instead writes the real script to $WORKDIR/.slurm_jobs/
#  (durable, shared storage), and sbatch submits that file directly. It also
#  means you can `cat .slurm_jobs/<name>.sh` any time to see exactly what a
#  step ran.
###############################################################################

JOBSCRIPT_DIR="${WORKDIR}/.slurm_jobs"
mkdir -p "$JOBSCRIPT_DIR"

write_job_script() {
  local name="$1" body="$2"
  local script_path="${JOBSCRIPT_DIR}/${name}.sh"
  cat > "$script_path" <<EOF
#!/bin/bash
set -euo pipefail
${body}
EOF
  chmod +x "$script_path"
  echo "$script_path"
}

# submit_job: writes a durable script that runs <command>, then self-validates,
# then touches its checkpoint (set -e means any failed line stops the chain).
# Returns (via stdout, ONLY the job id) so callers can chain --dependency.
submit_job() {
  local name="$1" command="$2" dep="${3:-}"
  local val_dir="${4:-SKIP}" val_pattern="${5:-SKIP}" val_minlines="${6:-2}"

  local body="$command
bash \"${SCRIPT_PATH}\" __validate_step__ \"$name\" \"$val_dir\" \"$val_pattern\" \"$val_minlines\"
touch \"${WORKDIR}/.checkpoints/${name}.done\""

  local script_path
  script_path=$(write_job_script "$name" "$body")

  local jobid
  if [[ -n "$dep" ]]; then
    jobid=$(sbatch --parsable --job-name="$name" --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=12:00:00 --output=logs/${name}_%j.log --dependency=afterok:${dep} "$script_path")
  else
    jobid=$(sbatch --parsable --job-name="$name" --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=12:00:00 --output=logs/${name}_%j.log "$script_path")
  fi
  echo "Submitted $name as $jobid (script: $script_path)" >&2
  echo "$jobid"
}

submit_phase_gate() {
  local phase_name="$1" dep="$2"

  local body="bash \"${SCRIPT_PATH}\" __validate_phase__ \"${phase_name}\"
touch \"${WORKDIR}/.checkpoints/${phase_name}_gate.done\""

  local script_path
  script_path=$(write_job_script "${phase_name}_gate" "$body")

  local jobid
  jobid=$(sbatch --parsable --job-name="${phase_name}_gate" \
    --partition=compute --cpus-per-task=1 --mem-per-cpu=2000 --time=00:15:00 \
    --output="logs/${phase_name}_gate_%j.log" \
    --dependency=afterok:${dep} \
    "$script_path")
  echo "Submitted ${phase_name}_gate as $jobid (script: $script_path)" >&2
  echo "$jobid"
}

########################################################################
# ============================ PHASE 1 ================================
#   Neojunction Calling (steps 01-10)
########################################################################

if ! step_done "01_tumor_purity"; then
  jobPrev=$(submit_job "01_tumor_purity" "Rscript ${NEOJUNCTION_DIR}/01_tumor_purity.R" "" \
    "${OUTPUT_DIR:-SKIP}" 'Patient_List_Post_TumorPurity_Filter_[0-9]+\.[0-9]{2}\.txt')
else echo "01_tumor_purity already done."; fi

if ! step_done "02_protein_coding_genes"; then
  jobPrev=$(submit_job "02_protein_coding_genes" "Rscript ${NEOJUNCTION_DIR}/02_protein_coding_genes.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'GTF_Protein_Coding_Genes_[0-9]{8}\.txt')
else echo "02_protein_coding_genes already done."; fi

if ! step_done "03_tpm_filter_10"; then
  jobPrev=$(submit_job "03_tpm_filter_10" "Rscript ${NEOJUNCTION_DIR}/03_tpm_filter_10.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'GTF_ProteinCoding_Filter[0-9]+_[0-9]{8}\.tsv')
else echo "03_tpm_filter_10 already done."; fi

if ! step_done "04_extract_annot_sj_to_analyze"; then
  jobPrev=$(submit_job "04_extract_annot_sj_to_analyze" "Rscript ${NEOJUNCTION_DIR}/04_extract_annot_sj_to_analyze.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_[0-9]{8}\.tsv')
else echo "04_extract_annot_sj_to_analyze already done."; fi

if ! step_done "05a_extract_sj_nonannot_hartwig.count.min10"; then
  jobPrev=$(submit_job "05a_extract_sj_nonannot_hartwig.count.min10" "Rscript ${NEOJUNCTION_DIR}/05a_extract_sj_nonannot_hartwig.count.min10.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'SJ_List_NonAnnotated_Candidates_Protein_Coding_[0-9]{8}\.tsv')
else echo "05a_extract_sj_nonannot_hartwig.count.min10 already done."; fi

if ! step_done "05b_build_irfinder_reference"; then
  echo "[INFO] Step 05b: Building IRFinder reference..."
  BODY_05B="bash ${NEOJUNCTION_DIR}/05b_build_irfinder_ref.sh
bash \"${SCRIPT_PATH}\" __validate_dir__ 05b_build_irfinder_reference \"${WORKDIR}/0_Input_Files/irfinder_reference_grch38\" 5
touch ${WORKDIR}/.checkpoints/05b_build_irfinder_reference.done"
  SCRIPT_05B=$(write_job_script "05b_build_irfinder_reference" "$BODY_05B")
  JOB_05B=$(sbatch --parsable --job-name="05b_build_irfinder_reference" \
    --partition=compute --cpus-per-task=16 --mem=64G --time=4:00:00 \
    --output=logs/05b_build_irfinder_ref_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_05B")
  echo "Submitted 05b_build_irfinder_reference as $JOB_05B (script: $SCRIPT_05B)"
  jobPrev="$JOB_05B"
else
  echo "05b_build_irfinder_reference already done."
fi

if ! step_done "05c_run_irfinder_hartwig_samples"; then
  echo "[INFO] Step 05c: Submitting IRFinder single sequential job (all samples)..."
  BODY_05C="bash ${NEOJUNCTION_DIR}/05c_run_irfinder_hartwig_samples.sh
bash \"${SCRIPT_PATH}\" __validate_irfinder__ 05c_run_irfinder_hartwig_samples \"${WORKDIR}/results/irfinder_hartwig_results\" \"${WORKDIR}/samples.txt\" 0.9
touch ${WORKDIR}/.checkpoints/05c_run_irfinder_hartwig_samples.done"
  SCRIPT_05C=$(write_job_script "05c_run_irfinder_hartwig_samples" "$BODY_05C")
  JOB_05C=$(sbatch --parsable --job-name="05c_run_irfinder_hartwig_samples" \
    --partition=compute --cpus-per-task=8 --mem=64G --time=48:00:00 \
    --output=logs/05c_irfinder_hartwig_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_05C")
  echo "Submitted 05c_run_irfinder_hartwig_samples as $JOB_05C (script: $SCRIPT_05C)"
  jobPrev="$JOB_05C"
else
  echo "05c_run_irfinder_hartwig_samples already done."
fi

if ! step_done "05d_aggregate_irfinder_results"; then
  jobPrev=$(submit_job "05d_aggregate_irfinder_results" "Rscript ${NEOJUNCTION_DIR}/05d_aggregate_irfinder_results.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'IR_Candidates_PSR10_High_Confidence\.tsv')
else echo "05d_aggregate_irfinder_results already done."; fi

if ! step_done "05e_extract_sj_nonannot_hartwig.count.min10_protein.coding_psr10"; then
  jobPrev=$(submit_job "05e_extract_sj_nonannot_hartwig.count.min10_protein.coding_psr10" "Rscript ${NEOJUNCTION_DIR}/05e_extract_sj_nonannot_hartwig.count.min10_protein.coding_psr10.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'SJ_PSR_NonAnnotated_Candidates_Protein_Coding_[0-9]{8}\.tsv')
else echo "05e_extract_sj_nonannot_hartwig.count.min10_protein.coding_psr10 already done."; fi

if ! step_done "05f_merge_SJ_IR_junctions"; then
  jobPrev=$(submit_job "05f_merge_SJ_IR_junctions" "Rscript ${NEOJUNCTION_DIR}/05f_merge_SJ_IR_junctions.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'SJ_Novel_Potential_Neojunctions\.tsv')
else echo "05f_merge_SJ_IR_junctions already done."; fi

if ! step_done "06_prep_overlap_table"; then
  jobPrev=$(submit_job "06_prep_overlap_table" "Rscript ${NEOJUNCTION_DIR}/06_prep_overlap_table.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'SJ_Overlap_Table_[0-9]{8}\.tsv')
else echo "06_prep_overlap_table already done."; fi

if ! step_done "07_count_depth_freq_judge"; then
  echo "[INFO] Step 07: Submitting with optimized memory allocation (8 CPUs, 128GB)..."
  BODY_07="Rscript ${NEOJUNCTION_DIR}/07_count_depth_freq_judge.R
bash \"${SCRIPT_PATH}\" __validate_step__ 07_count_depth_freq_judge \"${OUTPUT_DIR:-SKIP}\" 'PSR_Table_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/07_count_depth_freq_judge.done"
  SCRIPT_07=$(write_job_script "07_count_depth_freq_judge" "$BODY_07")
  JOB_07=$(sbatch --parsable --job-name="07_count_depth_freq_judge" \
    --partition=compute --cpus-per-task=8 --mem=128G --time=12:00:00 \
    --output=logs/07_count_depth_freq_judge_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_07")
  echo "Submitted 07_count_depth_freq_judge as $JOB_07 (cpus=8, mem=128G, script: $SCRIPT_07)"
  jobPrev="$JOB_07"
else
  echo "07_count_depth_freq_judge already done."
fi

if ! step_done "08_overlap_table_gtex"; then
  jobPrev=$(submit_job "08_overlap_table_gtex" "Rscript ${NEOJUNCTION_DIR}/08_overlap_table_gtex.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'GTEx_SJ_Overlap_Table_[0-9]{8}\.tsv')
else echo "08_overlap_table_gtex already done."; fi

if ! step_done "09_count_depth_freq_judge_psr_gtex"; then
  jobPrev=$(submit_job "09_count_depth_freq_judge_psr_gtex" "Rscript ${NEOJUNCTION_DIR}/09_count_depth_freq_judge_psr_gtex.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'PSR_Retained_Passed_GTEx_[0-9]{8}\.tsv')
else echo "09_count_depth_freq_judge_psr_gtex already done."; fi

if ! step_done "10_extract_neojunctions"; then
  jobPrev=$(submit_job "10_extract_neojunctions" "Rscript ${NEOJUNCTION_DIR}/10_extract_neojunctions.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" 'PSR_Neojunctions_[0-9]{8}\.tsv')
else echo "10_extract_neojunctions already done."; fi

if ! step_done "phase1_gate"; then
  jobPrev=$(submit_phase_gate "phase1" "$jobPrev")
else
  echo "phase1_gate already done."
fi

########################################################################
# ============================ PHASE 2 ================================
#   Neopeptide Prediction + Presentation Prediction (steps 11-15d)
########################################################################

# ---------------------------------------------------------------------
# Step 11: AA sequence prediction
# ---------------------------------------------------------------------
if ! step_done "11_aaseq_prediction"; then
  jobPrev=$(submit_job "11_aaseq_prediction" \
    "Rscript ${NEOPEPTIDE_DIR}/11_aaseq_prediction.R" "$jobPrev" \
    "${STEP11_OUTPUT_DIR:-SKIP}" 'Res_AA_Prediction_Confirmed_[0-9]{8}\.tsv')
else echo "11_aaseq_prediction already done."; fi

# ---------------------------------------------------------------------
# Step 12: N-mer generation
# ---------------------------------------------------------------------
if ! step_done "12_nmer_generation"; then
  jobPrev=$(submit_job "12_nmer_generation" \
    "Rscript ${NEOPEPTIDE_DIR}/12_nmer_generation.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" '2023_0812_complete_list_all_mers\.tsv')
else echo "12_nmer_generation already done."; fi

# ---------------------------------------------------------------------
# Step 13a: NetMHCpan 4.2 predictions (ALT + WT, 8-11mers)
# Resources: compute, 24 CPUs, 8042MB/cpu (~192G), 3 days
# Outputs:   netmhcpan_08mer_YYYY_MMDD.tsv  (and 09/10/11mer equivalents)
#            netmhcpan_08mer_wt_YYYY_MMDD.tsv (and 09/10/11mer equivalents)
# ---------------------------------------------------------------------
if ! step_done "13a_NetMHCPan_analysis"; then
  echo "[INFO] Step 13a: NetMHCpan -- installs on first run, skips if already installed..."
  BODY_13A="bash ${PRESENTATION_DIR}/13a_NetMHCPan_analysis.sh
bash \"${SCRIPT_PATH}\" __validate_step__ 13a_NetMHCPan_analysis \"${OUTPUT_DIR:-SKIP}\" 'netmhcpan_08mer_[0-9_]+\.tsv' 2
touch ${WORKDIR}/.checkpoints/13a_NetMHCPan_analysis.done"
  SCRIPT_13A=$(write_job_script "13a_NetMHCPan_analysis" "$BODY_13A")
  JOB_13A=$(sbatch --parsable --job-name="13a_NetMHCPan_analysis" \
    --partition=compute --cpus-per-task=24 --mem-per-cpu=8042 --time=3-00:00:00 \
    --output=logs/13a_NetMHCPan_analysis_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_13A")
  echo "Submitted 13a_NetMHCPan_analysis as $JOB_13A (partition=compute, cpus=24, 3 days, script: $SCRIPT_13A)"
  jobPrev="$JOB_13A"
else
  echo "13a_NetMHCPan_analysis already done."
fi

# ---------------------------------------------------------------------
# Step 13b: Select top allele per peptide from NetMHCpan SB candidates
# Resources: compute, 8 CPUs, 8042MB/cpu (~64G), 2 hours
#   - Loads all 4 n-mer lengths (ALT + WT)
#   - Filters to binder == "<= SB", picks lowest EL rank per peptide
# Outputs:   netmhcpan_08mer_selected_alleles_YYYYMMDD.tsv  (+ 09/10/11)
#            netmhcpan_08mer_wt_selected_alleles_YYYYMMDD.tsv (+ 09/10/11)
# ---------------------------------------------------------------------
if ! step_done "13b_select_top_alleles"; then
  echo "[INFO] Step 13b: Selecting top NetMHCpan alleles per peptide (SB only)..."
  BODY_13B="Rscript ${PRESENTATION_DIR}/13b_select_top_alleles.R
bash \"${SCRIPT_PATH}\" __validate_step__ 13b_select_top_alleles \"${STEP13_OUTPUT_DIR:-SKIP}\" 'netmhcpan_09mer_selected_alleles_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/13b_select_top_alleles.done"
  SCRIPT_13B=$(write_job_script "13b_select_top_alleles" "$BODY_13B")
  JOB_13B=$(sbatch --parsable --job-name="13b_select_top_alleles" \
    --partition=compute --cpus-per-task=8 --mem-per-cpu=8042 --time=2:00:00 \
    --output=logs/13b_select_top_alleles_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_13B")
  echo "Submitted 13b_select_top_alleles as $JOB_13B (cpus=8, ~64G, 2h, script: $SCRIPT_13B)"
  jobPrev="$JOB_13B"
else
  echo "13b_select_top_alleles already done."
fi

# ---------------------------------------------------------------------
# Step 13c: Generate NetMHCpan prediction figures
# Resources: compute, 4 CPUs, 8042MB/cpu (~32G), 1 hour
#   - Histograms and pie charts from 13b selected-allele files
#   - All 4 n-mer lengths, full distribution + top 10 percentile
# Outputs:   PDFs in STEP13_FIGURES_DIR
# ---------------------------------------------------------------------
if ! step_done "13c_generate_figures"; then
  echo "[INFO] Step 13c: Generating NetMHCpan prediction figures..."
  BODY_13C="Rscript ${PRESENTATION_DIR}/13c_generate_figures.R
bash \"${SCRIPT_PATH}\" __validate_step__ 13c_generate_figures \"${OUTPUT_DIR}/figures/step13\" 'histogram_all_09mer_n[0-9]+_[0-9]{8}\.pdf' 2
touch ${WORKDIR}/.checkpoints/13c_generate_figures.done"
  SCRIPT_13C=$(write_job_script "13c_generate_figures" "$BODY_13C")
  JOB_13C=$(sbatch --parsable --job-name="13c_generate_figures" \
    --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=1:00:00 \
    --output=logs/13c_generate_figures_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_13C")
  echo "Submitted 13c_generate_figures as $JOB_13C (cpus=4, ~32G, 1h, script: $SCRIPT_13C)"
  jobPrev="$JOB_13C"
else
  echo "13c_generate_figures already done."
fi

# ---------------------------------------------------------------------
# Step 14a: MHCflurry input dataframe generation
# Resources: compute, 24 CPUs, 8042MB/cpu (~192G), 12 hours
#   - WT tables are ~196 million rows each -- high memory required
# Outputs:   08mer_mhcflurry_input_wt_YYYY_MMDD.csv (+ 09/10/11)
# ---------------------------------------------------------------------
if ! step_done "14a_mhcflurry2_input_df_generation"; then
  echo "[INFO] Step 14a: Submitting with higher memory (WT tables are ~196 million rows each)..."
  BODY_14A="Rscript ${PRESENTATION_DIR}/14a_mhcflurry2_input_df_generation.R
bash \"${SCRIPT_PATH}\" __validate_step__ 14a_mhcflurry2_input_df_generation \"${STEP14_OUTPUT_DIR:-SKIP}\" '08mer_mhcflurry_input_wt_[0-9_]+\.csv' 2
touch ${WORKDIR}/.checkpoints/14a_mhcflurry2_input_df_generation.done"
  SCRIPT_14A=$(write_job_script "14a_mhcflurry2_input_df_generation" "$BODY_14A")
  JOB_14A=$(sbatch --parsable --job-name="14a_mhcflurry2_input_df_generation" \
    --partition=compute --cpus-per-task=24 --mem-per-cpu=8042 --time=12:00:00 \
    --output=logs/14a_mhcflurry2_input_df_generation_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_14A")
  echo "Submitted 14a_mhcflurry2_input_df_generation as $JOB_14A (cpus=24, ~192G, 12h, script: $SCRIPT_14A)"
  jobPrev="$JOB_14A"
else
  echo "14a_mhcflurry2_input_df_generation already done."
fi

# ---------------------------------------------------------------------
# Step 14b: MHCflurry 2.0 predictions with flanks (GPU)
# Resources: gpu, 1 GPU, 6 CPUs, 90G, 7 days
#   - WT files are 17-18GB each -- high memory + GPU required
# Outputs:   YYYYMMDD_08mers_flank_mhcflurry.csv (+ 09/10/11, ALT+WT)
# ---------------------------------------------------------------------
if ! step_done "14b_mhcflurry2_analysis_with_flank"; then
  echo "[INFO] Step 14b: Submitting to GPU partition with high memory (WT files are 17-18GB each)..."
  BODY_14B="bash ${PRESENTATION_DIR}/14b_mhcflurry2_analysis_with_flank.sh
bash \"${SCRIPT_PATH}\" __validate_step__ 14b_mhcflurry2_analysis_with_flank \"${STEP14_OUTPUT_DIR:-SKIP}\" '08mers_flank_mhcflurry_wt_[0-9_]+\.csv' 2
touch ${WORKDIR}/.checkpoints/14b_mhcflurry2_analysis_with_flank.done"
  SCRIPT_14B=$(write_job_script "14b_mhcflurry2_analysis_with_flank" "$BODY_14B")
  JOB_14B=$(sbatch --parsable --job-name="14b_mhcflurry2_analysis_with_flank" \
    --partition=gpu --gpus-per-node=1 --cpus-per-task=6 --mem=90G --time=7-00:00:00 \
    --output=logs/14b_mhcflurry2_analysis_with_flank_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_14B")
  echo "Submitted 14b_mhcflurry2_analysis_with_flank as $JOB_14B (partition=gpu, 6 cpus/90G, single-GPU, 7 days, script: $SCRIPT_14B)"
  jobPrev="$JOB_14B"
else
  echo "14b_mhcflurry2_analysis_with_flank already done."
fi

# ---------------------------------------------------------------------
# Step 14c: Select top MHCflurry alleles
# Resources: compute, 24 CPUs, 8042MB/cpu (~192G), 2 days
#   - Also loads WT files (17-18GB each), processes one at a time
# Outputs:   mhcflurry_08mer_top_YYYYMMDD.tsv (+ 09/10/11)
# ---------------------------------------------------------------------
if ! step_done "14c_select_top_alleles"; then
  echo "[INFO] Step 14c: Submitting with high memory (now also loads WT files, 17-18GB each)..."
  BODY_14C="Rscript ${PRESENTATION_DIR}/14c_select_top_alleles.R
bash \"${SCRIPT_PATH}\" __validate_step__ 14c_select_top_alleles \"${STEP14_OUTPUT_DIR:-SKIP}\" 'mhcflurry_08mer_top_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/14c_select_top_alleles.done"
  SCRIPT_14C=$(write_job_script "14c_select_top_alleles" "$BODY_14C")
  JOB_14C=$(sbatch --parsable --job-name="14c_select_top_alleles" \
    --partition=compute --cpus-per-task=24 --mem-per-cpu=8042 --time=2-00:00:00 \
    --output=logs/14c_select_top_alleles_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_14C")
  echo "Submitted 14c_select_top_alleles as $JOB_14C (cpus=24, ~192G, 2 days, script: $SCRIPT_14C)"
  jobPrev="$JOB_14C"
else
  echo "14c_select_top_alleles already done."
fi

# ---------------------------------------------------------------------
# Step 14d: Generate MHCflurry figures
# Resources: compute, 4 CPUs, 8042MB/cpu (~32G), 1 hour
# ---------------------------------------------------------------------
if ! step_done "14d_generate_figures"; then
  jobPrev=$(submit_job "14d_generate_figures" \
    "Rscript ${PRESENTATION_DIR}/14d_generate_figures.R" "$jobPrev")
else echo "14d_generate_figures already done."; fi

# =====================================================================
# Steps 15a-15d: Cross-analysis (NetMHCpan x MHCflurry) -- Concordance
# NOTE: 13a + 14b must BOTH be complete before 15a can run.
#   13a raw TSVs  -> 15a (NetMHCpan all binding levels)
#   14b raw CSVs  -> 15a (MHCflurry full prediction universe)
# =====================================================================

# ---------------------------------------------------------------------
# Step 15a: Concordance cross-analysis (NetMHCpan + MHCflurry)
# Resources: compute, 8 CPUs, 8042MB/cpu (~64G), 4 hours
#   TWO-STEP LOGIC:
#   Step 1: Build WT exclusion set from ALL WT predictions (both tools,
#           any binding level) -- these peptides are removed from ALT.
#   Step 2: Assign concordance tiers to clean ALT set AND WT separately:
#     Tier 1 (High confidence):   NMP EL rank <0.5% AND MHCflurry <500nM
#     Tier 2 (Medium confidence): NMP EL rank <2.0% AND MHCflurry <500nM
#     Tier 3 (Discordant):        Tools disagree -- flagged, kept for ref
# Outputs:   alt_concordance_tier1_YYYYMMDD.tsv  -- both tools strong
#            alt_concordance_tier2_YYYYMMDD.tsv  -- moderate agreement
#            alt_concordance_tier3_YYYYMMDD.tsv  -- discordant
#            alt_concordance_all_YYYYMMDD.tsv    -- all tiers combined
#            wt_concordance_tier1_YYYYMMDD.tsv   -- WT native Tier 1
#            wt_concordance_all_YYYYMMDD.tsv     -- WT all tiers
#            wt_exclusion_peptides_YYYYMMDD.txt  -- excluded WT peptides
#            cross_alg_all_nmers_YYYYMMDD.tsv    -- full join (15b/c/d)
# ---------------------------------------------------------------------
if ! step_done "15a_generate_combined_crossanalysis_dataframes"; then
  echo "[INFO] Step 15a: Building concordance cross-analysis (NMP x MHCflurry)..."
  BODY_15A="Rscript ${PRESENTATION_DIR}/15a_generate_combined_crossanalysis_dataframes.R
bash \"${SCRIPT_PATH}\" __validate_step__ 15a_generate_combined_crossanalysis_dataframes \"${STEP15_OUTPUT_DIR:-SKIP}\" 'alt_concordance_all_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/15a_generate_combined_crossanalysis_dataframes.done"
  SCRIPT_15A=$(write_job_script "15a_generate_combined_crossanalysis_dataframes" "$BODY_15A")
  JOB_15A=$(sbatch --parsable --job-name="15a_crossanalysis_df" \
    --partition=compute --cpus-per-task=8 --mem-per-cpu=8042 --time=4:00:00 \
    --output=logs/15a_generate_combined_crossanalysis_dataframes_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_15A")
  echo "Submitted 15a_generate_combined_crossanalysis_dataframes as $JOB_15A (cpus=8, ~64G, 4h, script: $SCRIPT_15A)"
  jobPrev="$JOB_15A"
else
  echo "15a_generate_combined_crossanalysis_dataframes already done."
fi

# ---------------------------------------------------------------------
# Step 15b: Plot concordance tier figures (ALT + WT)
# Resources: compute, 4 CPUs, 8042MB/cpu (~32G), 1 hour
#   - Scatter: NMP EL score vs MHCflurry presentation score, by tier
#   - Scatter: NMP EL rank vs MHCflurry affinity (nM) with thresholds
#   - Bar chart: peptide counts per concordance tier
#   (Figures generated for both ALT and WT)
# Outputs:   alt_scatter_ELscore_vs_pres_YYYYMMDD.pdf
#            alt_scatter_ELrank_vs_affinity_YYYYMMDD.pdf
#            alt_barplot_tier_counts_YYYYMMDD.pdf
#            wt_scatter_*.pdf  wt_barplot_*.pdf
#            cross_analysis_summary_nmp_mf_YYYYMMDD.tsv
# ---------------------------------------------------------------------
if ! step_done "15b_plot_cross_analysis"; then
  echo "[INFO] Step 15b: Generating concordance tier figures..."
  BODY_15B="Rscript ${PRESENTATION_DIR}/15b_plot_cross_analysis.R
bash \"${SCRIPT_PATH}\" __validate_step__ 15b_plot_cross_analysis \"${OUTPUT_DIR}/figures/step15\" 'alt_scatter_ELscore_vs_pres_[0-9]{8}\.pdf' 2
touch ${WORKDIR}/.checkpoints/15b_plot_cross_analysis.done"
  SCRIPT_15B=$(write_job_script "15b_plot_cross_analysis" "$BODY_15B")
  JOB_15B=$(sbatch --parsable --job-name="15b_plot_cross_analysis" \
    --partition=compute --cpus-per-task=4 --mem-per-cpu=8042 --time=1:00:00 \
    --output=logs/15b_plot_cross_analysis_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_15B")
  echo "Submitted 15b_plot_cross_analysis as $JOB_15B (cpus=4, ~32G, 1h, script: $SCRIPT_15B)"
  jobPrev="$JOB_15B"
else
  echo "15b_plot_cross_analysis already done."
fi

# ---------------------------------------------------------------------
# Step 15c: Map n-mers back to neojunctions (ALT all tiers + WT Tier 1)
# Resources: compute, 8 CPUs, 8042MB/cpu (~64G), 2 hours
#   - Loads 2023_0812_complete_list_all_mers.tsv (junction reference)
#   - For each unique peptide: str_detect(aa.seq.alt, peptide) to find
#     which junction(s) it originated from -- gets junc.id, type, fs
#   - Runs for ALL ALT tiers (preserves tier label in output)
#   - Also maps WT Tier 1 for native immunopeptidome reference
# Outputs:   alt_neoA_to_neoJ_map_YYYYMMDD.tsv  -- peptide x junction
#            alt_immunogenic_njs_YYYYMMDD.tsv    -- per-junction summary
#            wt_neoA_to_neoJ_map_YYYYMMDD.tsv
#            wt_immunogenic_njs_YYYYMMDD.tsv
# ---------------------------------------------------------------------
if ! step_done "15c_map_nmers_to_nj"; then
  echo "[INFO] Step 15c: Mapping n-mers to neojunctions via aa.seq.alt substring search..."
  BODY_15C="Rscript ${PRESENTATION_DIR}/15c_map_nmers_to_nj.R
bash \"${SCRIPT_PATH}\" __validate_step__ 15c_map_nmers_to_nj \"${STEP15_OUTPUT_DIR:-SKIP}\" 'alt_neoA_to_neoJ_map_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/15c_map_nmers_to_nj.done"
  SCRIPT_15C=$(write_job_script "15c_map_nmers_to_nj" "$BODY_15C")
  JOB_15C=$(sbatch --parsable --job-name="15c_map_nmers_to_nj" \
    --partition=compute --cpus-per-task=8 --mem-per-cpu=8042 --time=2:00:00 \
    --output=logs/15c_map_nmers_to_nj_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_15C")
  echo "Submitted 15c_map_nmers_to_nj as $JOB_15C (cpus=8, ~64G, 2h, script: $SCRIPT_15C)"
  jobPrev="$JOB_15C"
else
  echo "15c_map_nmers_to_nj already done."
fi

# ---------------------------------------------------------------------
# Step 15d: Matrix + bed files + figures (Tier 1 only)
# Resources: compute, 8 CPUs, 8042MB/cpu (~64G), 2 hours
#   - Uses ONLY Tier 1 (both tools strong) from 15c mapped output
#   - Bed files: ENST_ID|AA_START|AA_END|PEPTIDE|HLA_ALLELES
#     -> immunogenic_peptides_YYYYMMDD.bed  (tumour-specific)
#     -> wt_native_immunogenic_peptides_YYYYMMDD.bed   (native reference)
#   - Matrix: neo_id (peptide|junc.id|allele) x sample (307 samples)
#   - HLA summary, sample summary, figures by NJ type + frameshift
# Outputs:   immunogenic_peptides_YYYYMMDD.bed
#            wt_native_immunogenic_peptides_YYYYMMDD.bed
#            top_neoantigens_matrix_YYYYMMDD.tsv
#            top_neoantigens_long_YYYYMMDD.tsv
#            hla_alleles_summary_YYYYMMDD.tsv
#            sample_neoantigen_summary_YYYYMMDD.tsv
#            PDFs in OUTPUT_DIR/figures/step15
# ---------------------------------------------------------------------
if ! step_done "15d_scores_for_nj_types"; then
  echo "[INFO] Step 15d: Building Tier 1 matrix + bed files + figures..."
  BODY_15D="Rscript ${PRESENTATION_DIR}/15d_scores_for_nj_types.R
bash \"${SCRIPT_PATH}\" __validate_step__ 15d_scores_for_nj_types \"${STEP15_OUTPUT_DIR:-SKIP}\" 'top_neoantigens_matrix_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/15d_scores_for_nj_types.done"
  SCRIPT_15D=$(write_job_script "15d_scores_for_nj_types" "$BODY_15D")
  JOB_15D=$(sbatch --parsable --job-name="15d_scores_for_nj_types" \
    --partition=compute --cpus-per-task=8 --mem-per-cpu=8042 --time=2:00:00 \
    --output=logs/15d_scores_for_nj_types_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_15D")
  echo "Submitted 15d_scores_for_nj_types as $JOB_15D (cpus=8, ~64G, 2h, script: $SCRIPT_15D)"
  jobPrev="$JOB_15D"
else
  echo "15d_scores_for_nj_types already done."
fi

# ---------------------------------------------------------------------
# Phase 2 gate
# ---------------------------------------------------------------------
if ! step_done "phase2_gate"; then
  jobPrev=$(submit_phase_gate "phase2" "$jobPrev")
else
  echo "phase2_gate already done."
fi

echo "=== Submission complete. phase1 (01-10) -> phase1_gate -> phase2 (11-15d) -> phase2_gate. ==="
