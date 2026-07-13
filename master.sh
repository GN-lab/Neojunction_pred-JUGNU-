#!/bin/bash
#SBATCH --job-name=ssnip_driver
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
#  file. Normal `sbatch run_ssnip_pipeline_master.sh` submissions never hit
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

echo "=== Submitting SSNIP Pipeline (dependency-driven, 2 phases, self-validating) ==="

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
#   Neopeptide Prediction + Presentation Prediction (steps 11-15e)
########################################################################

if ! step_done "11_aaseq_prediction"; then
  jobPrev=$(submit_job "11_aaseq_prediction" "Rscript ${NEOPEPTIDE_DIR}/11_aaseq_prediction.R" "$jobPrev" \
    "${STEP11_OUTPUT_DIR:-SKIP}" 'Res_AA_Prediction_Confirmed_[0-9]{8}\.tsv')
else echo "11_aaseq_prediction already done."; fi

if ! step_done "12_nmer_generation"; then
  jobPrev=$(submit_job "12_nmer_generation" "Rscript ${NEOPEPTIDE_DIR}/12_nmer_generation.R" "$jobPrev" \
    "${OUTPUT_DIR:-SKIP}" '2023_0812_complete_list_all_mers\.tsv')
else echo "12_nmer_generation already done."; fi

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
  echo "Submitted 14a_mhcflurry2_input_df_generation as $JOB_14A (cpus=24, ~64G total, script: $SCRIPT_14A)"
  jobPrev="$JOB_14A"
else
  echo "14a_mhcflurry2_input_df_generation already done."
fi

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
  echo "Submitted 14b_mhcflurry2_analysis_with_flank as $JOB_14B (partition=gpu, 6 cpus/90G, single-GPU share -- script: $SCRIPT_14B)"
  jobPrev="$JOB_14B"
else
  echo "14b_mhcflurry2_analysis_with_flank already done."
fi

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
  echo "Submitted 14c_select_top_alleles as $JOB_14C (cpus=24, ~8042MB/cpu = ~64G total -- 14c now processes WT files one at a time, so peak memory is much lower than before, script: $SCRIPT_14C)"
  jobPrev="$JOB_14C"
else
  echo "14c_select_top_alleles already done."
fi

if ! step_done "14d_generate_figures"; then
  jobPrev=$(submit_job "14d_generate_figures" "Rscript ${PRESENTATION_DIR}/14d_generate_figures.R" "$jobPrev")
else echo "14d_generate_figures already done."; fi

if ! step_done "15d_scores_for_nj_types"; then
  echo "[INFO] Step 15d: Submitting with 5-day time limit..."
  BODY_15D="Rscript ${PRESENTATION_DIR}/15d_scores_for_nj_types.R
bash \"${SCRIPT_PATH}\" __validate_step__ 15d_scores_for_nj_types \"${STEP15_OUTPUT_DIR:-SKIP}\" 'top_neoantigens_matrix_[0-9]{8}\.tsv' 2
touch ${WORKDIR}/.checkpoints/15d_scores_for_nj_types.done"
  SCRIPT_15D=$(write_job_script "15d_scores_for_nj_types" "$BODY_15D")
  JOB_15D=$(sbatch --parsable --job-name="15d_scores_for_nj_types" \
    --partition=compute --cpus-per-task=24 --mem-per-cpu=8042 --time=5-00:00:00 \
    --output=logs/15d_scores_for_nj_types_%j.log \
    ${jobPrev:+--dependency=afterok:${jobPrev}} \
    "$SCRIPT_15D")
  echo "Submitted 15d_scores_for_nj_types as $JOB_15D (time=5-00:00:00 = 5 days, script: $SCRIPT_15D)"
  jobPrev="$JOB_15D"
else
  echo "15d_scores_for_nj_types already done."
fi

if ! step_done "phase2_gate"; then
  jobPrev=$(submit_phase_gate "phase2" "$jobPrev")
else
  echo "phase2_gate already done."
fi

echo "=== Submission complete. phase1 (01-10) -> phase1_gate -> phase2 (11-15e) -> phase2_gate. ==="
