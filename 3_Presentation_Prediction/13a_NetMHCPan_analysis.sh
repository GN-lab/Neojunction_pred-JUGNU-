#!/bin/bash
#SBATCH --job-name=netmhcpan
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8042
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/13a_netmhcpan_%j.log
#SBATCH --error=logs/13a_netmhcpan_%j.err

set -euo pipefail

# Title: Step 13a: NetMHCpan 4.1 Predictions
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: NetMHCpan 4.1 predicts MHC class I binding affinity and eluted
# ligand (EL) scores for cancer-specific n-mers generated in Step 12.
# Runs in BA+EL mode for 8-11mers across all patient-specific HLA alleles.
# Replaces HLAthena in the JUGNU pipeline.
#
# All paths come from config.sh -- no hardcoded values in this script.
# Required config.sh variables (add these if not already present):
#   NETMHCPAN_INSTALL_DIR  -- where to install/find NetMHCpan binary
#   NETMHCPAN_TAR          -- path to netMHCpan-4.1b.Linux.tar
#
# Install behaviour:
#   First run  : extracts tar, configures NMHOME, tests, then predicts
#   Subsequent : binary found -> skips install, runs predictions directly

###########################################################################
#  Load config ------------------------------------------------------------
###########################################################################

source "${WORKDIR}/config.sh"

# Validate required config variables
: "${OUTPUT_DIR:?OUTPUT_DIR not set -- source config.sh first}"
: "${OPTITYPE_OUTPUT_DIR:?OPTITYPE_OUTPUT_DIR not set in config.sh}"
: "${NETMHCPAN_BIN:?NETMHCPAN_BIN not set in config.sh}"

RUN_DATE=$(date +%Y_%m%d)

echo "[CONFIG] WORKDIR:             ${WORKDIR}"
echo "[CONFIG] OUTPUT_DIR:          ${OUTPUT_DIR}"
echo "[CONFIG] OPTITYPE_OUTPUT_DIR: ${OPTITYPE_OUTPUT_DIR}"
echo "[CONFIG] NETMHCPAN_BIN:       ${NETMHCPAN_BIN}"

###########################################################################
#  Step 0: Validate NetMHCpan binary --------------------------------------
###########################################################################
# NetMHCpan must be installed manually by each user before running JUGNU.
# Installation steps:
#   1. Download netMHCpan-4.2estatic.Linux.tar from:
#      https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/
#   2. Extract: tar -xf netMHCpan-4.2estatic.Linux.tar
#   3. Configure NMHOME inside the netMHCpan script (use python3 replace,
#      not sed, due to tab characters):
#      python3 -c "
#        c=open('netMHCpan-4.2/netMHCpan').read()
#        c=c.replace('setenv\tNMHOME\t/tools/src/netMHCpan-4.2',
#                    'setenv NMHOME /your/install/path/netMHCpan-4.2')
#        open('netMHCpan-4.2/netMHCpan','w').write(c)"
#   4. Set NETMHCPAN_BIN in config.sh

if [[ ! -x "${NETMHCPAN_BIN}" ]]; then
    echo "[ERROR] NetMHCpan binary not found or not executable: ${NETMHCPAN_BIN}"
    echo "  Please install NetMHCpan 4.2 manually -- see instructions above."
    exit 1
fi

# Sanity check with a known strong binder
echo "GILGFVFTL" > "${TMPDIR}/test_netmhcpan.pep"
if ! "${NETMHCPAN_BIN}" -a HLA-A02:01 -p "${TMPDIR}/test_netmhcpan.pep" -l 9 \
    > "${TMPDIR}/test_netmhcpan.out" 2>&1; then
    echo "[ERROR] NetMHCpan sanity check failed. Output:"
    cat "${TMPDIR}/test_netmhcpan.out"
    exit 1
fi

echo "[INFO] NetMHCpan validated successfully: ${NETMHCPAN_BIN}"

###########################################################################
#  Step 1: Build HLA allele list from OptiType output --------------------
###########################################################################
# Reads all OptiType *_result.tsv files from OPTITYPE_OUTPUT_DIR
# Converts A*02:01 format -> HLA-A02:01 (NetMHCpan format)

ALLELE_FILE="${TMPDIR}/netmhcpan_alleles_${RUN_DATE}.txt"

python3 - <<PYEOF
import os, glob, re

optitype_dir = "${OPTITYPE_OUTPUT_DIR}"
outfile      = "${ALLELE_FILE}"

alleles = set()
for f in glob.glob(os.path.join(optitype_dir, "**", "*_result.tsv"), recursive=True):
    with open(f) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            for col in parts[1:7]:   # A1, A2, B1, B2, C1, C2
                col = col.strip()
                if re.match(r"[ABC]\*\d+:\d+", col):
                    locus  = col[0]
                    fields = col[2:]  # strip "A*" prefix
                    alleles.add(f"HLA-{locus}{fields}")

with open(outfile, "w") as out:
    for a in sorted(alleles):
        out.write(a + "\n")

print(f"[INFO] Wrote {len(alleles)} unique HLA alleles to {outfile}")
PYEOF

ALLELE_LIST=$(paste -sd "," "${ALLELE_FILE}")
N_ALLELES=$(wc -l < "${ALLELE_FILE}")
echo "[INFO] ${N_ALLELES} unique HLA alleles loaded"

###########################################################################
#  Step 2: Run NetMHCpan predictions (8-11mers) --------------------------
###########################################################################
# Input:  hlathenalist_msic_0Xmers.tsv from Step 12 (n_mer column)
# Output: netmhcpan_0Xmer_YYYY_MMDD.txt (raw)
#         netmhcpan_0Xmer_YYYY_MMDD.tsv (parsed, for downstream R scripts)

run_netmhcpan() {
    local LENGTH=$1
    local TYPE=$2
    local LENPAD
    LENPAD=$(printf "%02d" "${LENGTH}")

    if [[ "${TYPE}" == "wt" ]]; then
        local INPUT_TSV="${OUTPUT_DIR}/2023_0812_hlathenalist_msic_${LENPAD}mers_wt.tsv"
        local INPUT_PEP="${TMPDIR}/netmhcpan_${LENPAD}mer_wt_peptides.pep"
        local OUT_RAW="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_wt_${RUN_DATE}.txt"
        local OUT_TSV="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_wt_${RUN_DATE}.tsv"
    else
        local INPUT_TSV="${OUTPUT_DIR}/2023_0812_hlathenalist_msic_${LENPAD}mers.tsv"
        local INPUT_PEP="${TMPDIR}/netmhcpan_${LENPAD}mer_peptides.pep"
        local OUT_RAW="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_${RUN_DATE}.txt"
        local OUT_TSV="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_${RUN_DATE}.tsv"
    fi

    if [[ ! -f "${INPUT_TSV}" ]] || [[ ! -s "${INPUT_TSV}" ]]; then
        echo "[WARN] Input not found: ${INPUT_TSV} -- skipping ${TYPE} ${LENPAD}mers"
        return
    fi

    # ------------------------------------------------------------------
    # Skip NetMHCpan run if raw .txt already exists (e.g. previous run).
    # Always regenerate the .tsv from the .txt so the parser fix applies.
    # ------------------------------------------------------------------
    if [[ -f "${OUT_RAW}" ]] && [[ -s "${OUT_RAW}" ]]; then
        echo "[INFO] ${TYPE} ${LENPAD}mer: raw .txt already exists -- skipping NetMHCpan run"
        echo "[INFO]   -> ${OUT_RAW}"
    else
        awk -F'\t' 'NR>1 && $1!="" {print $1}' "${INPUT_TSV}" | sort -u > "${INPUT_PEP}"
        local N_PEPS
        N_PEPS=$(wc -l < "${INPUT_PEP}")
        echo "[INFO] ${TYPE} ${LENPAD}mers: ${N_PEPS} unique peptides"

        if [[ ${N_PEPS} -eq 0 ]]; then
            echo "[WARN] No peptides for ${TYPE} ${LENPAD}mers -- skipping"
            return
        fi

        # NetMHCpan -a flag has a 1024 char limit -- split alleles into batches of 20
        echo "[INFO] Running NetMHCpan ${TYPE} ${LENPAD}mers (batched alleles)..."
        > "${OUT_RAW}"  # create/clear output file

        local BATCH_SIZE=20
        local ALLELES_ARRAY
        IFS=',' read -ra ALLELES_ARRAY <<< "${ALLELE_LIST}"
        local TOTAL=${#ALLELES_ARRAY[@]}
        local BATCH_NUM=0

        for (( i=0; i<TOTAL; i+=BATCH_SIZE )); do
            BATCH_NUM=$((BATCH_NUM + 1))
            local BATCH_ALLELES
            BATCH_ALLELES=$(IFS=','; echo "${ALLELES_ARRAY[*]:$i:$BATCH_SIZE}")
            echo "[INFO]   Batch ${BATCH_NUM}: alleles $((i+1))-$((i+BATCH_SIZE < TOTAL ? i+BATCH_SIZE : TOTAL)) of ${TOTAL}"

            "${NETMHCPAN_BIN}" \
                -a "${BATCH_ALLELES}" \
                -p "${INPUT_PEP}" \
                -l "${LENGTH}" \
                -BA \
                >> "${OUT_RAW}"
        done

        echo "[INFO] ${TYPE} ${LENPAD}mer raw output: ${OUT_RAW}"
    fi

    # Parse raw output to TSV
    python3 - <<PYEOF
infile  = "${OUT_RAW}"
outfile = "${OUT_TSV}"
rows    = []

with open(infile) as fh:
    for line in fh:
        line = line.rstrip()
        if not line or line.startswith("#") or line.startswith("-") \
           or line.startswith(" Pos") or line.startswith("Protein") \
           or line.startswith("Error"):
            continue
        parts = line.split()
        if len(parts) < 11:
            continue
        try:
            allele   = parts[1]
            peptide  = parts[2]
            core     = parts[3]
            # Column positions (0-indexed):
            # 0=Pos 1=MHC 2=Peptide 3=Core 4=Of 5=Gp 6=Gl 7=Ip 8=Il 9=Icore
            # 10=Identity 11=Score_EL 12=%Rank_EL 13=Score_BA 14=%Rank_BA 15=BindLevel
            el_score = parts[11]
            el_rank  = parts[12]
            ba_score = parts[13] if len(parts) > 13 else "NA"
            ba_rank  = parts[14] if len(parts) > 14 else "NA"
            # NetMHCpan splits bind level across two tokens: "<=" and "SB"/"WB"
            # so parts[-1] == "SB" or "WB", and parts[-2] == "<=" for binders.
            # Non-binders have no bind-level tokens, so parts[-1] is a number.
            if len(parts) >= 16 and parts[-2] == "<=" and parts[-1] == "SB":
                binder = "SB"
            elif len(parts) >= 16 and parts[-2] == "<=" and parts[-1] == "WB":
                binder = "WB"
            else:
                binder = "NB"
            if not allele.startswith("HLA") or not peptide.isalpha():
                continue
            rows.append([allele, peptide, core,
                         el_score, el_rank, ba_score, ba_rank, binder])
        except (IndexError, ValueError):
            continue

with open(outfile, "w") as out:
    out.write("allele\tpeptide\tcore\t"
              "netmhcpan_EL_score\tnetmhcpan_EL_rank\t"
              "netmhcpan_BA_score\tnetmhcpan_BA_rank\tbinder\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

print(f"[INFO] Parsed {len(rows)} rows -> {outfile}")
PYEOF

    echo "[INFO] ${TYPE} ${LENPAD}mer TSV: ${OUT_TSV}"
}

# Run ALT and WT for all four lengths in parallel
echo "[INFO] Running ALT predictions..."
for LEN in 8 9 10 11; do
    run_netmhcpan ${LEN} "alt" &
done
wait
echo "[INFO] ALT predictions complete."

echo "[INFO] Running WT predictions..."
for LEN in 8 9 10 11; do
    run_netmhcpan ${LEN} "wt" &
done
wait
echo "[INFO] WT predictions complete."

###########################################################################
#  Summary ----------------------------------------------------------------
###########################################################################

echo ""
echo "=== NetMHCpan predictions complete ==="
echo "ALT (cancer-specific):"
for LEN in 8 9 10 11; do
    LENPAD=$(printf "%02d" "${LEN}")
    TSV="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_${RUN_DATE}.tsv"
    if [[ -f "${TSV}" ]]; then
        ROWS=$(wc -l < "${TSV}")
        echo "  ${LENPAD}mer: ${ROWS} rows -> ${TSV}"
    else
        echo "  ${LENPAD}mer: MISSING"
    fi
done
echo "WT (native self-peptides):"
for LEN in 8 9 10 11; do
    LENPAD=$(printf "%02d" "${LEN}")
    TSV="${OUTPUT_DIR}/netmhcpan_${LENPAD}mer_wt_${RUN_DATE}.tsv"
    if [[ -f "${TSV}" ]]; then
        ROWS=$(wc -l < "${TSV}")
        echo "  ${LENPAD}mer: ${ROWS} rows -> ${TSV}"
    else
        echo "  ${LENPAD}mer: MISSING"
    fi
done