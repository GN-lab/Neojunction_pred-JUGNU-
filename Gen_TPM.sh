#!/bin/bash
conda activate python_env

# ========= SETUP =========
WORKDIR="/data"
SALMON_DIR="/data/nextflow_rnaseq/output/star_salmon"
OUTPUT_DIR="${WORKDIR}/0_Input_Files"
mkdir -p "${OUTPUT_DIR}"

echo "=== CLEAN TPM MATRIX GENERATION (SINGLE DIR: SF3b1) ==="
echo ""

# Input TPM file (only one)
TPM_FILE="${SALMON_DIR}/salmon.merged.transcript_tpm.tsv"
OUTPUT_TPM="${OUTPUT_DIR}/transcript_tpm_matrix.tsv"

echo "1. Checking input file..."
ls -lh "$TPM_FILE" || { echo "Error: TPM file not found!"; exit 1; }

echo ""
echo "First 2 lines of TPM file:"
head -2 "$TPM_FILE" | cut -f1-5
echo "Number of columns: $(head -1 "$TPM_FILE" | tr '\t' '\n' | wc -l)"

# ========= CREATE TRANSCRIPT TPM MATRIX =========
cat > "${WORKDIR}/scripts/create_tpm.py" << 'EOF'
#!/usr/bin/env python3
import sys

tpm_file = "/data/nextflow_rnaseq/output/star_salmon/salmon.merged.transcript_tpm.tsv"
output_file = "0_Input_Files/transcript_tpm_matrix.tsv"

print("Reading TPM file...")
with open(tpm_file, 'r') as f:
    header = f.readline().strip().split('\t')
    samples = header[2:]  # skip tx, gene_id
    data = {}
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        tx, gene_id = parts[0], parts[1]
        tpm_values = parts[2:]
        data[tx] = {'gene_id': gene_id, 'tpm': dict(zip(samples, tpm_values))}

print(f"Found {len(data)} transcripts across {len(samples)} samples.")
print("Writing transcript-level TPM matrix...")

with open(output_file, 'w') as out:
    out.write("tx\tgene_id\t" + "\t".join(samples) + "\n")
    for tx in sorted(data.keys()):
        gene_id = data[tx]['gene_id']
        tpm_values = [data[tx]['tpm'].get(s, '0') for s in samples]
        out.write(f"{tx}\t{gene_id}\t" + "\t".join(tpm_values) + "\n")

print(f"Done! Output written to: {output_file}")
EOF

chmod +x "${WORKDIR}/scripts/create_tpm.py"
python3 "${WORKDIR}/scripts/create_tpm.py"

echo ""
echo "Verifying transcript output..."
echo "File size: $(ls -lh "$OUTPUT_TPM" | awk '{print $5}')"
head -3 "$OUTPUT_TPM" | cut -f1-5
echo "Samples: $(head -1 "$OUTPUT_TPM" | tr '\t' '\n' | grep -v -E '^(tx|gene_id)$' | wc -l)"
echo "Transcripts: $(($(wc -l < "$OUTPUT_TPM") - 1))"
echo ""

# ========= GENE-LEVEL TPM =========
cat > "${WORKDIR}/scripts/create_gene_tpm.py" << 'EOF'
#!/usr/bin/env python3
from collections import defaultdict

transcript_file = "0_Input_Files/transcript_tpm_matrix.tsv"
gene_matrix_out = "0_Input_Files/gene_tpm_matrix.tsv"
gene_avg_out = "0_Input_Files/gene_tpm.tsv"

print("Reading transcript TPM matrix...")
gene_tpm = defaultdict(lambda: defaultdict(float))

with open(transcript_file, 'r') as f:
    header = f.readline().strip().split('\t')
    samples = header[2:]
    for i, line in enumerate(f, 1):
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        gene_id = parts[1]
        tpm_values = parts[2:]
        for j, s in enumerate(samples):
            try:
                gene_tpm[gene_id][s] += float(tpm_values[j])
            except ValueError:
                continue
        if i % 50000 == 0:
            print(f"Processed {i} transcripts...")

print(f"Total genes: {len(gene_tpm)}, Total samples: {len(samples)}")

# Write gene TPM matrix
print("Writing gene_tpm_matrix.tsv...")
with open(gene_matrix_out, 'w') as f:
    f.write("gene_id\t" + "\t".join(samples) + "\n")
    for gene_id in sorted(gene_tpm.keys()):
        vals = [str(round(gene_tpm[gene_id][s], 6)) for s in samples]
        f.write(f"{gene_id}\t" + "\t".join(vals) + "\n")

# Write averaged TPM
print("Writing gene_tpm.tsv...")
with open(gene_avg_out, 'w') as f:
    f.write("gene_id\tavg_tpm\n")
    for g in sorted(gene_tpm.keys()):
        vals = list(gene_tpm[g].values())
        avg = sum(vals) / len(vals)
        f.write(f"{g}\t{avg}\n")

print("Gene-level TPM generation complete.")
EOF

chmod +x "${WORKDIR}/scripts/create_gene_tpm.py"
python3 "${WORKDIR}/scripts/create_gene_tpm.py"

echo ""
echo "Verifying gene-level output..."
ls -lh "${OUTPUT_DIR}/gene_tpm_matrix.tsv" "${OUTPUT_DIR}/gene_tpm.tsv"
head -3 "${OUTPUT_DIR}/gene_tpm_matrix.tsv" | cut -f1-5
head -3 "${OUTPUT_DIR}/gene_tpm.tsv"
echo ""
echo "=== ALL TPM FILES GENERATED SUCCESSFULLY ==="
