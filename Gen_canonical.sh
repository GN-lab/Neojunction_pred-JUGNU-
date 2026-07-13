#!/usr/bin/env bash
#SBATCH --job-name=Canonical_Genome
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8042
#SBATCH --time=48:00:00
#SBATCH --output=logs/Canonical_Genome_%j.log

module load STAR

WORKDIR="$(pwd)"
export HOME="/data"

# Create under your current working directory
STAR_INDEX="${WORKDIR}/star_index_GRCh38"
mkdir -p "${STAR_INDEX}"

STAR \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles ../nextflow_rnaseq/data/Ref/genome.fa \
  --sjdbGTFfile     ../nextflow_rnaseq/data/Ref/genome.gtf \
  --sjdbOverhang    100 \
  --runThreadN      8
