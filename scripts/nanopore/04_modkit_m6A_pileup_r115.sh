#!/bin/bash
set -euo pipefail

# Nanopore m6A detection using modkit
# Reference: Ensembl GRCm39 r115
# Output: genome-wide m6A pileup BED

module load samtools
module load modkit

REF="/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# merged modBAM containing MM/ML tags
BAM="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/<MERGED_MODBAM>.bam"

OUTDIR="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/m6a_modkit_output"
mkdir -p "${OUTDIR}"

OUTBED="${OUTDIR}/mouse_cortex_merged_pileup_filt0.5_mod0.99_global0.9_allMods.bed"

modkit pileup \
  --ref "${REF}" \
  --filter-threshold 0.5 \
  --mod-thresholds a:0.99 \
  "${BAM}" \
  "${OUTBED}"

echo "[OK] m6A pileup written to ${OUTBED}"
