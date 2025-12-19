#!/bin/bash
set -euo pipefail

module load samtools

# ---- inputs ----
REF="/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# merged aligned modBAM (MM/ML tags present; same one used to generate your existing pileup bed)
BAM="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/<PUT_YOUR_MERGED_MODBAM_HERE>.bam"

OUTDIR="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/m6a_modkit_output"
mkdir -p "$OUTDIR"

OUTBED="${OUTDIR}/mouse_cortex_merged_pileup_filt0.5_mod0.99_global0.9_allMods.bed"

# ---- run ----
# thresholds are encoded in your output filename; keep them consistent
modkit pileup \
  --ref "$REF" \
  --filter-threshold 0.5 \
  --mod-thresholds a:0.99 \
  --output "$OUTBED" \
  "$BAM"

echo "[OK] Wrote: $OUTBED"
