#!/usr/bin/env bash
set -euo pipefail

# === r115 setup ===
export THREADS=28
export SALMON=/apps/salmon/1.10.1/bin/salmon

# refs (note the r115 cDNA fasta)
export GTF=/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz
export TXOME=/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.cdna.all.r115.fa
export GENOME=/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa

# new output root to keep r115 separate
export OUTDIR=/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115
mkdir -p "$OUTDIR" && cd "$OUTDIR"
