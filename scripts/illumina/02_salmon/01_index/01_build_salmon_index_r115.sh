#!/usr/bin/env bash
set -euo pipefail

# Work in the r115 output folder
cd "$OUTDIR"

catz() { case "$1" in *.gz) zcat "$1" ;; *) cat "$1" ;; esac; }

# 1) Decoy list
catz "$GENOME" | awk '/^>/{print substr($1,2)}' > decoys.txt

# 2) Build gentrome (r115 cDNA + genome)
rm -f gentrome.fa
catz "$TXOME"  > gentrome.fa
catz "$GENOME" >> gentrome.fa

# 3) Index
export INDEX_DIR=${OUTDIR}/salmon_index_mm39_r115_k25_fullDecoy
mkdir -p "$INDEX_DIR"

"$SALMON" index \
  -t gentrome.fa \
  -d decoys.txt \
  -p "$THREADS" \
  -k 25 \
  -i "$INDEX_DIR"
