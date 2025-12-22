#!/usr/bin/env bash
set -euo pipefail

SALMON=/apps/salmon/1.10.1/bin/salmon
THREADS=28

OUTDIR_BASE=/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115
INDEX_DIR=${OUTDIR_BASE}/salmon_index_mm39_r115_k25_fullDecoy

FASTQ_DIR=/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/fastq
GENEMAP_SYMBOL=${OUTDIR_BASE}/tx2gene_symbol.r115.from_fasta.tsv
GENEMAP_ENSID=${OUTDIR_BASE}/tx2gene_ensid.r115.from_fasta.tsv

RESULTS=${OUTDIR_BASE}/results
mkdir -p "$OUTDIR_BASE" "$RESULTS"
cd "$OUTDIR_BASE"

SAMPLES=(cortex-rep1 cortex-rep2 cortex-rep3 DRG-rep1)

for S in "${SAMPLES[@]}"; do
  R1=${FASTQ_DIR}/${S}_R1.trimmed.fastq.paired.fq.gz
  R2=${FASTQ_DIR}/${S}_R2.trimmed.fastq.paired.fq.gz

  echo "=== ${S}: transcript-level (quant.sf) ==="
  /usr/bin/time -v "$SALMON" quant \
    -i "$INDEX_DIR" -l ISR \
    -1 "$R1" -2 "$R2" \
    --validateMappings --seqBias --gcBias \
    -p "$THREADS" \
    -o ${OUTDIR_BASE}/${S}_transcripts

  echo "=== ${S}: gene-level by SYMBOL (quant.genes.sf) ==="
  /usr/bin/time -v "$SALMON" quant \
    -i "$INDEX_DIR" -l ISR \
    -1 "$R1" -2 "$R2" \
    --validateMappings --seqBias --gcBias \
    --geneMap "$GENEMAP_SYMBOL" \
    -p "$THREADS" \
    -o ${OUTDIR_BASE}/${S}_genes_symbol

  echo "=== ${S}: gene-level by ENSEMBL ID (quant.genes.sf) ==="
  /usr/bin/time -v "$SALMON" quant \
    -i "$INDEX_DIR" -l ISR \
    -1 "$R1" -2 "$R2" \
    --validateMappings --seqBias --gcBias \
    --geneMap "$GENEMAP_ENSID" \
    -p "$THREADS" \
    -o ${OUTDIR_BASE}/${S}_genes_ensid

  echo "--- Quick summary: ${S} ---"
  grep -E "Mapping rate|Num. fragments processed" \
    ${OUTDIR_BASE}/${S}_transcripts/logs/salmon_quant.log || true
  head -n 5 ${OUTDIR_BASE}/${S}_genes_symbol/quant.genes.sf || true
done

echo "All done. Outputs in: ${OUTDIR_BASE}/"
