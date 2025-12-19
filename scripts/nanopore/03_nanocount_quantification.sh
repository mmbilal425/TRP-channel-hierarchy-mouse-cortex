#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------
# NanoCount transcript-level quantification
# (Nanopore direct RNA; GRCm39 / Ensembl r115)
# --------------------------------------------

source /g/data/lf10/mb1232/apps/nanocount/nanocount_env/bin/activate

BASE="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics"
RESULTS="${BASE}/nanocount/results"
LOGS="${BASE}/nanocount/log"

mkdir -p "$LOGS"

for sid in cortex_rep1 cortex_rep2 cortex_rep3; do
    IN="${RESULTS}/${sid}/${sid}.txome.sorted.bam"
    OUT="${RESULTS}/${sid}/${sid}.nanocount_transcript.tsv"
    LOG="${LOGS}/${sid}.nanocount.log"

    echo "[INFO] Running NanoCount for ${sid}"
    python3 "$(which NanoCount)" \
        -i "$IN" \
        -o "$OUT" \
        2>&1 | tee "$LOG"
done

deactivate

echo "[OK] NanoCount quantification completed"
