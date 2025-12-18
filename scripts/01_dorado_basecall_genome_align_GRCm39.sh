#!/usr/bin/env bash
set -euo pipefail

# Purpose: Dorado basecalling (sup,m6A + polyA estimation) and genome alignment for cortex dRNA replicates.
# Input:   config/cortex_samples.tsv (sample_id, pod5_dir, out_dir, genome_fa, tx_fa)
# Output:  <out_dir>/{sample}.dorado_polyA_m6A.unmapped.bam
#          <out_dir>/{sample}.dorado_polyA_m6A.fastq
#          <out_dir>/{sample}.genome.bam
#          <out_dir>/{sample}.genome.sorted.bam (+ .bai)
#          <out_dir>/{sample}.genome.primary.sorted.bam (+ .bai)
# Notes:   Processing/QC only. NanoCount uses transcriptome mapping (see 02/03).

# ----------------------------- configuration -----------------------------
SAMPLES_TSV="config/cortex_samples.tsv"
LOGDIR="logs/01_dorado_genome_align"
mkdir -p "$LOGDIR"

# Dorado binary (preferred explicit path). If empty, script will try `dorado` from PATH.
DORADO_BIN="/g/data/lf10/mb1232/apps/dorado-0.7.0-linux-x64/bin/dorado"

# Threads
THREADS_ALIGN=28   # minimap2
THREADS_BAM=28     # samtools
KMER=14            # minimap2 splice k-mer (as used in your runs)

# ----------------------------- helpers -----------------------------------
have_cmd () { command -v "$1" >/dev/null 2>&1; }

load_modules_if_available () {
  # "C" behaviour: if modules exist, use them; otherwise rely on PATH.
  if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh || true
    if have_cmd module; then
      module -q load samtools || true
      module -q load minimap2 || true
    fi
  fi
}

resolve_dorado () {
  if [ -n "${DORADO_BIN}" ] && [ -x "${DORADO_BIN}" ]; then
    echo "${DORADO_BIN}"
    return 0
  fi
  if have_cmd dorado; then
    command -v dorado
    return 0
  fi
  echo "[ERROR] Dorado not found. Set DORADO_BIN to an executable dorado, or add dorado to PATH." >&2
  exit 1
}

# ----------------------------- main --------------------------------------
load_modules_if_available

if ! have_cmd samtools; then
  echo "[ERROR] samtools not found (module load failed and not in PATH)." >&2
  exit 1
fi
if ! have_cmd minimap2; then
  echo "[ERROR] minimap2 not found (module load failed and not in PATH)." >&2
  exit 1
fi

DORADO="$(resolve_dorado)"

# Read TSV robustly (skip header, keep tab separation)
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample pod5 outdir genome txfa; do
  # basic sanity checks
  if [ -z "${sample}" ] || [ -z "${pod5}" ] || [ -z "${outdir}" ] || [ -z "${genome}" ]; then
    echo "[WARN] Skipping incomplete row in $SAMPLES_TSV: sample='$sample'" >&2
    continue
  fi

  mkdir -p "$outdir"
  sample_log="$LOGDIR/${sample}.log"

  unmapped_bam="$outdir/${sample}.dorado_polyA_m6A.unmapped.bam"
  fastq="$outdir/${sample}.dorado_polyA_m6A.fastq"
  genome_bam="$outdir/${sample}.genome.bam"
  genome_sorted="$outdir/${sample}.genome.sorted.bam"
  genome_primary="$outdir/${sample}.genome.primary.sorted.bam"

  {
    echo "==== [$sample] START ===="
    echo "[INFO] pod5_dir: $pod5"
    echo "[INFO] out_dir:  $outdir"
    echo "[INFO] genome:   $genome"
    echo "[INFO] dorado:   $DORADO"
    echo "[INFO] minimap2: $(command -v minimap2)"
    echo "[INFO] samtools: $(command -v samtools)"

    echo "---- [$sample] Dorado basecall (sup,m6A + polyA) ----"
    "$DORADO" basecaller sup,m6A "$pod5" --estimate-poly-a > "$unmapped_bam"

    echo "---- [$sample] FASTQ export ----"
    samtools fastq -T "*" "$unmapped_bam" > "$fastq"

    echo "---- [$sample] Genome alignment (splice) ----"
    minimap2 -a -x splice -y -k "$KMER" -t "$THREADS_ALIGN" "$genome" "$fastq" \
      | samtools view -@ "$THREADS_BAM" -b -o "$genome_bam" -

    echo "---- [$sample] Sort + index ----"
    samtools sort -@ "$THREADS_BAM" -o "$genome_sorted" "$genome_bam"
    samtools index -@ "$THREADS_BAM" "$genome_sorted"

    echo "---- [$sample] Primary alignments (-F 260) ----"
    samtools view -@ "$THREADS_BAM" -b -F 260 "$genome_bam" \
      | samtools sort -@ "$THREADS_BAM" -o "$genome_primary" -
    samtools index -@ "$THREADS_BAM" "$genome_primary"

    echo "---- [$sample] Alignment counts ----"
    echo "total_alignments:   $(samtools view -c "$genome_bam")"
    echo "primary_alignments: $(samtools view -c "$genome_primary")"

    echo "==== [$sample] DONE ===="
  } 2>&1 | tee "$sample_log"

done
