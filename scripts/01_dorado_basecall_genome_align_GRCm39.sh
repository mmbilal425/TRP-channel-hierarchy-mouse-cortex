#!/usr/bin/env bash
set -euo pipefail

# Purpose: Dorado basecalling (sup,m6A + polyA estimation) and genome alignment for cortex dRNA replicates.
# Input:  config/cortex_samples.tsv (sample_id, pod5_dir, out_dir, genome_fa)
# Output: <out_dir>/{sample}.dorado.unmapped.bam, .fastq, genome.bam, genome.sorted.bam, genome.primary.bam (+ indexes)
# Notes:  This script is for processing/QC; NanoCount uses transcriptome mapping (see 02/03).

source /etc/profile.d/modules.sh
module load samtools
module load minimap2

# ---- EDIT THESE PATHS ON YOUR SYSTEM ----
DORADO_BIN="/g/data/lf10/mb1232/apps/dorado-0.7.0-linux-x64/bin/dorado"
THREADS_ALIGN=28   # for minimap2+samtools
THREADS_BAM=28
KMER=14

SAMPLES_TSV="config/cortex_samples.tsv"

# ----------------------------------------
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample pod5 outdir genome txfa; do
  mkdir -p "$outdir"

  unmapped_bam="$outdir/${sample}.dorado_polyA_m6A.unmapped.bam"
  fastq="$outdir/${sample}.dorado_polyA_m6A.fastq"
  genome_bam="$outdir/${sample}.genome.bam"
  genome_sorted="$outdir/${sample}.genome.sorted.bam"
  genome_primary="$outdir/${sample}.genome.primary.sorted.bam"

  echo "==== [$sample] Dorado basecall (sup,m6A + polyA) ===="
  "$DORADO_BIN" basecaller sup,m6A "$pod5" --estimate-poly-a > "$unmapped_bam"

  echo "==== [$sample] FASTQ export ===="
  samtools fastq -T "*" "$unmapped_bam" > "$fastq"

  echo "==== [$sample] Genome align (splice) ===="
  minimap2 -a -x splice -y -k "$KMER" -t "$THREADS_ALIGN" "$genome" "$fastq" \
    | samtools view -@ "$THREADS_BAM" -b -o "$genome_bam" -

  echo "==== [$sample] Sort + index genome alignments ===="
  samtools sort -@ "$THREADS_BAM" -o "$genome_sorted" "$genome_bam"
  samtools index -@ "$THREADS_BAM" "$genome_sorted"

  echo "==== [$sample] Primary genome alignments (-F 260) ===="
  samtools view -@ "$THREADS_BAM" -b -F 260 "$genome_bam" \
    | samtools sort -@ "$THREADS_BAM" -o "$genome_primary" -
  samtools index -@ "$THREADS_BAM" "$genome_primary"

  echo "==== [$sample] Counts ===="
  echo "total_alignments: $(samtools view -c "$genome_bam")"
  echo "primary_alignments: $(samtools view -c "$genome_primary")"
done
