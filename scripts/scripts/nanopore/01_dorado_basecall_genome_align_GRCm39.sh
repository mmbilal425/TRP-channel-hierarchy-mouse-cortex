#!/bin/bash
set -euo pipefail

# Dorado basecalling (sup,m6A + polyA estimation) and genome alignment (GRCm39)
# Input: config/cortex_samples.tsv (sample_id, pod5_dir, out_dir, genome_fa, tx_fa)

SAMPLES_TSV="config/cortex_samples.tsv"

source /etc/profile.d/modules.sh
module load samtools
module load minimap2

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample data_dir wd genome txfa; do
  mkdir -p "$wd" 2>/dev/null

  ##################################################
  # basecall
  ##################################################
  dorado basecaller sup,m6A "$data_dir" --estimate-poly-a > "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam"

  samtools view "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam" | head

  samtools fastq -T "*" "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam" > "${wd}/${sample}_dorado_polyA_m6A_all_reads.fastq"

  ##################################################
  # map to genome
  ##################################################
  minimap2 -a -x splice -y -k 14 -t 104 "${genome}" "${wd}/${sample}_dorado_polyA_m6A_all_reads.fastq" \
    | samtools view -@ 104 -b > "${wd}/${sample}_total_genome_alignments.bam"

  samtools view -@ 104 -c "${wd}/${sample}_total_genome_alignments.bam"

  samtools sort "${wd}/${sample}_total_genome_alignments.bam" > "${wd}/${sample}_total_genome_alignments_sorted.bam"
  samtools index "${wd}/${sample}_total_genome_alignments_sorted.bam"

  samtools view -@ 104 -b -F 260 "${wd}/${sample}_total_genome_alignments.bam" | samtools sort -@ 104 > "${wd}/${sample}_primary_genome_alignments.bam"
  samtools index "${wd}/${sample}_primary_genome_alignments.bam"

  samtools view -@ 104 -c "${wd}/${sample}_primary_genome_alignments.bam"
done
