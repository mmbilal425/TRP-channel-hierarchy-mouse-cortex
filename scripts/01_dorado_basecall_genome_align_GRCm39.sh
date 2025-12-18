#!/bin/bash
set -euo pipefail

SAMPLES_TSV="config/cortex_samples.tsv"

module load samtools
module load minimap2

# dorado should be available on PATH (as you used with export PATH in your runs)
# If you want to hardcode: export PATH="$PATH:/g/data/lf10/mb1232/apps/dorado-0.7.0-linux-x64/bin"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample data_dir wd genome txfa; do
  mkdir -p "$wd" 2>/dev/null

  ##################################################
  # basecall
  ##################################################
  dorado basecaller sup,m6A "$data_dir" --estimate-poly-a > "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam"

  # preview
  samtools view "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam" | head

  # bam -> fastq
  samtools fastq -T "*" "${wd}/${sample}_dorado_polyA_m6A_unmapped.bam" > "${wd}/${sample}_dorado_polyA_m6A_all_reads.fastq"

  ##################################################
  # map to genome
  ##################################################
  minimap2 -a -x splice -y -k 14 -t 104 "${genome}" "${wd}/${sample}_dorado_polyA_m6A_all_reads.fastq" \
    | samtools view -@ 104 -b > "${wd}/${sample}_total_genome_alignments.bam"

  # count
  samtools view -@ 104 -c "${wd}/${sample}_total_genome_alignments.bam"

  # sort + index
  samtools sort "${wd}/${sample}_total_genome_alignments.bam" > "${wd}/${sample}_total_genome_alignments_sorted.bam"
  samtools index "${wd}/${sample}_total_genome_alignments_sorted.bam"

  # primary alignments
  samtools view -@ 104 -b -F 260 "${wd}/${sample}_total_genome_alignments.bam" | samtools sort -@ 104 > "${wd}/${sample}_primary_genome_alignments.bam"
  samtools index "${wd}/${sample}_primary_genome_alignments.bam"
  samtools view -@ 104 -c "${wd}/${sample}_primary_genome_alignments.bam"

done
