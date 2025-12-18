#!/bin/bash
set -euo pipefail

# Map Dorado FASTQs to Ensembl r115 transcriptome (cDNA) for NanoCount
# Input:  config/cortex_samples.tsv (sample_id, pod5_dir, out_dir, genome_fa, tx_fa)
# Output: <out_dir>/${sample}.txome.sorted.bam (+ .bai)

SAMPLES_TSV="config/cortex_samples_nanopore.tsv"

source /etc/profile.d/modules.sh
module load minimap2
module load samtools

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample data_dir wd genome txfa; do

  FQ="${wd}/${sample}_dorado_polyA_m6A_all_reads.fastq"
  BAM="${wd}/${sample}.txome.sorted.bam"

  minimap2 -t 28 -ax map-ont -p 0.8 -N 50 "$txfa" "$FQ" \
    | samtools sort -@ 28 -o "$BAM" -

  samtools index -@ 28 "$BAM"
done
