#!/bin/bash
set -euo pipefail

# Dorado basecalling (sup + m6A + polyA estimation) and genome alignment (GRCm39)
# Input:  config/cortex_samples.tsv (sample_id, pod5_dir, out_dir, genome_fa)
# Output: <out_dir>/<sample>.dorado_polyA_m6A_unmapped.bam
#         <out_dir>/<sample>.dorado_polyA_m6A_all_reads.fastq
#         <out_dir>/<sample>.genome.primary.sorted.bam (+ index)

SAMPLES_TSV="config/cortex_samples.tsv"

source /etc/profile.d/modules.sh
module load samtools
module load minimap2

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample pod5_dir out_dir genome_fa; do
    mkdir -p "$out_dir"

    dorado basecaller sup,m6A "$pod5_dir" --estimate-poly-a \
        > "${out_dir}/${sample}.dorado_polyA_m6A_unmapped.bam"

    samtools fastq -T "*" \
        "${out_dir}/${sample}.dorado_polyA_m6A_unmapped.bam" \
        > "${out_dir}/${sample}.dorado_polyA_m6A_all_reads.fastq"

    minimap2 -a -x splice -y -k 14 -t 104 "$genome_fa" \
        "${out_dir}/${sample}.dorado_polyA_m6A_all_reads.fastq" \
        | samtools view -@ 104 -b \
        | samtools sort -@ 104 \
        -o "${out_dir}/${sample}.genome.primary.sorted.bam"

    samtools index "${out_dir}/${sample}.genome.primary.sorted.bam"
done
