#!/usr/bin/env bash
set -euo pipefail

# Purpose:
#   Dorado basecalling (sup,m6A + polyA estimation) and genome alignment
#   for Nanopore direct RNA-seq mouse cortex replicates.
#
# Input:
#   config/cortex_samples.tsv
#
# Output:
#   <out_dir>/{sample}.dorado_polyA_m6A.unmapped.bam
#   <out_dir>/{sample}.dorado_polyA_m6A.fastq
#   <out_dir>/{sample}.genome.sorted.bam
#   <out_dir>/{sample}.genome.primary.sorted.bam
#
# Notes:
#   Processing/QC only. Transcriptome mapping and NanoCount are handled separately.
