#!/usr/bin/env bash
set -euo pipefail

# mouse_cortex_rep3 trimming (fastx_toolkit + fastq-pair)

raw_data="/g/data/lf10/raw_data/GENEWIZ_2024-Jun_Ehsan"
analysis="/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina"
fastq_out="${analysis}/fastq"
logs="${analysis}/logs"

THREADS=25
lib="cortex-rep3"

mkdir -p "$fastq_out" "$logs"

export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin:$PATH"
export PATH="/g/data/lf10/tools/fastq-pair/bin:$PATH"

need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH"; exit 2; }; }
need_cmd fastx_clipper
need_cmd fastq_quality_filter
need_cmd fastq_pair
if command -v pigz >/dev/null 2>&1; then GZIP_CMD="pigz -f -p ${THREADS}"; else GZIP_CMD="gzip -f"; fi

log_file="${logs}/${lib}.log"

{
  echo "==== [$(date)] TRIM START ${lib} ===="

  fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -i <(zcat "${raw_data}/${lib}_L2_1.fq.gz" | fastq_quality_filter -q 33 -p 40) \
    > "${fastq_out}/${lib}_R1.trimmed.fastq"

  fastx_clipper -l 35 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
    -i <(zcat "${raw_data}/${lib}_L2_2.fq.gz" | fastq_quality_filter -q 33 -p 40) \
    > "${fastq_out}/${lib}_R2.trimmed.fastq"

  fastq_pair \
    "${fastq_out}/${lib}_R1.trimmed.fastq" \
    "${fastq_out}/${lib}_R2.trimmed.fastq"

  for TYPE in R1 R2; do
    ${GZIP_CMD} "${fastq_out}/${lib}_${TYPE}.trimmed.fastq.paired.fq"
    rm -f "${fastq_out}/${lib}_${TYPE}.trimmed.fastq.single.fq" \
          "${fastq_out}/${lib}_${TYPE}.trimmed.fastq"
  done

  echo "==== [$(date)] TRIM DONE ${lib} ===="
} 2>&1 | tee -a "$log_file"
