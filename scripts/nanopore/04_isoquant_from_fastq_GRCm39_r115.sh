#!/bin/bash
set -euo pipefail

# --------- paths (your inputs) ----------
REF="/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_GZ="/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz"

FQ1="/g/data/lf10/mb1232/analysis/mouse_cortex_rep1_analysis/mouse_cortex_rep1_dorado_polyA_m6A_all_reads.fastq"
FQ2="/g/data/lf10/mb1232/analysis/mouse_cortex_rep2_analysis/mouse_cortex_rep2_analysis_dorado_polyA_m6A_all_reads.fastq"
FQ3="/g/data/qq78/mb1232/2024_09_17_nanopore_rep3/analysis/mouse_cortex_rep3_analysis_dorado_polyA_m6A_all_reads.fastq"

OUT="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq"
PREFIX="Mouse_Cortex_dRNA"
THREADS=28

source /g/data/qq78/mb1232/apps/miniconda3/etc/profile.d/conda.sh
conda activate isoquant

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

mkdir -p "${OUT}/refs" "${OUT}/tmp"

GTF_LOCAL="${OUT}/refs/Mus_musculus.GRCm39.115.gtf"
if [[ "${GTF_GZ}" == *.gz ]]; then
  echo "Decompressing GTF → ${GTF_LOCAL}"
  gzip -cd "${GTF_GZ}" > "${GTF_LOCAL}"
else
  echo "Copying GTF → ${GTF_LOCAL}"
  cp -f "${GTF_GZ}" "${GTF_LOCAL}"
fi

rm -f "${OUT}/refs/Mus_musculus.GRCm39.115.db" || true

LOG="${OUT}/isoquant_run_$(date +%F_%H%M)_${THREADS}t.log"
cd "${OUT}"

python3 /g/data/qq78/mb1232/apps/miniconda3/envs/isoquant/bin/isoquant.py \
  --data_type nanopore \
  --stranded forward \
  --reference "${REF}" \
  --genedb "refs/$(basename "${GTF_LOCAL}")" \
  --complete_genedb \
  --fastq "${FQ1}" "${FQ2}" "${FQ3}" \
  --labels rep1 rep2 rep3 \
  --read_group file_name \
  --output "${OUT}" \
  --prefix "${PREFIX}" \
  --threads "${THREADS}" \
  --sqanti_output \
  --count_exons \
  --report_novel_unspliced true \
  --force \
  2>&1 | tee -a "${LOG}"

conda deactivate
echo "Done. Log: ${LOG}"
