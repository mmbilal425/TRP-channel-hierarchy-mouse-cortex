# Illumina RNA-seq processing (mouse cortex & DRG)

This directory contains the Illumina short-read RNA-seq processing
pipeline used in this study (GENEWIZ, NovaSeq 6000).

Reference genome and annotation:
- Mus musculus GRCm39
- Ensembl release 115

## Pipeline overview

- **01_trimming/** – Adapter and quality trimming (FASTX toolkit)
- **02_salmon/** – Transcript- and gene-level quantification (Salmon)
  - **01_index/** – Decoy-aware Salmon index (Ensembl r115)
  - **02_quant/** – Per-sample transcript and gene quantification
  - **03_merge/** – Merge gene-level TPMs across samples
