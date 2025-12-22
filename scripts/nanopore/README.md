# Nanopore processing (direct RNA; GRCm39 / Ensembl r115)

This folder contains the **primary Nanopore processing pipeline** used to generate:
- genome-aligned BAMs (QC / mapping)
- transcriptome-aligned BAMs (quantification)
- NanoCount transcript quantification (TPM / counts)
- IsoQuant isoform discovery
- modkit m6A pileup BED for downstream epitranscriptomics

## Run order

1) **01_dorado_basecall_genome_align_GRCm39.sh**
   - Input: POD5/FAST5 (raw Nanopore)
   - Output: genome-aligned BAM(s) (GRCm39)

2) **02_map_to_transcriptome_r115.sh**
   - Input: reads / BAM from step 1
   - Output: transcriptome-aligned BAM(s) (Ensembl r115 transcriptome)

3) **03_nanocount_quantification.sh**
   - Input: transcriptome-sorted BAM(s) from step 2
   - Output: NanoCount transcript tables (`*.nanocount_transcript.tsv`)
   - Note: these outputs are the source for all Nanopore-derived TPM figures.

4) **04_isoquant_from_fastq_GRCm39_r115.sh**
   - Input: FASTQ (or BAM depending on your workflow) + r115 GTF
   - Output: isoform annotations and classification tables

5) **05_modkit_m6A_pileup_r115.sh**
   - Input: aligned BAM(s) + reference
   - Output: modkit pileup BED (used by `downstream/m6A/`)

## Downstream analyses

- TPM tables/figures: `downstream/expression/`
- m6A analysis + figures: `downstream/m6A/`
- Isoform summaries/figures: `downstream/isoforms/`
