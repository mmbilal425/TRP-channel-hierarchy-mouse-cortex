# TRP_mouse_cortex

Reproducible analysis pipeline for Transcriptomic and Proteomic Profiling of TRP Channels in Mouse Cortex and Other Peripheral Organs 

This repository accompanies a manuscript describing integrated transcriptomic, isoform-level, epitranscriptomic, and qPCR validation analyses of TRP channels using:

- Illumina RNA-seq (Salmon)
- Nanopore direct RNA sequencing (NanoCount / IsoQuant)
- m6A modification analysis (modkit)
- qPCR validation
- Downstream manuscript figure generation

Large raw sequencing data and intermediate files are stored on HPC systems and are not included in this repository.

---

## Overview

This project performs multi-platform profiling of TRP channel expression in adult mouse cortex and DRG, including:

- Gene-level TPM and read-count quantification
- Cross-platform comparison (Illumina vs Nanopore)
- Isoform reconstruction and quantification
- m6A site detection and motif analysis
- qPCR validation
- Manuscript-ready figure generation

All manuscript figures are generated directly from scripts contained in this repository.

---

## Data availability

- Illumina RNA-seq and Nanopore direct RNA-seq: GEO (accession to be added)
- Proteomics data: PRIDE (accession to be added)

---

## Repository structure

- `config/`
  - TRP gene ID ↔ gene symbol mappings
  - Sample metadata
  - `project_paths.template.sh` (template for defining local paths)

- `scripts/`
  - Upstream processing pipelines
  - `scripts/illumina/`
    - Adapter trimming
    - Salmon index generation
    - Salmon quantification
    - Merge replicate-level gene TPM tables
  - `scripts/nanopore/`
    - Alignment (minimap2)
    - Gene-level quantification (NanoCount)
    - Isoform reconstruction (IsoQuant)

- `downstream/`
  - Manuscript-facing analyses and figure-generation scripts

  - `downstream/expression/`
    - `tables/`
      - Scripts generating NanoCount gene-level TPM and read-count tables (cortex)
      - Scripts generating merged Salmon TPM tables (cortex mean + DRG)
    - `figures/`
      - Main expression panels (Fig 1)
      - qPCR validation panels (Fig 2 and Fig 4)
      - Correlation plots
      - Supplementary expression figures (S2–S4)

  - `downstream/m6A/`
    - modkit pileup import
    - Motif annotation (DRACH, RAC)
    - Modified site calling
    - Gene assignment
    - Supplementary Figure S6 panels

---

## Key inputs

Most downstream scripts require:

- Ensembl r115 annotation:
  `Mus_musculus.GRCm39.115.gtf.gz`

- TRP gene ID + gene symbol mapping

- Merged Salmon gene-level TPM tables

- NanoCount transcript-level quantification outputs

- IsoQuant isoform models (for isoform panels)

- Annotated m6A site tables (for Supplementary Fig S6)

---

## Expression analyses

### NanoCount (Nanopore, cortex)

Located in:
`downstream/expression/tables/`

- `00_nanocount_allgenes_TPM_cortex_r115.py`
  Builds gene-level TPM table across cortex replicates + mean

- `01_nanocount_allgenes_readcounts_cortex_r115.py`
  Builds gene-level read-count table across cortex replicates + mean

Used in:
- Fig 1B
- Supp Fig S2a
- Supp Fig S3a
- Supp Fig S4 (isoform panels)

---

### Salmon (Illumina, cortex + DRG)

Merged outputs include:

- Cortex replicate TPMs
- Cortex mean TPM
- DRG TPM
- Gene ID annotation (Ensembl r115)

Used in:
- Fig 1C–D
- Supp Fig S2b,c,d
- Supp Fig S3b,d

---

### Isoform analysis (IsoQuant)

- Known TRP isoforms
- Novel TRP isoforms
- Long-read splice validation
- Isoform TPM panels (Supp Fig S4)

---

## m6A analysis (Nanopore direct RNA)

Located in:
`downstream/m6A/`

Pipeline steps:

1. Import modkit pileup output
2. Attach motif annotation (DRACH, RAC)
3. Call modified sites
4. Assign sites to genes
5. Generate Supplementary Fig S6 panels

Figure outputs include:

- Stoichiometry distribution (S6a)
- Ranked motif enrichment (S6b)
- Modified m6A sites per TRP gene (S6c)

---

## Software

Analyses were performed using:

- Dorado v0.7.0
- minimap2 v2.28
- Salmon v1.10.1
- NanoCount v1.2.1
- IsoQuant v3.7.0
- modkit (Nanopore modification calling)
- Python ≥ 3.9
- R (for selected plotting scripts)

Exact versions are documented in the manuscript.

---

## Reproducibility and path management

This repository contains:

- Scripts
- Metadata
- Configuration templates

Large raw data and intermediate files are stored on HPC systems.

Local paths must be defined in:

`config/project_paths.sh`

Example template:

RAW_DATA="/path/to/raw_fastq"
ANALYSIS="/path/to/analysis_output"
REFS="/path/to/reference_genomes"
TOOLS_BASE="/path/to/tools"
THREADS=25

This file is excluded from version control to prevent exposure of institutional storage paths.

---

## Notes

- All figure scripts produce publication-ready vector PDFs.
- `.py` scripts are the authoritative reproducible source.
- Large intermediate data files are intentionally excluded.
- The repository structure mirrors the logical progression of the manuscript.
