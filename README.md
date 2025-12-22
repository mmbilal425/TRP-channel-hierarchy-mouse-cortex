# TRP_mouse_cortex

Scripts used to analyse TRP channel expression in adult mouse cortex and DRG using Illumina RNA-seq (Salmon), Nanopore direct RNA sequencing (NanoCount / IsoQuant), and downstream plotting for main and supplementary figures.

## Overview
This repository accompanies a manuscript describing a multi-omic analysis of TRP channel transcripts (gene- and isoform-level) in adult mouse cortex and dorsal root ganglion (DRG). It contains:
- sample metadata/config files
- scripts to generate gene-level TPM and read-count tables
- plotting scripts used for manuscript figures

## Data availability
- Illumina and Nanopore sequencing data: GEO (accession to be added)
- Proteomics data: PRIDE (accession to be added)

## Repository structure
- `config/`
  - Small config/metadata files (e.g., sample sheets, TRP gene lists)
- `scripts/`
  - Upstream / processing scripts (Illumina and Nanopore pipelines)
- `downstream/`
  - Downstream analyses and figure-generation scripts
  - `downstream/expression/`
    - `tables/` – scripts that generate gene-level tables (TPM, read counts)
    - `figures/` – plotting scripts for main + supplementary expression figures

## Key inputs
Most downstream scripts expect:
- Ensembl r115 GTF: `Mus_musculus.GRCm39.115.gtf.gz`
- TRP gene IDs + symbols 

## Downstream expression (current)
### Tables
Located in `downstream/expression/tables/`

- **NanoCount (Nanopore, cortex)**
  - `00_nanocount_allgenes_TPM_cortex_r115.py`  
    Builds gene-level TPM table across cortex replicates + mean
  - `01_nanocount_allgenes_readcounts_cortex_r115.py`  
    Builds gene-level read-count table across cortex replicates + mean

- **Salmon (Illumina, cortex + DRG)**
  - Scripts that merge Salmon outputs and produce cortex mean + DRG columns
  - Outputs include merged TPM and read-count TSVs used by figure scripts

### Figures
Located in `downstream/expression/figures/`

- Main Figure 1 (expression panels)
  - Fig 1B: NanoCount cortex TPM
  - Fig 1C: Salmon cortex mean TPM
  - Fig 1D: Salmon DRG TPM

- Supplementary figures (examples)
  - Supp Fig S1: read-count comparisons (NanoCount + Salmon)
  - Supp Fig S2: TPM family pie, Venn overlap, log2FC heatmap
  - (additional scripts may be added as analysis expands)

## Software
Analyses were performed using:
- Salmon
- NanoCount
- IsoQuant
- Python (3.x) with pandas / numpy / matplotlib / seaborn
- R (where applicable)

## Notes
This repository stores **scripts + small config files**. Large intermediate results and raw data are kept on HPC storage and referenced by absolute paths inside scripts where appropriate.
