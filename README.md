# TRP-channel-hierarchy-mouse-cortex

Reproducible analysis pipeline for multi-omic profiling of Transient Receptor Potential (TRP) channel expression in adult mouse cortex and dorsal root ganglion (DRG).

---

## Preprint

This repository accompanies the following manuscript:

Bilal M, et al.  
*Integrated transcriptomics and proteomics define the TRP channel hierarchy in mouse cortex*.  
bioRxiv. 2026. https://doi.org/10.64898/2026.04.07.716663  

This repository contains all scripts used to generate the analyses and figures presented in the manuscript.

---

## Overview

This project performs integrated multi-omic profiling of TRP channel expression across mouse nervous system tissues using:

- Illumina RNA-seq (Salmon)  
- Nanopore direct RNA sequencing (NanoCount / IsoQuant)  
- m6A epitranscriptomic analysis (modkit)  
- qPCR validation  
- Proteomics (LC–MS/MS-based protein quantification)  
- Downstream manuscript figure generation  

The repository provides all scripts required to reproduce gene-level, isoform-level, and protein-level analyses presented in the manuscript.

Large raw sequencing and proteomics datasets are stored on HPC systems and are not included here.

---

## Data availability

RNA sequencing data are available in GEO under accession **GSE315051**, Link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315051.  
Proteomics data are available in PRIDE under accession **PXD074939**.

---

## Repository structure

### `config/`
- TRP gene ID ↔ gene symbol mappings  
- Sample metadata  
- `project_paths.template.sh` (template for local path configuration)

---

### `scripts/` (upstream pipelines)

#### `scripts/illumina/`
- Adapter trimming  
- Salmon index generation  
- Salmon quantification  
- Merging replicate-level TPM tables  

#### `scripts/nanopore/`
- Alignment (minimap2)  
- Gene-level quantification (NanoCount)  
- Isoform reconstruction (IsoQuant)  

---

### `downstream/` (manuscript analyses)

#### `downstream/expression/`

**tables/**
- Gene-level TPM and read-count generation  

**figures/**
- Main figures (Fig. 1)  
- qPCR validation (Fig. 2, Fig. 4)  
- Correlation analyses  
- Supplementary Figures 2–4  

---

#### `downstream/m6A/`

- modkit pileup import  
- Motif annotation (DRACH, RAC)  
- Modified site calling  
- Gene assignment  
- Supplementary Figure 6 panels  

---

#### `downstream/proteomics/`

Proteomics analysis is maintained as a modular workflow independent of transcriptomics.

**GO_analysis/**
- GO:CC enrichment and bubble plots  

**IP_analysis/**
- TRP immunoprecipitation rerun analysis and peptide heatmaps  

**heatmap/**
- Protein intensity heatmaps across cortex and peripheral organs  

**proteomics_figures/**
- Additional manuscript figures  

**heatmap_log/**
- Log10-transformed TRP intensity heatmaps (Fig. 3E, Supplementary Fig. 8D)  

**protein_per_group/**
- Protein counts per sample and PCA legend generation  

---

## Key inputs

Most downstream scripts require:

- Ensembl r115 annotation  
  - `Mus_musculus.GRCm39.115.gtf.gz`  
- TRP gene ID + gene symbol mapping  
- Merged Salmon TPM tables (cortex + DRG)  
- NanoCount transcript-level outputs  
- IsoQuant isoform models  
- Annotated m6A site tables  
- Processed proteomics intensity tables  

---

## Expression analyses

### NanoCount (Nanopore, cortex)

Located in: `downstream/expression/tables/`

- `00_nanocount_allgenes_TPM_cortex_r115.py`  
  → Gene-level TPM (replicates + mean)  

- `01_nanocount_allgenes_readcounts_cortex_r115.py`  
  → Gene-level read counts (replicates + mean)  

Used in:
- Fig. 1B  
- Supplementary Fig. 2a  
- Supplementary Fig. 3a  
- Supplementary Fig. 4  

---

### Salmon (Illumina, cortex + DRG)

Merged outputs include:

- Cortex replicate TPMs  
- Cortex mean TPM  
- DRG TPM  
- Gene ID annotation  

Used in:
- Fig. 1C–D  
- Supplementary Fig. 2b–d  
- Supplementary Fig. 3b–d  

---

### Isoform analysis (IsoQuant)

- Known TRP isoforms  
- Novel TRP isoforms  
- Long-read splice validation  

Used in:
- Supplementary Figure 4  

---

## m6A analysis

Located in: `downstream/m6A/`

Pipeline:
- Import modkit pileup  
- Annotate motifs  
- Call modified sites  
- Assign to genes  
- Generate figures  

Outputs:
- Stoichiometry distribution (S6a)  
- Motif enrichment (S6b)  
- Sites per gene (S6c)  

---

## Proteomics analysis

Located in: `downstream/proteomics/`

Includes:

- Membrane-aware protein extraction comparisons  
- TRP protein detection across tissues  
- PCA and protein-group analysis  
- GO enrichment (cellular component)  
- TRP-focused heatmaps  

All proteomics figures are generated from processed intensity tables included in the repository.

---

## Software

- Dorado v0.7.0  
- minimap2 v2.28  
- Salmon v1.10.1  
- NanoCount v1.2.1  
- IsoQuant v3.7.0  
- modkit  
- Python ≥ v3.9  
- R  

---

## Citation and licensing

This repository contains code associated with the manuscript:

Bilal M, et al.
*Integrated transcriptomics and proteomics define the TRP channel hierarchy in mouse cortex*.  
bioRxiv. 2026. https://doi.org/10.64898/2026.04.07.716663  

Repository: https://github.com/mmbilal425/TRP-channel-hierarchy-mouse-cortex  

If you use this code or analysis workflow, please cite the manuscript above and this repository.

All original code in this repository is released under the MIT License.  
Third-party tools and dependencies retain their respective licenses.
