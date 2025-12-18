# TRP_mouse_cortex

Scripts used to analyse TRP channel expression in adult mouse cortex using
Illumina RNA-seq, Nanopore direct RNA sequencing, qPCR, and membrane-aware
proteomics.

## Overview
This repository accompanies the manuscript describing a multi-omic analysis
of TRP channel transcripts and proteins in adult mouse cortex.

## Data availability
Raw sequencing data are available from GEO (accession to be added).
Proteomics data will be deposited in PRIDE.

## Scripts
- Illumina_RNAseq_TRP_analysis.sh  
  Quantification of TRP transcript abundance from Illumina RNA-seq data.

- Nanopore_dRNA_TRP_isoform_analysis.sh  
  Isoform-level analysis of TRP transcripts from Nanopore direct RNA sequencing.

- NanoCount_TRP_gene_TPM.py  
  Gene-level TPM summarisation from NanoCount output.

- plot_Fig1_TRP_expression.py  
  Generates Figure 1.

- plot_Fig2_qPCR_validation.py  
  Generates Figure 2.

- plot_Fig3_proteomics.R  
  Generates Figure 3.

## Software
Analyses were performed using Salmon, NanoCount, IsoQuant, Python 3.10,
and R 4.3.
