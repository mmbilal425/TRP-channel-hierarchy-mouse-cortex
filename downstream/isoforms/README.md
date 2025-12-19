# Isoform-level analysis of TRP transcripts (IsoQuant)

This folder contains downstream scripts used to identify and characterise
novel TRP channel isoforms from IsoQuant long-read RNA-seq analysis
(adult mouse cortex, Nanopore direct RNA).

All scripts operate on **precomputed IsoQuant outputs** and were executed
sequentially as described below. No re-alignment or IsoQuant reruns are
performed here.

## Input data
- IsoQuant outputs generated from Nanopore direct RNA-seq FASTQ files
- Ensembl GRCm39 release 115 gene annotation
- Curated list of TRP gene Ensembl IDs (`trp_gene_ids.txt`)

## Script order and purpose

**00_filter_novel_TRP_isoforms.py**  
Filters IsoQuant SQANTI-like annotations to retain novel TRP isoforms
classified as `novel_in_catalog` (NIC) or `novel_not_in_catalog` (NNC).

**01_attach_TPM_to_novel_TRP_isoforms.py**  
Joins replicate-level transcript TPMs to the filtered novel TRP isoforms
and computes mean TPM across replicates.

**02_export_novel_TRP_GTF_and_IDs.py**  
Exports transcript IDs and a subset GTF containing only novel TRP isoforms
for genome browser visualisation (IGV).

**03_add_gene_names_to_novel_TRP_table.py**  
Annotates Ensembl gene IDs with gene symbols using the same reference GTF
used for IsoQuant.

## Outputs
- Tabulated summary of novel TRP isoforms with mean TPM
- List of novel TRP transcript IDs
- GTF file containing only novel TRP transcript models

Plots derived from these outputs are located in `downstream/tpm/figures/`.
