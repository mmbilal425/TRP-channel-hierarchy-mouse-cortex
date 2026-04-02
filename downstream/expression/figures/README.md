# Expression figures – TRP mouse cortex

This directory contains plotting scripts used to generate expression-based
main and supplementary figures for the TRP mouse cortex study except m6A.

---

## Main figures

### Figure 1 – Transcriptomic overview
- **Fig. 1B** – Nanopore-derived TRP TPMs in adult mouse cortex  
  (`Fig1B_nanocount_TRP_TPM_cortex_bar.py`)

- **Fig. 1C** – Illumina (Salmon) TRP TPMs in cortex (mean across replicates)  
  (`Fig1C_salmon_TRP_TPM_cortex_mean_r115.py`)

- **Fig. 1D** – Illumina (Salmon) TRP TPMs in DRG  
  (`Fig1D_salmon_TRP_TPM_DRG_r115.py`)

### Figure 2 – qPCR validation and cross-platform comparison
- **Fig. 2A** – Mean Ct values for TRP genes  
  (`Fig2A_mean_Ct_value_trp_genes.py`)

- **Fig. 2B** – Mean Ct bar plots for TRP genes  
  (`Fig2B_mean_Ct_value_bar_plots_trp_genes.py`)

- **Fig. 2C** – Correlation plots across transcriptomic platforms  
  (`Fig2C_correlation_plots.py`)

### Figure 4 – Targeted validation in bulk cortex and FACS samples
- **Fig. 4A** – Bulk cortex qPCR values for Trpa1 and Trpv1  
  (`Fig4A_ct_value_trpa1_trpv1_bulk_cortex.py`)

- **Fig. 4B** – FACS qPCR values for Trpa1 and Trpv1  
  (`Fig4B_facs_ct_value_trpa1_trpv1.py`)

---

## Supplementary Figure 2 – Read-count-based comparisons
- **Supp Fig. 2a** – Nanopore TRP read counts (cortex)  
  (`SuppFig2a_nanocount_TRP_readcounts_cortex_bar.py`)

- **Supp Fig. 2b** – Illumina TRP read counts (cortex mean)  
  (`SuppFig2b_salmon_TRP_readcounts_cortexMean_r115_panel.py`)

- **Supp Fig. 2c** – Illumina TRP read counts (DRG)  
  (`SuppFig2c_salmon_TRP_readcounts_DRG_r115_panel.py`)

- **Supp Fig. 2d** – Cortex vs DRG TRP expression (dumbbell plot, log₁₀ TPM)  
  (`SuppFig2d_salmon_TRP_TPM_cortex_vs_DRG_dumbbell_log10_r115_panel.py`)

---

## Supplementary Figure 3 – Expression distribution and contrast
- **Supp Fig. 3a** – TRP family TPM distribution (Nanopore)  
  (`SuppFig3a_nanocount_TRP_family_TPM_piechart_cortex.py`)

- **Supp Fig. 3b** – TRP family TPM distribution (Illumina)  
  (`SuppFig3b_salmon_TRP_family_TPM_pie_cortexMean_r115.py`)

- **Supp Fig. 3c** – Overlap of expressed TRP genes (Cortex vs DRG)  
  (`SuppFig3c_salmon_TRP_overlap_venn_TPMge0p5_r115.py`)

- **Supp Fig. 3d** – TRP log₂ fold-change heatmap (Cortex − DRG)  
  (`SuppFig3d_salmon_TRP_log2FC_heatmap_cortex_minus_drg_r115.py`)

---

## Supplementary Figure 4 – Isoform-level expression
- **Supp Fig. 4a** – Known TRP isoform TPMs (IsoQuant)  
  (`SuppFig4a_isoquant_known_TRP_isoforms_TPM_bar.py`)

- **Supp Fig. 4b** – Novel TRP isoform TPMs (IsoQuant)  
  (`SuppFig4b_isoquant_novel_TRP_isoforms_TPM_bar.py`)
