# Proteomics

This folder contains proteomics analysis code added separately from the transcriptomic workflows.

This section was added to avoid changing the existing transcriptomic analysis structure.

# Notes
Large raw proteomics data files are not included in this repository.
Input tables provided here are processed subsets used for figure generation.
Folder structure is designed to keep proteomics analyses modular and independent from transcriptomic workflows.

## Structure

### GO_analysis/
Gene Ontology (GO:CC) enrichment analysis and bubble plot generation.

- `code/`: scripts for GO enrichment visualization  
- `data/`: GO result tables (`*_GO*.csv`) used as input  
- `output/`: generated bubble plots and summary tables  

---

### IP_analysis/
TRP immunoprecipitation (IP) rerun analysis and peptide-level heatmap.

- `code/`: TRP IP rerun analysis script  
- `data/`: input Excel file containing peptide intensities  
- `output/`: heatmap figures (PDF/PNG)  

---

### heatmap/
TRP protein intensity heatmaps across cortex and peripheral organs.

- `code/`: heatmap generation scripts  
- `data/`: input tables for cortex and organ comparisons  
- `output/`: publication-ready heatmaps  

---

### proteomics_figures/
Additional proteomics figures used in the manuscript.

#### heatmap_log/
Log10-transformed TRP intensity heatmaps.

- `code/`: log-transformed heatmap script  
- `data/`: processed protein intensity table  
- `output/`: Figure 3E and Supplementary Figure S8D  

#### protein_per_group/
Protein counts per sample and PCA legend generation.

- `code/`: protein count and PCA legend scripts  
- `data/`: proteinGroups and experimental design tables  
- `output/`: Figure 3C, Supplementary Figure S8A, PCA legends  

---

## Reproducibility

All scripts are designed to run locally using relative paths.  
Each analysis folder contains:

- `code/`: executable scripts  
- `data/`: required input files  
- `output/`: generated figures  

To reproduce a figure, navigate to the corresponding folder and run the script, for example:

```bash
python code/TRP_rerun.py
