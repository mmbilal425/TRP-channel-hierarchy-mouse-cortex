# m6A downstream analysis (Nanopore direct RNA; modkit)

This folder contains the downstream processing pipeline used to quantify m6A sites from a **modkit pileup BED** generated on merged Nanopore direct RNA reads, annotate sequence motifs (DRACH/RAC), call modified sites, assign sites to genes, and generate **Supplementary Fig. 6 (a–c)**.

Pipeline overview:
modkit pileup → motif annotation → modified site calling → gene assignment → figure generation

---

## Inputs

The following inputs are required (paths should be defined locally or via config files):

- **modkit pileup BED (merged BAM)**
- **Reference genome (GRCm39 primary assembly FASTA)**
- **GTF annotation (Ensembl r115)**
- **Gene interval BED file (for gene assignment)**

Example HPC paths used during analysis:
- `/g/data/lf10/.../m6a_modkit_output/...bed`
- `/g/data/lf10/.../reference_genomes/...fa`
- `/g/data/lf10/.../reference_genomes/...gtf`
- `/g/data/.../annotation_with_clusters/updated_gene.bed`

---

## Outputs

- Intermediate processed tables (Parquet/TSV)
- Publication-ready plots (PDF)

Example output directories:
- Tables: `/data/.../nanopore_data/.../data`
- Figures: `/data/.../nanopore_data/.../plots`

---

## Run order (processing)

### 1) Import modkit pileup BED → compact parquet
**Script:** `00_import_modkit_pileup_to_parquet.py`  
**Output:**
- `merged_m6a_pileups.parquet`  
  (chr, start, strand, coverage, stoichiometry; filtered for mod_type = "a")

---

### 2) Attach 5-mer motifs + DRACH/RAC classification
**Script:** `01_attach_motifs.py`  
**Output:**
- `merged_m6a_pileups_with_motifs.parquet`  
  (adds: motif, is_drach, is_rac)

---

### 3) Call modified sites + export BED
**Script:** `02_call_modified_sites_and_export_bed.py`  
**Output:**
- `is_modified_sites.bed`
- summary tables (TSV/Parquet)

---

### 4) Assign sites to genes by overlap
**Script:** `03_assign_sites_to_genes.py`  
**Output:**
- `m6a_sites_with_gene_id.tsv` (or `.parquet`)

---

## Figures (Supplementary Fig. 6)

Figure scripts are located in: **`figures/`**

These panels focus on **TRP genes only**.

- **6a — m6A stoichiometry distribution**
  - Script: `SuppFig_6a_m6A_stoichiometry_TRP.py`
  - Outputs:
    - `stoichiometry_histogram.pdf`
    - `stoichiometry_histogram.tsv`

- **6b — Ranked motif enrichment**
  - Script: `SuppFig_6b_m6A_ranked_motifs_TRP.py`
  - Outputs:
    - `motif_modification_rates_m6A.pdf`
    - `motif_modification_rates_m6A_legend.pdf`
    - `motif_modification_rates_m6A.tsv`

- **6c — m6A sites per TRP gene**
  - Script: `SuppFig_6c_m6A_sites_per_TRP_gene.py`
  - Outputs:
    - `modified_sites_per_TRP_gene.pdf`
    - `modified_sites_per_TRP_gene.tsv`

---

## Notes

- Analyses were performed on the **NCI Gadi HPC system**; example paths are provided for reference.
- All scripts are designed to be reproducible with user-defined paths.
- Figures are exported as **publication-ready vector PDFs** (300 dpi where applicable).
- Intermediate tables are saved to allow reuse without rerunning upstream steps.
