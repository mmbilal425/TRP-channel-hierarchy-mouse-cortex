# m6A downstream analysis (Nanopore direct RNA; modkit)

This folder contains the downstream processing steps used to quantify m6A sites from a **modkit pileup BED** generated on merged Nanopore direct RNA reads, annotate sequence motifs (DRACH/RAC), call modified sites, assign sites to genes, and generate **Supplementary Fig. S6 (a–c)**.

---

## Inputs (HPC paths)

- **modkit pileup BED (merged BAM):**  
  `/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/m6a_modkit_output/mouse_cortex_merged_pileup_filt0.5_mod0.99_global0.9_allMods.bed`

- **Reference genome (GRCm39 primary assembly FASTA):**  
  `/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa`

- **GTF annotation (Ensembl r115):**  
  `/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf`

- **Gene intervals for assignment (BED):**  
  `/g/data/qq78/as7425/rna_biogenesis_maps/dFORCE_prod_dec24/mouse/second_pass_index/annotation_with_clusters/updated_gene.bed`

---

## Outputs (HPC paths)

- **Intermediate tables:**  
  `/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data`

- **Plots:**  
  `/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/plots`

---

## Run order (processing)

### 1) Import modkit pileup BED → compact parquet
**Script:** `00_import_modkit_pileup_to_parquet.py`  
**Output (example):**
- `merged_m6a_pileups.parquet`  
  (chr, start, strand, coverage, stoichiometry for mod_type = "a", coverage ≥ 2)

### 2) Attach 5-mer motifs + DRACH/RAC classification
**Script:** `01_attach_motifs.py`  
**Output (example):**
- `merged_m6a_pileups_with_motifs.parquet`  
  (adds: motif, is_drach, is_rac)

### 3) Call modified sites + export BED
**Script:** `02_call_modified_sites_and_export_bed.py`  
**Output (example):**
- `is_modified_sites.bed`
- summary tables (TSV/Parquet)

### 4) Assign sites to genes by overlap
**Script:** `03_assign_sites_to_genes.py`  
**Output (example):**
- `m6a_sites_with_gene_id.tsv` (or `.parquet`)

---

## Figures (Supplementary Fig. S6)

Figure scripts are stored in: **`figures/`**

- **S6a:** `figures/SuppFig_S6a_m6A_stoichiometry_TRP.py`  
  → `stoichiometry_histogram.pdf` (+ `stoichiometry_histogram.tsv`)

- **S6b:** `figures/SuppFig_S6b_m6A_ranked_motifs_TRP.py`  
  → `motif_modification_rates_m6A.pdf`  
  → `motif_modification_rates_m6A_legend.pdf`  
  → `motif_modification_rates_m6A.tsv`

- **S6c:** `figures/SuppFig_S6c_m6A_sites_per_TRP_gene.py`  
  → `modified_sites_per_gene.pdf` (+ `modified_sites_per_gene.tsv`)

---

## Notes

- These scripts were developed and run on the NCI Gadi environment (HPC paths shown above).
- Output figure PDFs are exported at **publication-ready settings** (vector PDF; 300 dpi where relevant).
- For reproducibility, each step writes intermediate tables that can be reused without rerunning upstream steps.
