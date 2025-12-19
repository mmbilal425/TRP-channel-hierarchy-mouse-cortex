## m6A Figures (Nanopore direct RNA; modkit)

This folder contains scripts used to generate the publication-ready panels for the
m6A epitranscriptomic analysis (Supplementary Fig. S5).

### Inputs (from upstream m6A pipeline)

These figure scripts expect the annotated m6A tables produced by the upstream steps in:

- `downstream/m6A/00_import_modkit_pileup_to_parquet.py`
- `downstream/m6A/01_attach_motifs.py`
- `downstream/m6A/02_call_modified_sites_and_export_bed.py`
- `downstream/m6A/03_assign_sites_to_genes.py`

In practice, the figure scripts require a table that contains (at minimum):
- `mod_type` (m6A is `"a"`)
- `merged_bam_coverage`
- `merged_bam_stoich`
- `is_modified`
and for S5b:
- `motif`, `is_drach`, `is_rac`
and for S5c:
- `gene_id` (and a TRP gene list mapping `Geneid -> gene_name`)

### Supplementary Figure S5

| Panel | Script | Description | Outputs |
|------|--------|-------------|---------|
| S5a | `SuppFig_S5a_m6A_stoichiometry_TRP.py` | Histogram of m6A stoichiometry (%) for modified m6A sites (TRP-associated) | `stoichiometry_histogram.tsv`, `stoichiometry_histogram.pdf` |
| S5b | `SuppFig_S5b_m6A_ranked_motifs_TRP.py` | Ranked motif analysis of modified m6A sites; highlights canonical DRACH (red), RAC non-DRACH (blue), other (grey) | `motif_modification_rates_m6A.tsv`, `motif_modification_rates_m6A.pdf`, `motif_modification_rates_m6A_legend.pdf` |
| S5c | `SuppFig_S5c_m6A_sites_per_TRP_gene.py` | Count of modified m6A sites per TRP gene (bar plot) | `modified_sites_per_TRP_gene.tsv`, `SuppFig_S5c_m6A_sites_per_TRP_gene.pdf` |

### Notes

- Figures are exported as vector PDF (300 dpi) suitable for publication.
- If you previously uploaded `.ipynb` files and GitHub shows **Invalid Notebook**, those files are not valid JSON.
  Prefer the `.py` scripts in this folder as the reproducible source of truth.
- If you want the plots written into the repository for quick viewing, save the PDFs locally and commit them
  (recommended only for small final PDFs, not intermediate large files).
