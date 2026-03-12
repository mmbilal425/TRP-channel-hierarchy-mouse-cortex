## m6A Figures (Nanopore direct RNA; modkit)

This folder contains scripts used to generate the publication-ready panels for the
m6A epitranscriptomic analysis (Supplementary Fig. S6).

### Inputs (from upstream m6A pipeline)

These figure scripts expect the annotated m6A tables produced by the upstream steps in:

- `downstream/m6A/00_import_modkit_pileup_to_parquet.py`
- `downstream/m6A/01_attach_motifs.py`
- `downstream/m6A/02_call_modified_sites_and_export_bed.py`
- `downstream/m6A/03_assign_sites_to_genes.py`

In practice, the figure scripts require a table that contains at minimum:

- `mod_type` (m6A is `"a"`)
- `merged_bam_coverage`
- `merged_bam_stoich`
- `is_modified`

Additional required fields:
- for **S6b**: `motif`, `is_drach`, `is_rac`
- for **S6c**: `gene_id` (plus a TRP gene list mapping `Geneid -> gene_name`)

### Supplementary Figure S6

| Panel | Script | Description | Outputs |
|------|--------|-------------|---------|
| S6a | `SuppFig_S6a_m6A_stoichiometry_TRP.py` | Histogram of m6A stoichiometry (%) for modified m6A sites associated with TRP genes | `stoichiometry_histogram.tsv`, `stoichiometry_histogram.pdf` |
| S6b | `SuppFig_S6b_m6A_ranked_motifs_TRP.py` | Ranked motif analysis of modified m6A sites; highlights canonical DRACH, RAC non-DRACH, and other motifs | `motif_modification_rates_m6A.tsv`, `motif_modification_rates_m6A.pdf`, `motif_modification_rates_m6A_legend.pdf` |
| S6c | `SuppFig_S6c_m6A_sites_per_TRP_gene.py` | Number of modified m6A sites per TRP gene (bar plot) | `modified_sites_per_TRP_gene.tsv`, `SuppFig_S6c_m6A_sites_per_TRP_gene.pdf` |

### Notes

- Figures are exported as vector PDF suitable for publication.
- The `.py` scripts in this folder are the reproducible source of truth for figure generation.
- If final PDFs are added to the repository for quick viewing, commit only small final outputs rather than intermediate files.