# m6A downstream analysis (Nanopore direct RNA; modkit)

## Input
- modkit pileup BED:
  /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/merged_bam_files/m6a_modkit_output/mouse_cortex_merged_pileup_filt0.5_mod0.99_global0.9_allMods.bed
- reference genome:
  /g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa
- GTF (Ensembl r115):
  /g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf
- gene intervals (for assignment):
  /g/data/qq78/as7425/rna_biogenesis_maps/dFORCE_prod_dec24/mouse/second_pass_index/annotation_with_clusters/updated_gene.bed

## Outputs (HPC paths)
- data:
  /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data
- plots:
  /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/plots

## Run order
1) 00_import_modkit_pileup_to_parquet.py
   -> merged_m6a_pileups.parquet

2) 01_attach_motifs.py
   -> merged_m6a_pileups_with_motifs.parquet (or tsv/parquet output)

3) 02_call_modified_sites_and_export_bed.py
   -> is_modified_sites.bed (+ summary tables)

4) 03_assign_sites_to_genes.py
   -> m6a_sites_with_gene_id.tsv/parquet

5) 04_fig_motif_modification_rate.py
   -> motif_modification_rates_m6A.pdf (+ legend + tsv)

6) 05_fig_stoichiometry_histogram.py
   -> stoichiometry_histogram.pdf (+ tsv)

7) 06_fig_trp_modified_sites_bar.py
   -> modified_sites_per_gene.pdf (+ tsv)
