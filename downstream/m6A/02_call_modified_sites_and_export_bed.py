#!/usr/bin/env python3
"""
02_call_modified_sites_and_export_bed.py

Calls "is_modified" using your thresholds:
- coverage >= 2
- stoich >= 5
- modified_reads >= 2   where modified_reads = round(coverage * stoich/100)

Exports:
- is_modified_sites.bed  (chr, start, start+1)
- modified_sites_table.tsv (optional)
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data"
PLOTS_DIR = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/plots"

IN_PARQUET = "merged_m6a_pileups_with_motifs.parquet"

MIN_STOICH = 5
MIN_MOD_READS = 2
MIN_COV = 2

OUT_BED = "is_modified_sites.bed"
OUT_TSV = "modified_sites_table.tsv"
OUT_PARQUET = "modified_sites_table.parquet"


def main():
    data_dir = Path(DATA_DIR)
    plots_dir = Path(PLOTS_DIR)
    plots_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(data_dir / IN_PARQUET)

    # enforce coverage filter
    df = df[df["merged_bam_coverage"] >= MIN_COV].copy()

    # compute modified reads
    df["total_modified_reads"] = np.round(
        df["merged_bam_coverage"] * df["merged_bam_stoich"] / 100.0
    ).astype(int)

    df["is_modified"] = (
        (df["merged_bam_stoich"] >= MIN_STOICH) &
        (df["total_modified_reads"] >= MIN_MOD_READS)
    ).astype(np.int8)

    modified_df = df[df["is_modified"] == 1].copy()

    # save BED
    bed_df = modified_df[["chr", "start"]].copy()
    bed_df["end"] = bed_df["start"] + 1
    bed_df = bed_df.sort_values(["chr", "start"])

    bed_path = plots_dir / OUT_BED
    bed_df.to_csv(bed_path, sep="\t", header=False, index=False)

    # save table (optional)
    tsv_path = plots_dir / OUT_TSV
    modified_df.to_csv(tsv_path, sep="\t", index=False)

    pq_path = plots_dir / OUT_PARQUET
    modified_df.to_parquet(pq_path, index=False, compression="snappy")

    print(f"Modified sites: {len(modified_df):,}")
    print(f"Saved BED: {bed_path}")
    print(f"Saved TSV: {tsv_path}")
    print(f"Saved parquet: {pq_path}")


if __name__ == "__main__":
    main()
