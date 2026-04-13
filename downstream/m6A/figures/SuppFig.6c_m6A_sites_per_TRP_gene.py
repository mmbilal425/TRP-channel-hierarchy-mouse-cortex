#!/usr/bin/env python3
"""
Supplementary Figure S5c
Number of modified m6A sites per TRP gene (Nanopore direct RNA; modkit)

This script:
  1) loads the gene-assigned m6A table (from upstream m6A pipeline)
  2) filters to TRP genes
  3) counts modified sites per TRP gene
  4) saves a TSV and a publication-ready PDF bar plot

Inputs
------
1) Gene-assigned m6A table (TSV or Parquet) containing at least:
   - gene_id
   - is_modified
   - mod_type   (optional but recommended)
   - merged_bam_coverage (optional)
   - merged_bam_stoich   (optional)

2) TRP gene ID list:
   - a whitespace-delimited file with columns: Geneid  gene_name
     (your file: trp_gene_ids.txt)

Outputs
-------
- modified_sites_per_TRP_gene.tsv
- SuppFig_S5c_m6A_sites_per_TRP_gene.pdf
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

logging.getLogger("fontTools").setLevel(logging.ERROR)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")


# -----------------------------
# Helpers
# -----------------------------
def read_table(path: Path) -> pd.DataFrame:
    """Read TSV/CSV/Parquet based on file extension."""
    suffix = path.suffix.lower()
    if suffix in [".parquet"]:
        return pd.read_parquet(path)
    if suffix in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t", low_memory=False)
    if suffix in [".csv"]:
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported input format: {path}")


def plot_sites_per_gene(count_df: pd.DataFrame, out_pdf: Path, y_max: int = 75) -> None:
    """Make publication-style bar plot."""
    if count_df.empty:
        raise ValueError("No TRP genes with modified sites found after filtering.")

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 6,
        "pdf.fonttype": 42,
        "xtick.major.pad": 1.4,
        "ytick.major.pad": 1.4,
    })

    fig = plt.figure(figsize=(3, 1))
    fig.subplots_adjust(left=0.18, right=0.95, bottom=0.30, top=0.92)
    ax = fig.add_subplot(111)

    x = np.arange(len(count_df))
    ax.bar(x, count_df["modified_sites"].values, color="blue", edgecolor="none")

    ax.set_xlabel("Gene", fontsize=6, labelpad=1.4)
    ax.set_ylabel("Number of modified sites", fontsize=6, labelpad=1.4)

    # Spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    # Ticks
    tick_length = 0.06 * 28.35  # 0.06 cm in points
    ax.tick_params(axis="both", which="major",
                   length=tick_length, width=0.5,
                   direction="out", pad=1.4)

    ax.set_xticks(x)
    ax.set_xticklabels(count_df["gene_name"].values, rotation=90, ha="center")
    ax.margins(x=0.01)
    ax.tick_params(axis="x", pad=2.5)

    # Y axis
    ax.set_yticks([0, 25, 50, 75] if y_max >= 75 else np.linspace(0, y_max, 4, dtype=int))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{int(v)}"))
    ax.set_ylim(0, y_max)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Main
# -----------------------------
def main():
    p = argparse.ArgumentParser(description="Supp Fig S5c: count modified m6A sites per TRP gene.")
    p.add_argument("--input", required=True, help="Gene-assigned m6A table (TSV/CSV/Parquet).")
    p.add_argument("--trp-genes", required=True, help="TRP gene list with columns: Geneid gene_name.")
    p.add_argument("--outdir", required=True, help="Output directory for TSV+PDF.")
    p.add_argument("--ymax", type=int, default=75, help="Y-axis max for the plot (default: 75).")
    p.add_argument("--require-modtype-a", action="store_true",
                   help="If set, require mod_type == 'a' (m6A) before counting.")
    args = p.parse_args()

    input_path = Path(args.input)
    trp_path = Path(args.trp_genes)
    outdir = Path(args.outdir)

    logging.info(f"Reading input table: {input_path}")
    df = read_table(input_path)

    # Validate columns
    required = {"gene_id", "is_modified"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    # Load TRP genes
    logging.info(f"Reading TRP gene list: {trp_path}")
    gene_table = pd.read_csv(trp_path, sep=r"\s+", engine="python")
    if "Geneid" not in gene_table.columns or "gene_name" not in gene_table.columns:
        raise ValueError("TRP gene list must contain columns: Geneid and gene_name")

    genes_of_interest = gene_table["Geneid"].tolist()
    id_to_name = dict(zip(gene_table["Geneid"], gene_table["gene_name"]))

    # Filter to modified + TRP genes (+ optionally m6A only)
    mask = (df["is_modified"] == 1) & (df["gene_id"].isin(genes_of_interest))
    if args.require_modtype_a:
        if "mod_type" not in df.columns:
            raise ValueError("--require-modtype-a set but input table has no 'mod_type' column.")
        mask = mask & (df["mod_type"] == "a")

    bar_df = df.loc[mask].copy()
    logging.info(f"Modified TRP-associated sites: {len(bar_df):,}")

    # Count sites per gene_id
    site_counts = (
        bar_df.groupby("gene_id", observed=True)
              .size()
              .sort_values(ascending=False)
              .rename_axis("gene_id")
              .reset_index(name="modified_sites")
    )
    site_counts["gene_name"] = site_counts["gene_id"].map(id_to_name)
    site_counts = site_counts[["gene_name", "modified_sites"]]

    # Save TSV
    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / "modified_sites_per_TRP_gene.tsv"
    site_counts.to_csv(tsv_path, sep="\t", index=False)
    logging.info(f"Counts table saved: {tsv_path}")

    # Plot
    pdf_path = outdir / "SuppFig_S5c_m6A_sites_per_TRP_gene.pdf"
    plot_sites_per_gene(site_counts, pdf_path, y_max=args.ymax)
    logging.info(f"Plot saved: {pdf_path}")


if __name__ == "__main__":
    main()
