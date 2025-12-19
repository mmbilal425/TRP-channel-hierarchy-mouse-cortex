#!/usr/bin/env python3
"""
Supplementary Figure S5b
Ranked motif analysis for modified m6A sites (Nanopore direct RNA; modkit)

What this script does
---------------------
1) Loads an annotated m6A table (must contain motifs + DRACH/RAC flags).
2) Recomputes "is_modified" using:
      - min_stoichiometry (%)
      - min_modified_reads (coverage * stoich/100)
   and filters to coverage >= 2 and mod_type == "a"
3) Aggregates by 5-mer motif:
      - modified_sites (sum is_modified)
      - total_sites
      - modification_rate = modified_sites / total_sites
4) Plots top-N motifs ranked by modification_rate, colored by:
      - DRACH (red)
      - RAC but not DRACH (blue)
      - other (lavender)
   and saves:
      - motif_modification_rates_m6A.pdf
      - motif_modification_rates_m6A_legend.pdf
      - motif_modification_rates_m6A.tsv

Inputs (expected columns)
-------------------------
Required columns in the input table:
  - motif (5-mer string)
  - is_drach (0/1)
  - is_rac (0/1)
  - mod_type (e.g. "a" for m6A)
  - merged_bam_coverage (int)
  - merged_bam_stoich (% float)
  - total_protein_coding_modified_reads (int or float; if missing, we compute it)

Example upstream provenance:
  - 00_import_modkit_pileup_to_parquet.py
  - 01_attach_motifs.py
  - 02_call_modified_sites_and_export_bed.py (often where modified_reads / is_modified are computed)

Author: TRP_mouse_cortex pipeline
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Silence fontTools chatter in some environments
logging.getLogger("fontTools").setLevel(logging.ERROR)

# -----------------------------
# Plot styling (Nature-ish)
# -----------------------------
def set_plot_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 6,
        "pdf.fonttype": 42,
        "xtick.major.pad": 1.4,
        "ytick.major.pad": 1.4,
    })


def parse_args():
    p = argparse.ArgumentParser(
        description="Supp Fig S5b: Ranked m6A motifs (DRACH/RAC enrichment)."
    )
    p.add_argument(
        "--input",
        default="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/m6a_sites_annotated.parquet",
        help="Input table (parquet/tsv/csv) containing motif + is_drach/is_rac + coverage/stoich."
    )
    p.add_argument(
        "--outdir",
        default="/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/plots",
        help="Output directory for plots + TSV."
    )
    p.add_argument("--min-stoich", type=float, default=5.0, help="Min stoichiometry (%) to call modified.")
    p.add_argument("--min-mod-reads", type=int, default=2, help="Min modified reads to call modified.")
    p.add_argument("--min-coverage", type=int, default=2, help="Min coverage to keep a site.")
    p.add_argument("--top-n", type=int, default=40, help="Number of top motifs to plot.")
    return p.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported input format: {suffix}. Use parquet/tsv/csv.")


def require_columns(df: pd.DataFrame, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def main():
    args = parse_args()
    in_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    logging.info(f"Loading: {in_path}")
    df = load_table(in_path)

    # Required inputs for this figure
    require_columns(df, [
        "motif", "is_drach", "is_rac",
        "mod_type", "merged_bam_coverage", "merged_bam_stoich"
    ])

    # Compute modified reads if missing (your pipeline sometimes precomputes this)
    if "total_protein_coding_modified_reads" not in df.columns:
        logging.info("Column total_protein_coding_modified_reads missing; computing from coverage * stoich / 100.")
        df["total_protein_coding_modified_reads"] = np.round(
            df["merged_bam_coverage"].astype(float) * df["merged_bam_stoich"].astype(float) / 100.0
        ).astype(int)

    # Recompute is_modified using your thresholds
    df = df.copy()
    df["is_modified"] = (
        (df["merged_bam_stoich"] >= args.min_stoich) &
        (df["total_protein_coding_modified_reads"] >= args.min_mod_reads)
    ).astype(int)

    # Filters: coverage >= 2 and mod_type == 'a'
    df_f = df[
        (df["merged_bam_coverage"] >= args.min_coverage) &
        (df["mod_type"] == "a")
    ].copy()

    logging.info(f"After filters: {len(df_f):,} sites")

    # Aggregate by motif
    motif_stats = (
        df_f.groupby("motif", as_index=False)
            .agg(
                modified_sites=("is_modified", "sum"),
                total_sites=("is_modified", "size"),
                is_drach=("is_drach", "first"),
                is_rac=("is_rac", "first"),
            )
    )
    motif_stats["modification_rate"] = motif_stats["modified_sites"] / motif_stats["total_sites"]

    # Take top-N by modification rate
    motif_stats = motif_stats.sort_values("modification_rate", ascending=False).head(args.top_n).reset_index(drop=True)

    # Save table
    tsv_out = outdir / "motif_modification_rates_m6A.tsv"
    motif_stats.to_csv(tsv_out, sep="\t", index=False)
    logging.info(f"Wrote: {tsv_out}")

    # -----------------------------
    # Plot
    # -----------------------------
    set_plot_style()
    fig = plt.figure(figsize=(2.3 / 2.54, 2.7 / 2.54))
    fig.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.92)
    ax = fig.add_subplot(111)

    # Colors: DRACH red, RAC(non-DRACH) blue, Other lavender
    def motif_color(row):
        if int(row["is_drach"]) == 1:
            return "red"
        if int(row["is_rac"]) == 1 and int(row["is_drach"]) == 0:
            return "blue"
        return "lavender"

    ax.bar(
        np.arange(len(motif_stats)),
        motif_stats["modification_rate"].values,
        color=[motif_color(r) for _, r in motif_stats.iterrows()],
        width=1.0,
        align="edge",
        edgecolor="none",
    )

    # Median dashed line
    median = float(np.median(motif_stats["modification_rate"].values))
    ax.axhline(median, color="black", linestyle="--", linewidth=0.5)

    # Label top-5 motifs with curved connectors (same idea as your notebook)
    top_5 = motif_stats.head(5)
    max_height = float(motif_stats["modification_rate"].max())
    label_x_positions = np.linspace(8, max(10, args.top_n - 5), 5)
    label_y_positions = np.linspace(max_height * 1.05, max_height * 0.7, 5)

    for i, row in top_5.iterrows():
        x = i + 0.5
        y = float(row["modification_rate"])
        ax.annotate(
            f"{row['motif']}",
            xy=(x, y),
            xytext=(label_x_positions[i], label_y_positions[i]),
            ha="left",
            va="center",
            fontsize=4.5,
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                lw=0.5,
                connectionstyle="arc3,rad=+0.2",
            ),
        )

    ax.set_xlabel("Ranked motifs", fontsize=6, labelpad=1.4)
    ax.set_ylabel("m6A rate in mRNA", fontsize=6, labelpad=1.4)

    # Ticks and spines
    tick_length = 0.06 * 28.35
    ax.tick_params(axis="both", which="major", length=tick_length, width=0.5, direction="out", pad=1.4)

    ax.set_xticks([0, 20, 40])
    ax.set_xticklabels(["0", "20", "40"], fontsize=6)

    ax.set_yticks([0, 0.25, 0.50])
    ax.set_yticklabels(["0", "0.25", "0.5"], fontsize=6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    ax.set_xlim(-3, len(motif_stats))
    ax.set_ylim(0, 0.5)

    pdf_out = outdir / "motif_modification_rates_m6A.pdf"
    fig.savefig(pdf_out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Wrote: {pdf_out}")

    # -----------------------------
    # Legend figure (separate)
    # -----------------------------
    legend_fig = plt.figure(figsize=(1.5 / 2.54, 1.0 / 2.54))
    legend_ax = legend_fig.add_subplot(111)

    legend_elements = [
        Patch(facecolor="red", edgecolor="none", label="DRACH"),
        Patch(facecolor="blue", edgecolor="none", label="RAC (non-DRACH)"),
        Patch(facecolor="lavender", edgecolor="none", label="Other"),
    ]
    legend_ax.legend(
        handles=legend_elements,
        frameon=False,
        handletextpad=0.5,
        fontsize=6,
        loc="center",
        bbox_to_anchor=(0.5, 0.5),
    )
    legend_ax.set_axis_off()

    legend_out = outdir / "motif_modification_rates_m6A_legend.pdf"
    legend_fig.savefig(legend_out, dpi=300, bbox_inches="tight")
    plt.close(legend_fig)
    logging.info(f"Wrote: {legend_out}")

    # -----------------------------
    # Quick summary stats (stdout)
    # -----------------------------
    total_m6a = int(df_f["is_modified"].sum())
    drach_m6a = int(df_f[(df_f["is_drach"] == 1) & (df_f["is_modified"] == 1)].shape[0])
    rac_m6a = int(df_f[(df_f["is_rac"] == 1) & (df_f["is_drach"] == 0) & (df_f["is_modified"] == 1)].shape[0])
    other_m6a = int(df_f[(df_f["is_rac"] == 0) & (df_f["is_drach"] == 0) & (df_f["is_modified"] == 1)].shape[0])

    logging.info(f"Highest modification rate (top motifs table): {motif_stats['modification_rate'].max():.6f}")
    logging.info("m6A site distribution:")
    if total_m6a > 0:
        logging.info(f"  DRACH: {drach_m6a:,} ({drach_m6a/total_m6a*100:.1f}%)")
        logging.info(f"  RAC (non-DRACH): {rac_m6a:,} ({rac_m6a/total_m6a*100:.1f}%)")
        logging.info(f"  Other: {other_m6a:,} ({other_m6a/total_m6a*100:.1f}%)")
        logging.info(f"  Total modified m6A sites: {total_m6a:,}")
    else:
        logging.info("  No modified sites detected under the current thresholds.")


if __name__ == "__main__":
    main()
