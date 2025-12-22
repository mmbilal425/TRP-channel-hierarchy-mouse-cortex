#!/usr/bin/env python3
"""
Supplementary Fig. S2a — Nanopore (NanoCount) Cortex TRP family TPM pie chart

This script:
  1) Reads the ordered TRP TPM table produced by Fig1B pipeline:
       TRP_TPM_CortexMean_NanoCount_ordered.tsv
  2) Aggregates TPM_mean by TRP family
  3) Saves a family TPM summary table
  4) Generates a TPM-weighted pie chart (publication PDF)

Inputs:
  --tpm_file: TRP_TPM_CortexMean_NanoCount_ordered.tsv
Outputs:
  TRP_TPM_CortexMean_NanoCount_family_TPM_piechart.pdf
  TRP_TPM_CortexMean_NanoCount_family_TPM_totals.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm


def pick_font() -> str:
    available = {f.name for f in fm.fontManager.ttflist}
    return "Arial" if "Arial" in available else "DejaVu Sans"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Supp Fig S2a: Nanopore (NanoCount) Cortex TRP family TPM pie chart."
    )
    parser.add_argument(
        "--tpm_file",
        type=Path,
        required=True,
        help="Ordered TPM TSV (TRP_TPM_CortexMean_NanoCount_ordered.tsv).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for PDF + TSV.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Export DPI for PDF (default: 600).",
    )
    args = parser.parse_args()

    font_choice = pick_font()
    plt.rcParams.update(
        {
            "font.family": font_choice,
            "font.size": 14,
            "pdf.fonttype": 42,
            "text.color": "black",
            "axes.labelcolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "axes.edgecolor": "black",
        }
    )

    tpm_file: Path = args.tpm_file
    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    out_pdf = outdir / "TRP_TPM_CortexMean_NanoCount_family_TPM_piechart.pdf"
    out_tsv = outdir / "TRP_TPM_CortexMean_NanoCount_family_TPM_totals.tsv"

    # ---------- load ----------
    df = pd.read_csv(tpm_file, sep="\t")

    required = {"GeneName", "family", "TPM_mean"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in {tpm_file}: {sorted(missing)}\n"
            f"Found: {list(df.columns)}"
        )

    df = df[["GeneName", "family", "TPM_mean"]].copy()
    df["GeneName"] = df["GeneName"].astype(str)
    df["family"] = df["family"].astype(str)
    df["TPM_mean"] = pd.to_numeric(df["TPM_mean"], errors="coerce").fillna(0.0)

    # ---------- aggregate ----------
    family_tpm = (
        df.groupby("family", as_index=False)["TPM_mean"]
        .sum()
        .rename(columns={"TPM_mean": "TPM"})
        .sort_values("TPM", ascending=False)
    )

    family_tpm.to_csv(out_tsv, sep="\t", index=False)

    # ---------- palette (consistent with your TRP family scheme) ----------
    palette = {
        "TRPML": "#4682b4",  # Mcoln
        "TRPP": "#dda0dd",   # Pkd
        "TRPA": "#696969",
        "TRPC": "#fc8d62",
        "TRPM": "#66a61e",
        "TRPV": "#8da0cb",
    }
    colors = [palette.get(f, "#999999") for f in family_tpm["family"]]

    # ---------- plot ----------
    plt.figure(figsize=(7, 7))
    wedges, texts, autotexts = plt.pie(
        family_tpm["TPM"].values,
        labels=family_tpm["family"].values,
        autopct="%1.1f%%",
        startangle=90,
        counterclock=False,
        colors=colors,
        textprops={"fontsize": 12},
    )

    # If TRPA slice is tiny, nudge percent text slightly upward
    fam_list = list(family_tpm["family"])
    if "TRPA" in fam_list:
        idx = fam_list.index("TRPA")
        x, y = autotexts[idx].get_position()
        autotexts[idx].set_position((x, y + 0.12))

    plt.title("TPM distribution by TRP gene family (Nanopore, cortex)", fontsize=18)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=args.dpi, format="pdf", bbox_inches="tight")
    plt.close()

    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_tsv}")


if __name__ == "__main__":
    main()
