#!/usr/bin/env python3
"""
00_filter_novel_TRP_isoforms.py

Filter IsoQuant SQANTI-like classification to TRP genes and novel categories (NIC/NNC).

Inputs:
  - SQANTI-like TSV (no header): <prefix>.novel_vs_known.SQANTI-like.tsv
  - TRP gene list: trp_gene_ids.txt (expects columns including 'Geneid' and 'gene_name')

Outputs:
  - TSV of novel TRP isoforms (classification only; no TPMs yet)
"""

from pathlib import Path
import pandas as pd


# ---------------- USER CONFIG ----------------
BASE = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/Mouse_Cortex_dRNA")
TRP_LIST = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")

SQANTI_FILE = BASE / "Mouse_Cortex_dRNA.novel_vs_known.SQANTI-like.tsv"

OUTDIR = BASE / "TRP_novel"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_FILTERED = OUTDIR / "novel_TRP_isoforms_classification.tsv"

NOVEL_SET = {"novel_in_catalog", "novel_not_in_catalog"}  # NIC/NNC
# ------------------------------------------------


def main() -> None:
    if not SQANTI_FILE.exists():
        raise FileNotFoundError(f"Missing SQANTI-like file: {SQANTI_FILE}")
    if not TRP_LIST.exists():
        raise FileNotFoundError(f"Missing TRP list file: {TRP_LIST}")

    trp = pd.read_csv(TRP_LIST, sep="\t", dtype=str)
    trp.columns = [c.strip() for c in trp.columns]

    if "Geneid" not in trp.columns:
        raise ValueError(f"TRP list must contain 'Geneid' column. Found: {list(trp.columns)}")

    trp_gene_ids = set(trp["Geneid"].dropna().astype(str))

    # SQANTI-like is usually headerless; we only need a few fields.
    # Based on your earlier working setup:
    sq_cols = [
        "transcript_id", "exons", "strand", "length", "TPM",
        "structural_category", "gene_id", "ref_transcript_id", "exon_count", "intron_count",
        "junction_count"
    ]
    sq = pd.read_csv(SQANTI_FILE, sep="\t", header=None, names=None, low_memory=False)

    # Defensive: ensure at least 7 columns exist (for transcript_id, structural_category, gene_id)
    if sq.shape[1] < 7:
        raise ValueError(f"SQANTI-like file has too few columns ({sq.shape[1]}). Check file: {SQANTI_FILE}")

    # Map required columns by position (as in your earlier workflow):
    # 0 transcript_id, 5 structural_category, 6 gene_id, and later isoform_type (often around col 27)
    sq_df = pd.DataFrame({
        "transcript_id": sq.iloc[:, 0].astype(str),
        "structural_category": sq.iloc[:, 5].astype(str),
        "gene_id": sq.iloc[:, 6].astype(str),
    })

    # isoform_type position can vary; try to pull if present (your earlier used ~27)
    if sq.shape[1] > 27:
        sq_df["isoform_type"] = sq.iloc[:, 27].astype(str)
    else:
        sq_df["isoform_type"] = ""

    sq_df = sq_df.dropna(subset=["transcript_id", "gene_id"])
    before = len(sq_df)

    sq_df = sq_df[
        sq_df["structural_category"].isin(NOVEL_SET) &
        sq_df["gene_id"].isin(trp_gene_ids)
    ].copy()

    print(f"[00] SQANTI rows: {before:,}")
    print(f"[00] Novel TRP isoforms (NIC/NNC): {len(sq_df):,}")

    sq_df.to_csv(OUT_FILTERED, sep="\t", index=False)
    print(f"[00] Wrote: {OUT_FILTERED}")


if __name__ == "__main__":
    main()

