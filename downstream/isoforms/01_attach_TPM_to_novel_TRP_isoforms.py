#!/usr/bin/env python3
"""
01_attach_TPM_to_novel_TRP_isoforms.py

Attach replicate TPMs from IsoQuant discovered transcript TPM table to the filtered novel TRP isoforms,
and compute mean TPM.

Inputs:
  - novel_TRP_isoforms_classification.tsv (from step 00)
  - <prefix>.discovered_transcript_grouped_tpm.tsv

Outputs:
  - novel_TRP_isoforms_with_TPM.tsv
"""

from pathlib import Path
import pandas as pd
import numpy as np


# ---------------- USER CONFIG ----------------
BASE = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/Mouse_Cortex_dRNA")

TPM_FILE = BASE / "Mouse_Cortex_dRNA.discovered_transcript_grouped_tpm.tsv"
IN_FILTERED = BASE / "TRP_novel/novel_TRP_isoforms_classification.tsv"

OUTDIR = BASE / "TRP_novel"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_WITH_TPM = OUTDIR / "novel_TRP_isoforms_with_TPM.tsv"
# ------------------------------------------------


def main() -> None:
    if not IN_FILTERED.exists():
        raise FileNotFoundError(f"Missing input: {IN_FILTERED}")
    if not TPM_FILE.exists():
        raise FileNotFoundError(f"Missing TPM file: {TPM_FILE}")

    df = pd.read_csv(IN_FILTERED, sep="\t", dtype=str)
    if "transcript_id" not in df.columns:
        raise ValueError(f"Expected 'transcript_id' in {IN_FILTERED}. Found: {list(df.columns)}")

    tpm = pd.read_csv(TPM_FILE, sep="\t")
    # The ID column may be "#feature_id" or "feature_id" depending on IsoQuant version
    if "#feature_id" in tpm.columns:
        tpm = tpm.rename(columns={"#feature_id": "transcript_id"})
    elif "feature_id" in tpm.columns:
        tpm = tpm.rename(columns={"feature_id": "transcript_id"})
    else:
        # fallback: first column
        tpm = tpm.rename(columns={tpm.columns[0]: "transcript_id"})

    # Keep only numeric TPM columns (replicate columns)
    rep_cols = [c for c in tpm.columns if c != "transcript_id"]
    rep_cols = [c for c in rep_cols if np.issubdtype(tpm[c].dtype, np.number)]

    if len(rep_cols) == 0:
        raise ValueError(f"No numeric TPM replicate columns found in {TPM_FILE}. Columns: {list(tpm.columns)}")

    merged = df.merge(tpm[["transcript_id"] + rep_cols], on="transcript_id", how="left")

    for c in rep_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0.0)

    merged["mean_TPM"] = merged[rep_cols].mean(axis=1, skipna=True)

    # Order for readability
    keep_cols = ["gene_id", "transcript_id", "structural_category", "isoform_type", "mean_TPM"] + rep_cols
    for c in keep_cols:
        if c not in merged.columns:
            merged[c] = "" if c in ["gene_id", "structural_category", "isoform_type"] else 0.0

    out = merged[keep_cols].copy()
    out = out.sort_values(["gene_id", "structural_category", "mean_TPM"], ascending=[True, True, False])

    out.to_csv(OUT_WITH_TPM, sep="\t", index=False)
    print(f"[01] Wrote: {OUT_WITH_TPM}")
    print(f"[01] Rows: {len(out):,}")


if __name__ == "__main__":
    main()

