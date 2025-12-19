#!/usr/bin/env python3
"""
02_export_novel_TRP_GTF_and_IDs.py

Export transcript IDs and subset IsoQuant transcript_models.gtf for IGV validation.

Inputs:
  - novel_TRP_isoforms_with_TPM.tsv (from step 01)
  - <prefix>.transcript_models.gtf (IsoQuant output)

Outputs:
  - novel_TRP_transcript_ids.txt
  - novel_TRP_transcript_models.gtf
"""

from pathlib import Path
import pandas as pd
import re


# ---------------- USER CONFIG ----------------
BASE = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/Mouse_Cortex_dRNA")

IN_TSV = BASE / "TRP_novel/novel_TRP_isoforms_with_TPM.tsv"
MODELS_GTF = BASE / "Mouse_Cortex_dRNA.transcript_models.gtf"

OUTDIR = BASE / "TRP_novel"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_IDS = OUTDIR / "novel_TRP_transcript_ids.txt"
OUT_GTF = OUTDIR / "novel_TRP_transcript_models.gtf"
# ------------------------------------------------


def main() -> None:
    if not IN_TSV.exists():
        raise FileNotFoundError(f"Missing input TSV: {IN_TSV}")
    if not MODELS_GTF.exists():
        raise FileNotFoundError(f"Missing transcript models GTF: {MODELS_GTF}")

    df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
    if "transcript_id" not in df.columns:
        raise ValueError(f"Expected 'transcript_id' in {IN_TSV}. Found: {list(df.columns)}")

    ids = df["transcript_id"].dropna().astype(str).unique().tolist()
    OUT_IDS.write_text("\n".join(ids) + "\n")
    print(f"[02] Wrote IDs ({len(ids)}): {OUT_IDS}")

    tid_re = re.compile(r'transcript_id "([^"]+)"')
    kept = 0

    with open(MODELS_GTF, "r") as fin, open(OUT_GTF, "w") as fout:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            m = tid_re.search(line)
            if m and m.group(1) in set(ids):
                fout.write(line)
                kept += 1

    print(f"[02] Wrote subset GTF: {OUT_GTF} (lines kept={kept})")


if __name__ == "__main__":
    main()

