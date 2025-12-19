#!/usr/bin/env python3
"""
03_add_gene_names_to_novel_TRP_table.py

Add gene symbols (gene_name) to the novel TRP isoforms table using the same Ensembl r115 GTF used for IsoQuant.

Inputs:
  - novel_TRP_isoforms_with_TPM.tsv (from step 01)
  - refs/Mus_musculus.GRCm39.115.gtf (decompressed local copy used for IsoQuant)

Outputs:
  - novel_TRP_isoforms_with_TPM_and_names.tsv
"""

from pathlib import Path
import pandas as pd
import re


# ---------------- USER CONFIG ----------------
BASE = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq")

IN_TSV  = BASE / "Mouse_Cortex_dRNA/TRP_novel/novel_TRP_isoforms_with_TPM.tsv"
REF_GTF = BASE / "refs/Mus_musculus.GRCm39.115.gtf"

OUT_TSV = BASE / "Mouse_Cortex_dRNA/TRP_novel/novel_TRP_isoforms_with_TPM_and_names.tsv"
# ------------------------------------------------


def build_geneid_to_name(gtf_path: Path) -> dict:
    gid2name = {}
    attr_re = re.compile(r'(\w+)\s+"([^"]+)"')

    with open(gtf_path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != "gene":
                continue
            attrs = dict(attr_re.findall(cols[8]))
            gid = attrs.get("gene_id")
            gnm = attrs.get("gene_name")
            if gid and gnm:
                gid2name[gid] = gnm
    return gid2name


def main() -> None:
    if not IN_TSV.exists():
        raise FileNotFoundError(f"Missing input TSV: {IN_TSV}")
    if not REF_GTF.exists():
        raise FileNotFoundError(f"Missing reference GTF: {REF_GTF}")

    gid2name = build_geneid_to_name(REF_GTF)
    print(f"[03] gene_id→gene_name mappings: {len(gid2name):,}")

    df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
    if "gene_id" not in df.columns:
        raise ValueError(f"Expected 'gene_id' in {IN_TSV}. Found: {list(df.columns)}")

    gene_names = df["gene_id"].map(gid2name).fillna(df["gene_id"])
    if "gene_name" in df.columns:
        df["gene_name"] = gene_names
    else:
        df.insert(1, "gene_name", gene_names)

    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[03] Wrote: {OUT_TSV}")
    print(f"[03] Rows: {len(df):,}")


if __name__ == "__main__":
    main()

