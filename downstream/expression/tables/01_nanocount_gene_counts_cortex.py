#!/usr/bin/env python3
"""
01_nanocount_gene_counts_cortex.py

Build gene-level read-count tables from NanoCount transcript-level outputs (Nanopore direct RNA).

What it does
------------
1) Reads Ensembl GTF (r115) to map transcript_id -> gene_id / gene_name
2) Loads NanoCount transcript TSVs for cortex replicates (rep1–rep3)
3) Collapses transcript counts to gene-level counts per replicate + mean
4) Writes:
   - NanoCount_gene_Counts_all.tsv  (all genes)
   - TRP_ReadCounts_cortex_reps_with_mean_FROM_NanoCount.tsv
     (TRP-only; ordered by TRP list; zero-padded)

Inputs (defaults)
-----------------
- GTF: /g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz
- NanoCount transcript TSVs:
  /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep{1,2,3}/cortex_rep{1,2,3}.nanocount_transcript.tsv
- TRP list (optional): /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt
  Expected columns include: Geneid, gene_name

Outputs
-------
- /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/NanoCount_gene_Counts_all.tsv
- /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/TRP_ReadCounts_cortex_reps_with_mean_FROM_NanoCount.tsv
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


# ----------------------------
# Paths (edit if needed)
# ----------------------------
GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")

REP_TX = {
    "rep1": Path(
        "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/"
        "nanocount/results/cortex_rep1/cortex_rep1.nanocount_transcript.tsv"
    ),
    "rep2": Path(
        "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/"
        "nanocount/results/cortex_rep2/cortex_rep2.nanocount_transcript.tsv"
    ),
    "rep3": Path(
        "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/"
        "nanocount/results/cortex_rep3/cortex_rep3.nanocount_transcript.tsv"
    ),
}

TRP_LIST = Path(
    "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt"
)

OUT_DIR = Path(
    "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "NanoCount_gene_Counts_all.tsv"
OUT_TRP = OUT_DIR / "TRP_ReadCounts_cortex_reps_with_mean_FROM_NanoCount.tsv"


# ----------------------------
# Helpers
# ----------------------------
def norm_tx_id(tx_id: str) -> str:
    """Strip transcript version suffix (e.g., ENSMUST... .2 -> ENSMUST...)."""
    if not isinstance(tx_id, str):
        return tx_id
    return re.sub(r"\.\d+$", "", tx_id)


def read_gtf_tx2gene(gtf_path: Path) -> pd.DataFrame:
    """Extract transcript_id -> gene_id, gene_name from GTF."""
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    rows: list[tuple[str, str, str]] = []
    op = gzip.open if gtf_path.suffix == ".gz" else open

    with op(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "transcript":
                continue

            attrs = cols[8]
            m_tx = re.search(r'transcript_id "([^"]+)"', attrs)
            m_gid = re.search(r'gene_id "([^"]+)"', attrs)
            m_gnm = re.search(r'gene_name "([^"]+)"', attrs)

            if not (m_tx and m_gid):
                continue

            tx_id = norm_tx_id(m_tx.group(1))
            gene_id = m_gid.group(1)
            gene_name = m_gnm.group(1) if m_gnm else gene_id

            rows.append((tx_id, gene_id, gene_name))

    df = pd.DataFrame(rows, columns=["transcript_id", "gene_id", "gene_name"]).drop_duplicates()
    print(f"[GTF] transcript→gene rows: {len(df):,}")
    return df


def load_nc_counts(tsv_path: Path, label: str) -> pd.DataFrame:
    """Load NanoCount transcript TSV and return transcript_id + Counts_{label}."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"[{label}] NanoCount TSV not found: {tsv_path}")

    df = pd.read_csv(tsv_path, sep="\t")
    cols = {c.lower(): c for c in df.columns}

    tcol = cols.get("transcript_name") or cols.get("transcript") or cols.get("target")
    ccol = cols.get("est_count") or cols.get("count") or cols.get("reads") or cols.get("numreads")

    if tcol is None or ccol is None:
        raise ValueError(
            f"[{label}] Required columns not found in {tsv_path}\n"
            f"Have columns: {list(df.columns)}\n"
            f"Need: transcript_name (or transcript/target) and est_count (or count/reads)"
        )

    out = pd.DataFrame(
        {
            "transcript_id": df[tcol].map(lambda x: norm_tx_id(str(x)) if pd.notna(x) else x),
            f"Counts_{label}": pd.to_numeric(df[ccol], errors="coerce").fillna(0.0),
        }
    )
    print(f"[{label}] transcripts loaded: {len(out):,}")
    return out


def load_trp_order(trp_path: Path) -> list[str]:
    """Return TRP gene_name order from TRP list TSV."""
    df = pd.read_csv(trp_path, sep="\t", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    if "gene_name" not in df.columns:
        raise ValueError(f"'gene_name' not found in TRP list; columns={list(df.columns)}")
    return df["gene_name"].dropna().astype(str).tolist()


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    tx2gene = read_gtf_tx2gene(GTF)
    rep_dfs = [load_nc_counts(REP_TX[k], k) for k in ["rep1", "rep2", "rep3"]]

    # Merge replicate counts by transcript_id
    tx = rep_dfs[0]
    for d in rep_dfs[1:]:
        tx = tx.merge(d, on="transcript_id", how="outer")
    tx = tx.fillna(0.0)

    # Map transcript -> gene
    tx = tx.merge(tx2gene, on="transcript_id", how="left")
    mapped = int(tx["gene_id"].notna().sum())
    print(f"[MAP] merged rows={len(tx):,} | mapped={mapped:,} | unmapped={len(tx)-mapped:,}")
    tx = tx.dropna(subset=["gene_id"])

    # Collapse to gene-level counts
    count_cols = [c for c in tx.columns if c.startswith("Counts_")]
    gene_counts = (
        tx.groupby(["gene_id", "gene_name"], as_index=False)[count_cols]
        .sum()
    )
    gene_counts["Counts_mean"] = gene_counts[count_cols].mean(axis=1)

    # Save ALL genes
    gene_counts.sort_values(["gene_name", "gene_id"], kind="stable").to_csv(
        OUT_ALL, sep="\t", index=False
    )
    print(f"[OK] Saved ALL genes → {OUT_ALL}")

    # TRP-only (ordered, zero-padded)
    if TRP_LIST.exists():
        trp_order = load_trp_order(TRP_LIST)

        want = pd.DataFrame({"Gene": trp_order})
        have = gene_counts.rename(columns={"gene_name": "Gene"}).copy()

        trp = want.merge(have[["Gene"] + count_cols + ["Counts_mean"]], on="Gene", how="left")
        for c in count_cols + ["Counts_mean"]:
            trp[c] = pd.to_numeric(trp[c], errors="coerce").fillna(0.0)

        trp.to_csv(OUT_TRP, sep="\t", index=False)
        print(f"[OK] Saved TRP counts (ordered, zero-padded) → {OUT_TRP}")
    else:
        print(f"[WARN] TRP list not found; skipped TRP-only export: {TRP_LIST}")

    print("[DONE]")


if __name__ == "__main__":
    main()

