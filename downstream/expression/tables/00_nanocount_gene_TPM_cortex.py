#!/usr/bin/env python3
"""
00_nanocount_gene_TPM_cortex.py

Build gene-level TPM tables from NanoCount transcript-level outputs (Nanopore direct RNA).

What it does
------------
1) Reads Ensembl GTF (r115) to map transcript_id -> gene_id / gene_name / gene_biotype
2) Loads NanoCount transcript TSVs for cortex replicates (rep1–rep3)
3) Collapses transcript TPMs to gene-level TPMs per replicate + mean
4) Writes:
   - NanoCount_gene_TPM_all.tsv  (all genes)
   - TRP_TPM_CortexMean_NanoCount.tsv (TRP subset; optional if TRP list exists)

Inputs (HPC paths are defaults)
-------------------------------
- GTF: /g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz
- NanoCount transcript TSVs:
  /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep{1,2,3}/cortex_rep{1,2,3}.nanocount_transcript.tsv
- TRP list (optional): /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt
  Expected columns include: Geneid, gene_name (tab-separated)

Outputs
-------
- /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/NanoCount_gene_TPM_all.tsv
- /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/TRP_TPM_CortexMean_NanoCount.tsv
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

OUT_DIR = Path(
    "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRP_LIST = Path(
    "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt"
)

OUT_TPM_ALL = OUT_DIR / "NanoCount_gene_TPM_all.tsv"
OUT_TPM_TRP = OUT_DIR / "TRP_TPM_CortexMean_NanoCount.tsv"


# ----------------------------
# Helpers
# ----------------------------
def norm_tx_id(tx_id: str) -> str:
    """Strip transcript version suffix (e.g., ENSMUST... .2 -> ENSMUST...)."""
    if not isinstance(tx_id, str):
        return tx_id
    return re.sub(r"\.\d+$", "", tx_id)


def read_gtf_tx2gene(gtf_path: Path) -> pd.DataFrame:
    """Extract transcript_id -> gene_id, gene_name, gene_biotype from GTF."""
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    rows: list[tuple[str, str, str, str]] = []
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
            m_gbt = re.search(r'gene_biotype "([^"]+)"', attrs) or re.search(
                r'gene_type "([^"]+)"', attrs
            )

            if not (m_tx and m_gid):
                continue

            tx_id = norm_tx_id(m_tx.group(1))
            gene_id = m_gid.group(1)
            gene_name = m_gnm.group(1) if m_gnm else gene_id
            gene_biotype = m_gbt.group(1) if m_gbt else ""

            rows.append((tx_id, gene_id, gene_name, gene_biotype))

    df = (
        pd.DataFrame(rows, columns=["transcript_id", "gene_id", "gene_name", "gene_biotype"])
        .drop_duplicates()
    )
    print(f"[GTF] transcript→gene rows: {len(df):,}")
    return df


def load_nc_tpm(tsv_path: Path, label: str) -> pd.DataFrame:
    """Load NanoCount transcript TSV and return transcript_id + TPM_{label}."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"[{label}] NanoCount TSV not found: {tsv_path}")

    df = pd.read_csv(tsv_path, sep="\t")
    cols = {c.lower(): c for c in df.columns}

    tcol = cols.get("transcript_name") or cols.get("transcript") or cols.get("target")
    tpm_col = cols.get("tpm")

    if tcol is None or tpm_col is None:
        raise ValueError(
            f"[{label}] Required columns not found in {tsv_path}\n"
            f"Have columns: {list(df.columns)}\n"
            f"Need: transcript_name (or transcript/target) and TPM"
        )

    out = pd.DataFrame(
        {
            "transcript_id": df[tcol].map(norm_tx_id),
            f"TPM_{label}": pd.to_numeric(df[tpm_col], errors="coerce").fillna(0.0),
        }
    )
    print(f"[{label}] transcripts loaded: {len(out):,}")
    return out


def load_trp_list(trp_path: Path) -> tuple[set[str], set[str]]:
    """
    Return (TRP gene_ids set, TRP gene_names_upper set) from a TRP list TSV.

    Expected columns include: Geneid, gene_name
    """
    df = pd.read_csv(trp_path, sep="\t", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    gene_ids = set(df.get("geneid", pd.Series([], dtype=str)).dropna())
    gene_names_upper = set(
        df.get("gene_name", pd.Series([], dtype=str)).dropna().astype(str).str.upper()
    )
    return gene_ids, gene_names_upper


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    # Map transcript -> gene
    tx2gene = read_gtf_tx2gene(GTF)

    # Load replicate TPMs
    rep_dfs = [load_nc_tpm(REP_TX[k], k) for k in ["rep1", "rep2", "rep3"]]

    # Merge replicate TPMs by transcript_id
    tx = rep_dfs[0]
    for d in rep_dfs[1:]:
        tx = tx.merge(d, on="transcript_id", how="outer")
    tx = tx.fillna(0.0)

    # Attach gene metadata
    tx = tx.merge(tx2gene, on="transcript_id", how="left")
    tx = tx.dropna(subset=["gene_id"])

    # Collapse to gene-level TPM
    tpm_cols = [c for c in tx.columns if c.startswith("TPM_")]
    gene_tpm = (
        tx.groupby(["gene_id", "gene_name", "gene_biotype"], as_index=False)[tpm_cols]
        .sum()
    )
    gene_tpm["TPM_mean"] = gene_tpm[tpm_cols].mean(axis=1)

    # Save ALL genes table
    gene_tpm.sort_values(["gene_name", "gene_id"], kind="stable").to_csv(
        OUT_TPM_ALL, sep="\t", index=False
    )
    print(f"[OK] Saved ALL genes TPM → {OUT_TPM_ALL}")

    # Optional TRP subset
    if TRP_LIST.exists():
        trp_ids, trp_names_upper = load_trp_list(TRP_LIST)

        tmp = gene_tpm.copy()
        tmp["_nameU"] = tmp["gene_name"].fillna("").astype(str).str.upper()

        trp_tpm = tmp[
            tmp["gene_id"].isin(trp_ids) | tmp["_nameU"].isin(trp_names_upper)
        ].drop(columns="_nameU")

        trp_tpm = trp_tpm.rename(columns={"gene_name": "Gene"})
        trp_tpm = trp_tpm.sort_values("Gene", kind="stable")

        trp_tpm.to_csv(OUT_TPM_TRP, sep="\t", index=False)
        print(f"[OK] Saved TRP TPM → {OUT_TPM_TRP}")
    else:
        print(f"[WARN] TRP list not found, skipping TRP export: {TRP_LIST}")

    print("[DONE]")


if __name__ == "__main__":
    main()

