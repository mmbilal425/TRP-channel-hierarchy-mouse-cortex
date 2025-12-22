# 04_build_known_TRP_isoforms_TPM_table.py
# Build table of known/reference TRP transcript isoforms with TPM from IsoQuant outputs

BASE      = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/Mouse_Cortex_dRNA"

# Prefer the per-replicate table; fall back to single-TPM table if needed
TPM_CAND  = [
    f"{BASE}/Mouse_Cortex_dRNA.transcript_grouped_tpm.tsv",  # rep1/rep2/rep3
    f"{BASE}/Mouse_Cortex_dRNA.transcript_tpm.tsv",          # single TPM
]

# Use the exact GTF you fed to IsoQuant (115, non-chr) copied into refs/
GTF_PATH  = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/refs/Mus_musculus.GRCm39.115.gtf"

TRP_LIST  = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt"

OUT_DIR   = f"{BASE}/TRP_known"
OUT_TSV   = f"{OUT_DIR}/known_TRP_isoforms_TPM.tsv"

# =================== imports ===================
import os, re
import numpy as np
import pandas as pd

os.makedirs(OUT_DIR, exist_ok=True)

# ------------------- helper: load TPM in a robust way -------------------
def load_transcript_tpm(paths):
    """
    Return DataFrame with columns:
      transcript_id, mean_TPM, and (if present) rep* columns.
    Works with transcript_grouped_tpm.tsv (rep columns) or transcript_tpm.tsv (single TPM).
    """
    last_err = None
    for p in paths:
        if not os.path.exists(p):
            continue
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
            # normalize header (remove leading #, trim)
            df.columns = [c.strip().lstrip("#") for c in df.columns]

            # find transcript id column
            id_col = None
            for cand in ("feature_id", "transcript_id", "id"):
                if cand in df.columns:
                    id_col = cand
                    break
            if id_col is None:
                raise ValueError(f"No transcript id column in {p}; columns={df.columns.tolist()}")

            # two possible formats
            if "TPM" in df.columns:
                # single TPM column
                out = df[[id_col, "TPM"]].copy()
                out = out.rename(columns={id_col: "transcript_id"})
                out["mean_TPM"] = pd.to_numeric(out["TPM"], errors="coerce")
                out = out[["transcript_id", "mean_TPM"]]
            else:
                # per-replicate columns (keep numeric only)
                rep_cols = [c for c in df.columns if c != id_col]
                num_cols = [c for c in rep_cols if pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.9]
                if not num_cols:
                    raise ValueError(f"No numeric TPM columns in {p}; columns={df.columns.tolist()}")

                num_df = df[num_cols].apply(pd.to_numeric, errors="coerce")
                out = df.rename(columns={id_col: "transcript_id"})[["transcript_id"] + num_cols].copy()
                out["mean_TPM"] = num_df.mean(axis=1)
                out = out[["transcript_id", "mean_TPM"] + num_cols]

            out["transcript_id"] = out["transcript_id"].astype(str)
            print(f"[TPM] loaded: {p} (rows={len(out)})")
            return out
        except Exception as e:
            last_err = e
            print(f"[warn] failed to read {p}: {e}")
    raise last_err if last_err else FileNotFoundError("No TPM file found in TPM_CAND")

# ------------------- helper: map transcript → gene (+ name) --------------
def build_tx_maps_from_gtf(gtf_path):
    gtf = pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, dtype=str,
        names=["chr","source","feature","start","end","score","strand","frame","attributes"]
    )
    # transcript -> gene_id
    tx = gtf[gtf["feature"] == "transcript"].copy()
    tx["transcript_id"] = tx["attributes"].str.extract(r'transcript_id "([^"]+)"')[0]
    tx["gene_id"]       = tx["attributes"].str.extract(r'gene_id "([^"]+)"')[0]
    tx = tx.dropna(subset=["transcript_id","gene_id"]).drop_duplicates(subset=["transcript_id"])

    # gene_id -> gene_name
    g = gtf[gtf["feature"] == "gene"].copy()
    g["gene_id"]   = g["attributes"].str.extract(r'gene_id "([^"]+)"')[0]
    g["gene_name"] = g["attributes"].str.extract(r'gene_name "([^"]+)"')[0]
    gid2name = dict(zip(g["gene_id"], g["gene_name"]))

    return tx[["transcript_id","gene_id"]], gid2name

# ------------------- helper: TRP list ------------------------------------
def load_trp_gene_ids(path):
    trp = pd.read_csv(path, sep="\t", dtype=str)
    col = "Geneid" if "Geneid" in trp.columns else trp.columns[0]
    return set(trp[col].astype(str))

# =================== build table ===================
tpm = load_transcript_tpm(TPM_CAND)
txmap, gid2name = build_tx_maps_from_gtf(GTF_PATH)
trp_set = load_trp_gene_ids(TRP_LIST)

# Keep only reference/known transcript IDs (ignore IsoQuant "transcriptXXXX" models)
KNOWN_PREFIXES = ("ENSMUST", "ENST", "NM_", "XM_", "XR_")
tpm_known = tpm[tpm["transcript_id"].str.startswith(KNOWN_PREFIXES)].copy()

# Add gene_id and gene_name
tpm_known = tpm_known.merge(txmap, on="transcript_id", how="left")
tpm_known["gene_name"] = tpm_known["gene_id"].map(gid2name).fillna(tpm_known["gene_id"])

# Filter to TRP genes and TPM > 0
known_trp = tpm_known[tpm_known["gene_id"].isin(trp_set)].copy()
known_trp = known_trp[known_trp["mean_TPM"] > 0]

# Tidy columns + sort
rep_cols = [c for c in known_trp.columns if c.lower().startswith("rep")]
cols_out = ["gene_id","gene_name","transcript_id","mean_TPM"] + rep_cols
known_trp = known_trp[cols_out].sort_values("mean_TPM", ascending=False)

# Save + summary
known_trp.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\n# known TRP isoforms with TPM > 0: {len(known_trp)}")
print(f"[saved] {OUT_TSV}")
print("\nTop rows:")
print(known_trp.head(10))
