
import re, gzip
import pandas as pd
from pathlib import Path

# --------- INPUTS ----------
REP_TX = {
    "rep1": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep1/cortex_rep1.nanocount_transcript.tsv"),
    "rep2": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep2/cortex_rep2.nanocount_transcript.tsv"),
    "rep3": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep3/cortex_rep3.nanocount_transcript.tsv"),
}
# Ensembl r115 (full, not ".chr.gtf")
GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")

# TRP list (2 cols: Geneid,gene_name) - used for order + zero-padding
TRP_LIST = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")

OUT_DIR = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_ALL = OUT_DIR / "NanoCount_gene_Counts_all.tsv"
OUT_TRP = OUT_DIR / "TRP_ReadCounts_cortex_reps_with_mean_FROM_NanoCount.tsv"

# --------- helpers ----------
def norm_tx_id(x: str) -> str:
    return re.sub(r"\.\d+$", "", str(x)) if pd.notna(x) else x

def read_gtf_tx2gene(gtf_path: Path) -> pd.DataFrame:
    rows = []
    op = gzip.open if gtf_path.suffix == ".gz" else open
    with op(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"): 
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            a = cols[8]
            m_tx  = re.search(r'transcript_id "([^"]+)"', a)
            m_gid = re.search(r'gene_id "([^"]+)"', a)
            m_gnm = re.search(r'gene_name "([^"]+)"', a)
            if not (m_tx and m_gid):
                continue
            rows.append((norm_tx_id(m_tx.group(1)),
                         m_gid.group(1),
                         m_gnm.group(1) if m_gnm else m_gid.group(1)))
    df = pd.DataFrame(rows, columns=["transcript_id","gene_id","gene_name"]).drop_duplicates()
    print(f"[GTF] transcript→gene rows: {len(df):,}")
    return df

def load_nc_counts(path: Path, rep_label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    tcol = cols.get("transcript_name") or cols.get("transcript") or cols.get("target")
    ccol = cols.get("est_count") or cols.get("count") or cols.get("reads") or cols.get("numreads")
    if tcol is None or ccol is None:
        raise ValueError(f"[{rep_label}] Missing required columns in {path}\n"
                         f"Have: {list(df.columns)}\nNeed: transcript_name + est_count")
    out = pd.DataFrame({
        "transcript_id": df[tcol].map(norm_tx_id),
        f"Counts_{rep_label}": pd.to_numeric(df[ccol], errors="coerce").fillna(0.0)
    })
    print(f"[{rep_label}] transcripts: {len(out):,}")
    return out

# --------- tx→gene map + load reps ----------
tx2gene = read_gtf_tx2gene(GTF)
rep_dfs = [load_nc_counts(REP_TX[k], k) for k in ["rep1","rep2","rep3"]]

# Merge transcripts across reps (outer), fill missing with 0
tx = rep_dfs[0]
for d in rep_dfs[1:]:
    tx = tx.merge(d, on="transcript_id", how="outer")
tx = tx.fillna(0.0)

# Map transcript → gene; drop unmapped
tx = tx.merge(tx2gene, on="transcript_id", how="left")
mapped = tx["gene_id"].notna().sum()
print(f"[MAP] merged rows={len(tx):,} | mapped={mapped:,} | unmapped={len(tx)-mapped:,}")
tx = tx.dropna(subset=["gene_id"])

# --------- collapse to gene ----------
count_cols = [c for c in tx.columns if c.startswith("Counts_")]
gene_counts = (tx.groupby(["gene_id","gene_name"], as_index=False)[count_cols]
                 .sum())
gene_counts["Counts_mean"] = gene_counts[count_cols].mean(axis=1)

# Save ALL genes
gene_counts.sort_values(["gene_name","gene_id"], kind="stable").to_csv(OUT_ALL, sep="\t", index=False)
print(f"[OK] Saved ALL genes → {OUT_ALL}")

# --------- TRP-only, ordered + zero-padded ----------
if TRP_LIST.exists():
    trp_df = pd.read_csv(TRP_LIST, sep="\t", dtype=str)
    trp_df.columns = [c.strip().lower() for c in trp_df.columns]
    if "gene_name" not in trp_df.columns:
        raise ValueError(f"'gene_name' not found in TRP list; columns={list(trp_df.columns)}")
    trp_order = trp_df["gene_name"].dropna().astype(str).tolist()

    want = pd.DataFrame({"Gene": trp_order})
    have = gene_counts.rename(columns={"gene_name": "Gene"}).copy()

    trp = want.merge(have[["Gene"] + count_cols + ["Counts_mean"]], on="Gene", how="left")
    for c in count_cols + ["Counts_mean"]:
        trp[c] = trp[c].fillna(0.0)

    trp = trp[["Gene"] + count_cols + ["Counts_mean"]]
    trp.to_csv(OUT_TRP, sep="\t", index=False)
    print(f"[OK] Saved TRP counts (ordered, zero-padded) → {OUT_TRP}")
else:
    print(f"[WARN] TRP list not found at {TRP_LIST}; skipped TRP-only export.")

print("[DONE]")
