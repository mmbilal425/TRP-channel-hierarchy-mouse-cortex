import pandas as pd, re, gzip
from pathlib import Path

# ---- Paths ----
GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")
REP_TX = {
    "rep1": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep1/cortex_rep1.nanocount_transcript.tsv"),
    "rep2": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep2/cortex_rep2.nanocount_transcript.tsv"),
    "rep3": Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/cortex_rep3/cortex_rep3.nanocount_transcript.tsv"),
}
OUT_DIR = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRP_LIST = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")
USE_TRP = TRP_LIST.exists()

OUT_TPM_ALL = OUT_DIR / "NanoCount_gene_TPM_all.tsv"
OUT_TPM_TRP = OUT_DIR / "TRP_TPM_CortexMean_NanoCount.tsv"

# ---- Helpers ----
def norm_tx_id(x: str) -> str:
    """Strip transcript version suffix (e.g., ENSMUST... .2 → ENSMUST...)."""
    return re.sub(r"\.\d+$", "", x) if isinstance(x, str) else x

def read_gtf_tx2gene(gtf_path: Path) -> pd.DataFrame:
    """Extract transcript_id → gene_id, gene_name, gene_biotype."""
    rows = []
    op = gzip.open if gtf_path.suffix == ".gz" else open
    with op(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"): 
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = cols[8]
            m_tx  = re.search(r'transcript_id "([^"]+)"', attrs)
            m_gid = re.search(r'gene_id "([^"]+)"', attrs)
            m_gnm = re.search(r'gene_name "([^"]+)"', attrs)
            m_gbt = re.search(r'gene_biotype "([^"]+)"', attrs) or re.search(r'gene_type "([^"]+)"', attrs)
            if not (m_tx and m_gid):
                continue
            rows.append((
                norm_tx_id(m_tx.group(1)),
                m_gid.group(1),
                (m_gnm.group(1) if m_gnm else m_gid.group(1)),
                (m_gbt.group(1) if m_gbt else "")
            ))
    df = pd.DataFrame(rows, columns=["transcript_id", "gene_id", "gene_name", "gene_biotype"]).drop_duplicates()
    print(f"[GTF] transcript→gene rows: {len(df):,}")
    return df

def load_nc_tpm(path: Path, label: str) -> pd.DataFrame:
    """Load NanoCount transcript table and extract TPM only."""
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    tcol = cols.get("transcript_name") or cols.get("transcript") or cols.get("target")
    tpm_c = cols.get("tpm")
    if tcol is None or tpm_c is None:
        raise ValueError(f"[{label}] TPM column not found in {path}")
    return pd.DataFrame({
        "transcript_id": df[tcol].map(norm_tx_id),
        f"TPM_{label}": pd.to_numeric(df[tpm_c], errors="coerce").fillna(0.0),
    })

# ---- Load GTF and replicates ----
tx2gene = read_gtf_tx2gene(GTF)
rep_dfs = [load_nc_tpm(REP_TX[k], k) for k in ["rep1", "rep2", "rep3"]]

# Merge replicates and map to gene
tx = rep_dfs[0]
for d in rep_dfs[1:]:
    tx = tx.merge(d, on="transcript_id", how="outer")
tx = tx.fillna(0.0).merge(tx2gene, on="transcript_id", how="left").dropna(subset=["gene_id"])

# ---- Collapse to gene-level TPM ----
tpm_cols = [c for c in tx.columns if c.startswith("TPM_")]
gene_tpm = (
    tx.groupby(["gene_id", "gene_name", "gene_biotype"], as_index=False)[tpm_cols]
      .sum()
)
gene_tpm["TPM_mean"] = gene_tpm[tpm_cols].mean(axis=1)

# ---- Save all genes ----
gene_tpm.sort_values(["gene_name", "gene_id"]).to_csv(OUT_TPM_ALL, sep="\t", index=False)
print(f"[OK] Saved ALL genes TPM → {OUT_TPM_ALL}")

# ---- Optional: TRP subset ----
if USE_TRP:
    trp_df = pd.read_csv(TRP_LIST, sep="\t", dtype=str)
    trp_df.columns = [c.strip().lower() for c in trp_df.columns]
    trp_names = set(trp_df.get("gene_name", pd.Series([], dtype=str)).dropna().str.upper())
    trp_ids = set(trp_df.get("geneid", pd.Series([], dtype=str)).dropna())
    gene_tpm["_nmU"] = gene_tpm["gene_name"].fillna("").astype(str).str.upper()
    trp_tpm = gene_tpm[gene_tpm["gene_id"].isin(trp_ids) | gene_tpm["_nmU"].isin(trp_names)].drop(columns="_nmU")
    trp_tpm.rename(columns={"gene_name": "Gene"}, inplace=True)
    trp_tpm.sort_values("Gene", inplace=True)
    trp_tpm.to_csv(OUT_TPM_TRP, sep="\t", index=False)
    print(f"[OK] Saved TRP TPM → {OUT_TPM_TRP}")

print("[DONE]")
