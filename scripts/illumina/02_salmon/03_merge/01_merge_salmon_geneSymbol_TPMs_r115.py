"""
PIPELINE STEP — Illumina RNA-seq (Salmon, Ensembl r115)

Purpose:
Merge per-sample Salmon gene-level TPMs (quant.genes.sf, gene symbols)
into a single matrix with cortex replicates and DRG.

Output:
salmon_output_r115/results/Merged_Salmon_GeneSymbol_TPMs_r115.tsv

NOTE:
Downstream analysis, filtering, and figure-specific tables are handled in:
downstream/expression/tables/
"""

import pandas as pd
import glob, gzip, re
from pathlib import Path

# --- Paths (r115) ---
base = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115")
results = base / "results"
results.mkdir(parents=True, exist_ok=True)

GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")
out_tsv = results / "Merged_Salmon_GeneSymbol_TPMs_r115.tsv"

# --- Collect all gene-level quant files (symbols) ---
files = sorted(glob.glob(str(base / "*_genes_symbol" / "quant.genes.sf")))
if not files:
    raise FileNotFoundError("No files matched *_genes_symbol/quant.genes.sf under salmon_output_r115")

# --- Build gene_name -> gene_id map from GTF ---
def build_name_to_id_map(gtf_path: Path):
    name_to_id = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            m_id = re.search(r'gene_id "([^"]+)"', attrs)
            m_nm = re.search(r'gene_name "([^"]+)"', attrs)
            if m_id and m_nm:
                name_to_id.setdefault(m_nm.group(1), m_id.group(1))
    return name_to_id

print("[gtf] building gene_name → gene_id map …")
name_to_geneid = build_name_to_id_map(GTF)

# --- Load TPMs per sample ---
dfs = []
for f in files:
    sample = Path(f).parent.name.replace("_genes_symbol", "")
    df = pd.read_csv(f, sep="\t", usecols=["Name", "TPM"])
    df = df.rename(columns={"Name": "gene_name", "TPM": f"TPM_{sample}"})
    dfs.append(df)

# --- Merge ---
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="gene_name", how="outer")

merged = merged.fillna(0.0)

# --- Add gene_id ---
merged["gene_id"] = merged["gene_name"].map(name_to_geneid).fillna("")

# --- Cortex mean ---
cortex_cols = [c for c in merged.columns if c.startswith("TPM_cortex-rep")]
merged["TPM_Cortex_mean"] = merged[cortex_cols].mean(axis=1)

# --- DRG ---
if "TPM_DRG-rep1" in merged.columns:
    merged = merged.rename(columns={"TPM_DRG-rep1": "TPM_DRG"})

# --- Final ---
ordered_cols = ["gene_id", "gene_name"] + cortex_cols + ["TPM_Cortex_mean"]
if "TPM_DRG" in merged.columns:
    ordered_cols += ["TPM_DRG"]

final = merged[ordered_cols].sort_values("gene_name", kind="stable")
final.to_csv(out_tsv, sep="\t", index=False)

print(f"[write] {out_tsv}")
