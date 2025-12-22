"""
Supplementary Table: Salmon (Illumina) gene-level read counts for adult mouse cortex.
Inputs: quant.genes.sf from cortex rep1 - rep3 (gene symbols), r115 GTF for gene_id mapping.
Output: AllGenes_ReadCounts_Cortex_reps_withMean_Salmon_r115.tsv
"""

from pathlib import Path
import pandas as pd
import gzip
import re

# ---------------------------- Paths ----------------------------
BASE = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115")
RES  = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

CORTEX_DIRS = [
    BASE / "cortex-rep1_genes_symbol",
    BASE / "cortex-rep2_genes_symbol",
    BASE / "cortex-rep3_genes_symbol",
]

GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")

OUT_TSV = RES / "AllGenes_ReadCounts_Cortex_reps_withMean_Salmon_r115.tsv"

# ---------------------------- Build gene_name -> gene_id map ----------------------------
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
            if not m_id or not m_nm:
                continue

            # keep first occurrence
            name_to_id.setdefault(m_nm.group(1), m_id.group(1))

    return name_to_id

print("[gtf] building gene_name → gene_id map …")
name_to_geneid = build_name_to_id_map(GTF)
print(f"[gtf] mapped {len(name_to_geneid):,} gene names")

# ---------------------------- Load Salmon read counts per rep ----------------------------
dfs = []
for i, d in enumerate(CORTEX_DIRS, start=1):
    qfile = d / "quant.genes.sf"
    if not qfile.exists():
        raise FileNotFoundError(f"Missing: {qfile}")

    tmp = pd.read_csv(qfile, sep="\t", usecols=["Name", "NumReads"])
    tmp = tmp.rename(columns={"Name": "gene_name", "NumReads": f"Counts_rep{i}"})
    dfs.append(tmp)

# ---------------------------- Merge reps ----------------------------
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="gene_name", how="outer")

count_cols = [c for c in merged.columns if c.startswith("Counts_rep")]
for c in count_cols:
    merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)

# ---------------------------- Add gene_id ----------------------------
merged["gene_id"] = merged["gene_name"].map(name_to_geneid).fillna("")

# ---------------------------- Compute mean (LAST column) ----------------------------
merged["Counts_mean"] = merged[count_cols].mean(axis=1)

# ---------------------------- Final table ----------------------------
final = merged[["gene_id", "gene_name"] + count_cols + ["Counts_mean"]].copy()

# optional rounding: keep mean as float (recommended)
final = final.sort_values("gene_name", kind="stable")

# ---------------------------- Save ----------------------------
final.to_csv(OUT_TSV, sep="\t", index=False)

mapped = (final["gene_id"] != "").sum()
print(f"[write] {OUT_TSV}")
print(f"[map] gene_id mapped for {mapped:,}/{len(final):,} genes ({mapped/len(final):.1%})")
print(final.head(10))
