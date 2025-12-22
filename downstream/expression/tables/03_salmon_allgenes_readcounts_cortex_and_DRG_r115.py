from pathlib import Path
import pandas as pd
import gzip
import re

# ---------------------------- Paths ----------------------------
BASE = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115")
RES  = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

CORTEX_QSFS = [
    BASE / "cortex-rep1_genes_symbol" / "quant.genes.sf",
    BASE / "cortex-rep2_genes_symbol" / "quant.genes.sf",
    BASE / "cortex-rep3_genes_symbol" / "quant.genes.sf",
]
DRG_QSF = BASE / "DRG-rep1_genes_symbol" / "quant.genes.sf"

GTF = Path("/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf.gz")

OUT_TSV = RES / "Merged_Salmon_GeneSymbol_ReadCounts_r115.tsv"

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
            if m_id and m_nm:
                name_to_id.setdefault(m_nm.group(1), m_id.group(1))
    return name_to_id

print("[gtf] building gene_name → gene_id map …")
name_to_geneid = build_name_to_id_map(GTF)
print(f"[gtf] mapped {len(name_to_geneid):,} gene names")

# ---------------------------- Load cortex read counts ----------------------------
dfs = []
for i, qsf in enumerate(CORTEX_QSFS, start=1):
    if not qsf.exists():
        raise FileNotFoundError(f"Missing: {qsf}")
    tmp = pd.read_csv(qsf, sep="\t", usecols=["Name", "NumReads"])
    tmp = tmp.rename(columns={"Name": "gene_name", "NumReads": f"Counts_cortex-rep{i}"})
    dfs.append(tmp)

cortex = dfs[0]
for d in dfs[1:]:
    cortex = cortex.merge(d, on="gene_name", how="outer")

cortex_cols = [c for c in cortex.columns if c.startswith("Counts_cortex-rep")]
for c in cortex_cols:
    cortex[c] = pd.to_numeric(cortex[c], errors="coerce").fillna(0.0)

# Cortex mean (ONLY cortex reps)
cortex["Counts_Cortex_mean"] = cortex[cortex_cols].mean(axis=1)

# ---------------------------- Load DRG read counts ----------------------------
if not DRG_QSF.exists():
    raise FileNotFoundError(f"Missing: {DRG_QSF}")

drg = pd.read_csv(DRG_QSF, sep="\t", usecols=["Name", "NumReads"])
drg = drg.rename(columns={"Name": "gene_name", "NumReads": "Counts_DRG"})
drg["Counts_DRG"] = pd.to_numeric(drg["Counts_DRG"], errors="coerce").fillna(0.0)

# ---------------------------- Merge cortex + DRG ----------------------------
merged = cortex.merge(drg, on="gene_name", how="outer")

# fill missing values with 0
for c in cortex_cols + ["Counts_Cortex_mean", "Counts_DRG"]:
    merged[c] = pd.to_numeric(merged.get(c, 0), errors="coerce").fillna(0.0)

# add gene_id
merged["gene_id"] = merged["gene_name"].map(name_to_geneid).fillna("")

# ---------------------------- Final column order ----------------------------
final = merged[["gene_id", "gene_name"] + cortex_cols + ["Counts_Cortex_mean", "Counts_DRG"]].copy()

# optional: keep counts as floats (Salmon NumReads can be non-integer)
# If you want integers, uncomment next line:
# final[cortex_cols + ["Counts_Cortex_mean", "Counts_DRG"]] = final[cortex_cols + ["Counts_Cortex_mean", "Counts_DRG"]].round().astype(int)

final = final.sort_values("gene_name", kind="stable")

# ---------------------------- Save ----------------------------
final.to_csv(OUT_TSV, sep="\t", index=False)

mapped = (final["gene_id"] != "").sum()
print(f"[write] {OUT_TSV}")
print(f"[map] gene_id mapped for {mapped:,}/{len(final):,} genes ({mapped/len(final):.1%})")
print(final.head(10))
