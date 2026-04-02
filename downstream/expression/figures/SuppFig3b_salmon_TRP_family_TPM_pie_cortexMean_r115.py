# === TRP family TPM pie chart (Illumina, cortex mean) — r115 ===
# Input : Merged_Salmon_GeneSymbol_TPMs_r115.tsv
# Filter: strictly by trp_gene_ids.txt
# Output: TRP_family_TPM_distribution_piechart_r115.pdf (600 dpi)

import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

# -------------------- paths --------------------
RES = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
MERGED = RES / "Merged_Salmon_GeneSymbol_TPMs_r115.tsv"
TRP_LIST = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")
OUT_PDF = RES / "TRP_family_TPM_distribution_piechart_r115.pdf"

# -------------------- font settings --------------------
available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial"
plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
})

# -------------------- helpers --------------------
def load_trp_symbols(path: Path):
    syms = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("geneid"):
                continue
            if "→" in line:
                sym = line.split("→", 1)[1].strip()
            else:
                sym = re.split(r"[,\t ]+", line)[-1].strip()
            sym = re.split(r"[,\t ]+", sym)[0]
            if sym:
                syms.append(sym)
    return {s.upper() for s in syms}

def get_family(name: str) -> str:
    n = name.lower()
    if n.startswith("mcoln"): return "TRPML"
    if n.startswith("pkd"):   return "TRPP"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpv"):  return "TRPV"
    return "Other"

palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}

# -------------------- load TRP allowlist --------------------
TRP_UP = load_trp_symbols(TRP_LIST)
if not TRP_UP:
    raise ValueError("TRP gene list is empty")

# -------------------- load merged TPM table --------------------
if not MERGED.exists():
    raise FileNotFoundError(f"Missing input file: {MERGED}")

df = pd.read_csv(MERGED, sep="\t")

required = {"gene_name", "TPM_Cortex_mean"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

df = df.rename(columns={"gene_name": "GeneName"})
df["GeneName"] = df["GeneName"].astype(str)
df["TPM_Cortex_mean"] = pd.to_numeric(df["TPM_Cortex_mean"], errors="coerce").fillna(0.0)

# -------------------- strict TRP filtering --------------------
df = df[df["GeneName"].str.upper().isin(TRP_UP)].copy()

# -------------------- assign family + aggregate --------------------
df["family"] = df["GeneName"].apply(get_family)
df = df[df["family"] != "Other"].copy()

family_tpm = (
    df.groupby("family", as_index=False)["TPM_Cortex_mean"]
      .sum()
      .rename(columns={"TPM_Cortex_mean": "TPM"})
)

family_order = ["TRPML", "TRPM", "TRPC", "TRPV", "TRPP", "TRPA"]
family_tpm["family"] = pd.Categorical(
    family_tpm["family"], categories=family_order, ordered=True
)
family_tpm = family_tpm.sort_values("family").dropna()

colors = [palette[f] for f in family_tpm["family"]]

# -------------------- plot --------------------
plt.figure(figsize=(7, 7))
wedges, texts, autotexts = plt.pie(
    family_tpm["TPM"],
    labels=family_tpm["family"],
    autopct="%1.1f%%",
    startangle=90,
    counterclock=False,
    colors=colors,
    textprops={"fontsize": 16},
)

# nudge tiny TRPA label
if "TRPA" in family_tpm["family"].values:
    i = list(family_tpm["family"]).index("TRPA")
    x, y = autotexts[i].get_position()
    autotexts[i].set_position((x, y + 0.12))

plt.title("TPM distribution by TRP gene family (Illumina, cortex)", fontsize=18)
plt.tight_layout()
plt.savefig(OUT_PDF, dpi=600, format="pdf", bbox_inches="tight")
plt.show()

print(f"Saved: {OUT_PDF}")
