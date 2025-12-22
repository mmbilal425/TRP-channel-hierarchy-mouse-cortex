"""
Figure 1D - Illumina RNA-seq (Salmon) TRP gene expression in adult mouse DRG.
TPM values from DRG replicate 1 (r115 annotation).
"""


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm

# -------------------- paths --------------------
BASE = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115")
RES  = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

IN_TSV   = RES / "Merged_Salmon_GeneSymbol_TPMs_r115.tsv"
TRP_LIST = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")

OUT_PDF = RES / "TRP_TPM_DRG_Salmon_r115_panel.pdf"
YTICKS_STEP = 10

# -------------------- fonts / rcParams (same as your layout) --------------------
available = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available else "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "axes.labelsize": 18,
    "xtick.labelsize": 12,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
})
sns.set_style("white")

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
    return syms

def family(name: str) -> str:
    n = str(name).lower()
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpv"):  return "TRPV"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("pkd"):   return "TRPP"
    if n.startswith("mcoln"): return "TRPML"
    return "Other"

gene_order = [
    "Mcoln1","Mcoln2","Mcoln3",
    "Pkd2","Pkd2l1","Pkd2l2",
    "Trpa1",
    "Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7",
    "Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8",
    "Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6",
]
ord_map = {g.upper(): i for i, g in enumerate(gene_order)}

palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}
hue_order = ["TRPML","TRPP","TRPA","TRPC","TRPM","TRPV"]

# -------------------- load merged table --------------------
df = pd.read_csv(IN_TSV, sep="\t")

required = {"gene_name", "TPM_DRG"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing expected columns in {IN_TSV}:\n  {sorted(missing)}")

df["gene_name"] = df["gene_name"].astype(str)
df["TPM_DRG"] = pd.to_numeric(df["TPM_DRG"], errors="coerce").fillna(0.0)

# -------------------- filter to TRP genes --------------------
trp_syms = load_trp_symbols(TRP_LIST)
TRP_UP = {s.upper() for s in trp_syms}

q = df[df["gene_name"].str.upper().isin(TRP_UP)].copy()

q["GeneName"] = q["gene_name"].str.strip().str.capitalize()
q["_ord"] = q["GeneName"].str.upper().map(ord_map)
q = q.dropna(subset=["_ord"]).sort_values("_ord", kind="stable").copy()

q["family"] = q["GeneName"].apply(family)
q = q[q["family"] != "Other"].copy()

# -------------------- plot (EXACT same layout) --------------------
max_val = float(q["TPM_DRG"].max()) if len(q) else 0.0
top = float(np.ceil(max_val / YTICKS_STEP) * YTICKS_STEP) if max_val > 0 else YTICKS_STEP
yticks = np.arange(0, top + YTICKS_STEP, YTICKS_STEP)

fig = plt.figure(figsize=(8.5, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[7.2, 1.2], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    x="GeneName", y="TPM_DRG", hue="family",
    data=q, dodge=False, palette=palette,
    hue_order=hue_order, ax=ax
)

ax.margins(x=0.01)
ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("TPM (DRG)", labelpad=10)

ax.set_yticks(yticks)
ax.set_ylim(0, top)
ax.margins(y=0.02)

ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey")
ax.set_axisbelow(True)

for spine in ["left", "bottom"]:
    ax.spines[spine].set_linewidth(1.2)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.tick_params(axis="both", which="both", direction="out",
               length=3.8, width=1.0, color="black",
               bottom=True, top=False, left=True, right=False)
ax.tick_params(axis="x", pad=2)

plt.setp(ax.get_xticklabels(), rotation=60, ha="right")

handles, labels = ax.get_legend_handles_labels()
if ax.get_legend() is not None:
    ax.get_legend().remove()

ax_leg.legend(
    handles, labels,
    title="Gene family",
    loc="center left",
    frameon=False,
    borderaxespad=0.0,
    title_fontsize=14,
    fontsize=12
)

fig.subplots_adjust(bottom=0.32, left=0.12, right=0.98, top=0.90)
fig.savefig(OUT_PDF, dpi=600, format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PDF}")
