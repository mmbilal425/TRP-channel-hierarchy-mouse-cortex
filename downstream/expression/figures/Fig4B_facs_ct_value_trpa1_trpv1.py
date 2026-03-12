# TRP genes data from FACS

# === FACS TaqMan → Relative expression TSVs (2^-ΔCt vs Gapdh), no plots ===
# Input:
#   /g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/TaqMan_FACS_TRP_raw data.csv
# Outputs (written to facs_results/):
#   1) TaqMan_FACS_TRP_normalized.tsv           (slim: Cell_Type, Gene, Relative_Expression)
#   2) TaqMan_FACS_TRP_relative_expression_full.tsv
#        (full: + Mean_Ct, Gapdh_Ct, Delta_Ct, ND)
#   3) TaqMan_FACS_TRP_relative_expression_QC.txt (GAPDH per cell, ND counts)

import os
import numpy as np
import pandas as pd

# ---------- paths ----------
INFILE = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/TaqMan_FACS_TRP_raw data.csv"
OUTDIR = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/facs_results"
os.makedirs(OUTDIR, exist_ok=True)

OUT_TSV_MIN   = os.path.join(OUTDIR, "TaqMan_FACS_TRP_normalized.tsv")
OUT_TSV_FULL  = os.path.join(OUTDIR, "TaqMan_FACS_TRP_relative_expression_full.tsv")
OUT_QC        = os.path.join(OUTDIR, "TaqMan_FACS_TRP_relative_expression_QC.txt")

# ---------- load ----------
df = pd.read_csv(INFILE)

# keep only the necessary columns and standardize names
colmap = {c: c.strip() for c in df.columns}
df.rename(columns=colmap, inplace=True)

needed = ["Cell_Type", "Gene", "Mean_Ct"]
missing = [c for c in needed if c not in df.columns]
if missing:
    raise ValueError(f"Missing required column(s) in input CSV: {missing}")

df = df[needed].copy()

# Normalize cell-type strings (keeps your labels but unifies capitalization/spacing)
cell_fix = {
    "Neuron": "Neurons",
    "Astrocyte": "Astrocytes",
    "Double positive": "Double_Positive",
    "Double-positive": "Double_Positive",
    "Double_Positive": "Double_Positive",
}
df["Cell_Type"] = df["Cell_Type"].astype(str).str.strip().replace(cell_fix)

# Tidy gene names, keep given casing
df["Gene"] = df["Gene"].astype(str).str.strip()

# Parse Ct; blanks or non-numeric → NaN (treated as ND)
df["Mean_Ct"] = pd.to_numeric(df["Mean_Ct"], errors="coerce")

# ---------- compute ΔCt & 2^-ΔCt vs GAPDH per cell type ----------
# Map GAPDH Ct per cell (average if multiple rows)
is_gapdh = df["Gene"].str.lower().eq("gapdh")
gapdh_by_cell = (
    df.loc[is_gapdh, ["Cell_Type", "Mean_Ct"]]
      .dropna()
      .groupby("Cell_Type", as_index=True)["Mean_Ct"]
      .mean()
      .to_dict()
)

# attach Gapdh_Ct to each row
df["Gapdh_Ct"] = df["Cell_Type"].map(gapdh_by_cell)

# ΔCt and relative expression
df["Delta_Ct"] = df["Mean_Ct"] - df["Gapdh_Ct"]
df["Relative_Expression"] = np.where(
    df[["Mean_Ct", "Gapdh_Ct"]].notna().all(axis=1),
    np.power(2.0, -df["Delta_Ct"]),
    np.nan
)

# ND flag: missing gene Ct OR missing GAPDH for that cell type
df["ND"] = df["Mean_Ct"].isna() | df["Gapdh_Ct"].isna()

# ---------- save ----------
# 1) Full table (helpful for checking)
full_cols = ["Cell_Type","Gene","Mean_Ct","Gapdh_Ct","Delta_Ct","Relative_Expression","ND"]
df.sort_values(["Gene","Cell_Type"], kind="stable")[full_cols].to_csv(
    OUT_TSV_FULL, sep="\t", index=False, float_format="%.8g"
)

# 2) Minimal plotting/analysis TSV (drop GAPDH rows + NaNs)
minimal = (
    df.loc[~is_gapdh, ["Cell_Type","Gene","Relative_Expression"]]
      .dropna(subset=["Relative_Expression"])
      .sort_values(["Gene","Cell_Type"], kind="stable")
)
minimal.to_csv(OUT_TSV_MIN, sep="\t", index=False, float_format="%.8g")

# 3) QC summary text
lines = []
lines.append("=== FACS TaqMan Relative Expression QC ===")
lines.append(f"Input file: {INFILE}")
lines.append(f"Output (minimal TSV): {OUT_TSV_MIN}")
lines.append(f"Output (full TSV): {OUT_TSV_FULL}")
lines.append("")
lines.append("Per-cell-type GAPDH Ct (mean):")
for ct, val in gapdh_by_cell.items():
    lines.append(f"  - {ct}: {val:.3f}")
lines.append("")
total = len(df); nd_total = int(df["ND"].sum())
lines.append(f"Rows total: {total}")
lines.append(f"ND (blank Ct or missing GAPDH): {nd_total} ({nd_total/total:.1%})")
lines.append("")
lines.append("ND by cell type:")
for ct, sub in df.groupby("Cell_Type"):
    n=len(sub); nd=int(sub['ND'].sum())
    frac = (nd/n if n else 0.0)
    lines.append(f"  - {ct}: {nd}/{n} ({frac:.1%})")

with open(OUT_QC, "w") as fh:
    fh.write("\n".join(lines))

print("Saved:", OUT_TSV_MIN)
print("Saved:", OUT_TSV_FULL)
print("Saved QC:", OUT_QC)



# Mean Ct value for TRPA1 & TRPV1 FACS

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib import font_manager as fm

STAR_SIZE = 12
LEGEND_SYMBOL_SIZE = 11
JITTER_WIDTH = 0.22

NOAMP_SYMBOL = "‡"

class TextSymbolHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        t = plt.Text(
            x0 + width / 2, y0 + height / 2, NOAMP_SYMBOL,
            fontsize=LEGEND_SYMBOL_SIZE, ha="center", va="center",
            color="black", transform=trans
        )
        return [t]

data_file = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/TaqMan_FACS_TRP_raw data.csv"
out_dir   = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/facs_results"
out_file  = os.path.join(out_dir, "Mean_Ct_TRPA1_TRPV1_facs_panel.pdf")
os.makedirs(out_dir, exist_ok=True)

available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"
plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "axes.labelsize": 18,
    "xtick.labelsize": 14,
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

df = pd.read_csv(data_file)
df.columns = [c.strip().replace(" ", "_") for c in df.columns]

df["Mean_Ct"] = pd.to_numeric(df["Mean_Ct"], errors="coerce")
df.loc[df["Mean_Ct"].fillna(0) == 0, "Mean_Ct"] = np.nan

cell_fix = {
    "Neuron": "Neurons",
    "Astrocyte": "Astrocytes",
    "Double positive": "Double_Positive",
    "Double-positive": "Double_Positive",
}
df["Cell_Type"] = df["Cell_Type"].astype(str).str.strip().replace(cell_fix)

df = df[df["Gene"].astype(str).str.lower().isin(["trpa1", "trpv1", "gapdh"])].copy()

cell_order = ["Neurons", "Astrocytes", "Double_Positive"]
palette = {"Neurons": "#1964B0", "Astrocytes": "#F4A637", "Double_Positive": "#DEDEDE"}
df["Cell_Type"] = pd.Categorical(df["Cell_Type"], categories=cell_order, ordered=True)

gapdh_map = (
    df[df["Gene"].str.lower() == "gapdh"]
    .groupby("Cell_Type", observed=True)["Mean_Ct"]
    .mean()
    .to_dict()
)

gene_order = ["TRPA1", "TRPV1"]
name_map = {"trpa1": "TRPA1", "trpv1": "TRPV1"}
df["GeneU"] = df["Gene"].astype(str).str.lower().map(name_map)

# Wider spacing between genes
gene_x = {"TRPA1": 0.0, "TRPV1": 1.6}

ct_offset = {
    ct: (i - (len(cell_order) - 1) / 2) * JITTER_WIDTH
    for i, ct in enumerate(cell_order)
}

fig, ax = plt.subplots(figsize=(5.4, 3.8))

ax.set_ylim(20, 45)
y_symbol = ax.get_ylim()[1] - 0.25   

for ct, y in gapdh_map.items():
    if ct in palette and np.isfinite(y):
        ax.axhline(y=y, color=palette[ct], linestyle="--", linewidth=1.2, alpha=0.9, zorder=1)

for ct in cell_order:
    sub = df[(df["GeneU"].isin(gene_order)) & (df["Cell_Type"] == ct) & df["Mean_Ct"].notna()]
    if sub.empty:
        continue
    x = sub["GeneU"].map(gene_x).to_numpy(dtype=float)
    ax.scatter(
        x + ct_offset[ct],
        sub["Mean_Ct"].to_numpy(dtype=float),
        s=90, color=palette[ct], edgecolor="none", linewidth=0, zorder=3
    )

need_symbol = (
    df[(df["GeneU"].isin(gene_order)) & (df["Cell_Type"].isin(cell_order))]
    .groupby(["GeneU", "Cell_Type"], observed=True)["Mean_Ct"]
    .apply(lambda s: s.isna().any())
    .reset_index(name="noamp")
)

for _, r in need_symbol[need_symbol["noamp"]].iterrows():
    g = r["GeneU"]
    ct = r["Cell_Type"]
    if (g not in gene_x) or (ct not in ct_offset):
        continue
    ax.text(
        gene_x[g] + ct_offset[ct], y_symbol, NOAMP_SYMBOL,
        ha="center", va="center",
        fontsize=STAR_SIZE, color="black",
        zorder=5, clip_on=False
    )

ax.set_xlabel("Gene", labelpad=4)
ax.set_ylabel("Mean Ct", labelpad=6)

display_labels = {"TRPA1": "Trpa1", "TRPV1": "Trpv1"}
ax.set_xticks([gene_x[g] for g in gene_order])
ax.set_xticklabels([display_labels[g] for g in gene_order], rotation=60, ha="right")
ax.set_yticks(np.arange(20, 50, 5))

ax.minorticks_off()
ax.tick_params(
    axis="both", which="major",
    direction="out", length=3.8, width=1.0, color="black",
    bottom=True, top=False, left=True, right=False
)
ax.tick_params(axis="x", pad=1)

ax.set_xlim(-0.7, 2.3)

for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.xaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.set_axisbelow(True)

cell_handles = [
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor=palette[ct], markeredgecolor="none",
           markersize=8, label=ct)
    for ct in cell_order
]
gapdh_proxy = Line2D([0], [0], color="#555555", linestyle="--", linewidth=1.2,
                     label="Gapdh Ct (per cell type)")
symbol_proxy = "NOAMP"  # string key for handler_map

ax.legend(
    cell_handles + [gapdh_proxy, symbol_proxy],
    ["Neurons", "Astrocytes", "Double positive",
     "Gapdh Ct (per cell type)", "No amplification"],
    frameon=False, title="Cell type", title_fontsize=14, fontsize=12,
    bbox_to_anchor=(1.02, 1.02), loc="upper left",
    handler_map={str: TextSymbolHandler()}
)

plt.tight_layout()
plt.savefig(out_file, dpi=600, bbox_inches="tight")
plt.close()
print("Saved:", out_file)