import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

# -------------------- style --------------------
available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

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

# -------------------- paths --------------------
RES = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
INFILE  = RES / "TRP_TPM_CortexMean_vs_DRG_r115.tsv"
OUTFILE = RES / "TRP_TPM_Cortex_vs_DRG_dumbbell_log10_r115_panel.pdf"

# -------------------- load --------------------
df = pd.read_csv(INFILE, sep="\t", dtype={"GeneName": str})
required = {"GeneName", "Cortex_TPM_mean", "DRG_TPM"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {INFILE}:\n  {sorted(missing)}")

df["GeneName"] = df["GeneName"].astype(str).str.strip().str.capitalize()
df["Cortex_TPM_mean"] = pd.to_numeric(df["Cortex_TPM_mean"], errors="coerce").fillna(0.0)
df["DRG_TPM"]         = pd.to_numeric(df["DRG_TPM"],         errors="coerce").fillna(0.0)

# -------------------- TRP family --------------------
def get_family(name: str) -> str:
    n = str(name).lower()
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpv"):  return "TRPV"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("pkd"):   return "TRPP"
    if n.startswith("mcoln"): return "TRPML"
    return "Other"

df["family"] = df["GeneName"].apply(get_family)
df = df[df["family"] != "Other"].copy()

palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb"
}

# log10(TPM+1)
df["y_cortex"] = np.log10(df["Cortex_TPM_mean"] + 1.0)
df["y_drg"]    = np.log10(df["DRG_TPM"] + 1.0)

# -------------------- enforce biological gene order --------------------
gene_order = [
    "Mcoln1","Mcoln2","Mcoln3",
    "Pkd2","Pkd2l1","Pkd2l2",
    "Trpa1",
    "Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7",
    "Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8",
    "Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6"
]
upper_map = {g.upper(): i for i, g in enumerate(gene_order)}
df["_ord"] = df["GeneName"].str.upper().map(upper_map)
df = df.dropna(subset=["_ord"]).sort_values("_ord", kind="stable").copy()

order = df["GeneName"].tolist()
xpos = pd.Series(range(len(order)), index=order)
df["_x"] = df["GeneName"].map(xpos)

# -------------------- plot --------------------
sns.set_style("ticks")

fig = plt.figure(figsize=(8.5, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[7.2, 1.2], wspace=0.05)
ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

# connecting lines
for _, r in df.iterrows():
    ax.plot([r["_x"], r["_x"]], [r["y_cortex"], r["y_drg"]],
            color="lightgrey", linewidth=1.1, zorder=1)

# points (DRG filled, Cortex hollow)
for fam, sub in df.groupby("family", sort=False):
    c = palette.get(fam, "grey")
    ax.scatter(sub["_x"], sub["y_drg"], s=40, color=c, edgecolors="none", zorder=3)
    ax.scatter(sub["_x"], sub["y_cortex"], s=40, facecolors="white",
               edgecolors=c, linewidths=1.1, zorder=4)

ax.set_xticks(list(xpos.values))
ax.set_xticklabels(list(xpos.index), rotation=60, ha="right")
ax.set_ylabel("log$_{10}$(TPM + 1)", labelpad=10)
ax.set_xlabel("Gene", labelpad=6)

ax.tick_params(axis="both", which="major",
               direction="out", length=3.8, width=1.0,
               bottom=True, left=True, top=False, right=False)
ax.tick_params(axis="x", pad=2)

ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey")
ax.set_axisbelow(True)
for s in ["left", "bottom"]:
    ax.spines[s].set_linewidth(1.2)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

# -------------------- legend: keep ONLY "Samples" --------------------
sample_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black',
           markerfacecolor='white', markeredgewidth=1.1, markersize=6, label='Cortex (hollow)'),
    Line2D([0], [0], marker='o', linestyle='', color='black',
           markerfacecolor='black', markersize=6, label='DRG (filled)')
]
ax_leg.legend(
    handles=sample_handles,
    title="Samples",
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(-0.2, 1.00),  
    borderaxespad=0.0,
    title_fontsize=14,
    fontsize=12
)

fig.subplots_adjust(bottom=0.32, left=0.12, right=0.98, top=0.90)
fig.savefig(str(OUTFILE), dpi=600, format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUTFILE}")
