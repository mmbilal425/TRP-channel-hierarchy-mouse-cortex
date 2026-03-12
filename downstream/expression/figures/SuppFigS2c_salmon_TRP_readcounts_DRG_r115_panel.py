import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from matplotlib import font_manager as fm

# -------------------- paths --------------------
counts_file = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results/Merged_Salmon_GeneSymbol_ReadCounts_r115.tsv")
outdir = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
outdir.mkdir(parents=True, exist_ok=True)
out_pdf = outdir / "TRP_ReadCounts_DRG_Salmon_r115_panel.pdf"

# -------------------- style (same layout) --------------------
font_choice = "Arial" if "Arial" in {f.name for f in fm.fontManager.ttflist} else "DejaVu Sans"
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

# -------------------- load --------------------
df = pd.read_csv(counts_file, sep="\t")

# Expect: gene_id, gene_name, Counts_DRG
need = {"gene_name", "Counts_DRG"}
missing = need - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {counts_file}:\n  {sorted(missing)}")

df["GeneName"] = df["gene_name"].astype(str).str.strip().str.capitalize()
df["Counts_DRG"] = pd.to_numeric(df["Counts_DRG"], errors="coerce").fillna(0.0)

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

# -------------------- canonical order --------------------
gene_order = [
    "mcoln1","mcoln2","mcoln3",
    "pkd2","pkd2l1","pkd2l2",
    "trpa1",
    "trpc1","trpc2","trpc3","trpc4","trpc5","trpc6","trpc7",
    "trpm1","trpm2","trpm3","trpm4","trpm5","trpm6","trpm7","trpm8",
    "trpv1","trpv2","trpv3","trpv4","trpv5","trpv6"
]
upper_map = {g.upper(): i for i, g in enumerate(gene_order)}
df["_order"] = df["GeneName"].str.upper().map(upper_map)
df = df.dropna(subset=["_order"]).sort_values("_order", kind="stable").copy()

def nice_case(s: str) -> str:
    s = str(s).lower()
    if s.startswith("trp") and len(s) >= 4:
        return s[:4].capitalize() + s[4:]
    if s.startswith("pkd"):
        return "Pkd" + s[3:]
    if s.startswith("mcoln"):
        return "Mcoln" + s[5:]
    return s.capitalize()

ordered_labels = [nice_case(g) for g in gene_order]
df["label"] = pd.Categorical(df["GeneName"].astype(str), categories=ordered_labels, ordered=True)
df = df.sort_values("label", kind="stable").copy()

# -------------------- palette --------------------
palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}
hue_order = ["TRPML","TRPP","TRPA","TRPC","TRPM","TRPV"]

value_col = "Counts_DRG"

# -------------------- y-axis ticks --------------------
vmax = float(df[value_col].max()) if len(df) else 0.0
yticks_step = 1000  # keep identical to your DRG plot
top = float(np.ceil(vmax / yticks_step) * yticks_step) if vmax > 0 else yticks_step
yticks = np.arange(0, top + yticks_step, yticks_step)

# -------------------- plot (EXACT same layout) --------------------
fig = plt.figure(figsize=(8.5, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[7.2, 1.2], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    x="label", y=value_col, hue="family",
    data=df, dodge=False, palette=palette,
    hue_order=hue_order, ax=ax
)

ax.margins(x=0.01)
ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("Read count (DRG)", labelpad=10)

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

# ---- OPTIONAL: bar labels (kept exactly like your old script) ----
for p in ax.patches:
    h = p.get_height()
    r = int(round(h))
    if r == 0:
        continue
    ax.annotate(
        f"{r:,}",
        (p.get_x() + p.get_width()/2., h),
        ha="center",
        va="bottom",
        fontsize=7,
        color="black",
        xytext=(0, 4),
        textcoords="offset points",
        clip_on=False
    )

fig.subplots_adjust(bottom=0.32, left=0.12, right=0.98, top=0.90)
fig.savefig(out_pdf, dpi=600, format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
