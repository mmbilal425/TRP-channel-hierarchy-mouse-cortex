# SuppFigS3a_isoquant_known_TRP_isoforms_TPM_bar.py
# Plot known TRP isoform TPM (Supplementary Figure)

import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm

BASE     = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq/Mouse_Cortex_dRNA"
IN_TSV   = f"{BASE}/TRP_known/known_TRP_isoforms_TPM.tsv"
OUT_PDF  = f"{BASE}/TRP_known/known_TRP_isoform_expression_TPM_sorted_by_family.pdf"
os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)

available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "axes.labelsize": 18,
    "xtick.labelsize": 10,
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

df = pd.read_csv(IN_TSV, sep="\t")

def get_family(name: str) -> str:
    n = (str(name) if name is not None else "").lower().strip()
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpv"):  return "TRPV"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("pkd"):   return "PKD"
    if n.startswith("mcoln"): return "MCOLN"
    return "Other"

df["family"] = df["gene_name"].astype(str).apply(get_family)
df = df[df["family"].isin(["MCOLN", "PKD", "TRPC", "TRPM", "TRPV"])].copy()

family_display_map = {"MCOLN": "TRPML", "PKD": "TRPP", "TRPC": "TRPC", "TRPM": "TRPM", "TRPV": "TRPV"}
df["family_disp"] = df["family"].map(family_display_map)

def gene_sort_key(gene: str, family: str):
    g = (str(gene) if gene is not None else "").lower()
    m = re.search(r"(\d+)", g)
    base = int(m.group(1)) if m else 1_000_000
    subrank = 0
    if family == "PKD":
        m2 = re.search(r"l(\d+)$", g)
        if m2:
            subrank = int(m2.group(1))
    return base, subrank, g

tmp = df.apply(lambda r: gene_sort_key(r["gene_name"], r["family"]), axis=1, result_type="expand")
df["gene_num"]     = tmp[0].astype(int)
df["gene_subrank"] = tmp[1].astype(int)

df["label"] = df["gene_name"].astype(str).str.strip() + " (" + df["transcript_id"].astype(str).str.strip() + ")"

family_order = ["TRPML", "TRPP", "TRPC", "TRPM", "TRPV"]
df["family_disp"] = pd.Categorical(df["family_disp"], categories=family_order, ordered=True)

df = (df
      .sort_values(["family_disp", "gene_num", "gene_subrank", "mean_TPM", "transcript_id"],
                   ascending=[True, True, True, False, True])
      .reset_index(drop=True))

labels_order = df["label"].tolist()
df["label"] = pd.Categorical(df["label"], categories=labels_order, ordered=True)

palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}

n = len(df)
fig_w = max(10.5, n * 0.42)
fig_h = 6.2

fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[7.6, 1.15], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    data=df,
    x="label", y="mean_TPM",
    hue="family_disp", hue_order=family_order,
    dodge=False, palette=palette, ax=ax
)

xpos = np.arange(n)
ax.set_xticks(xpos)
ax.set_xticklabels(labels_order, rotation=60, ha="right")

ax.set_xlabel("Gene (Transcript ID)", labelpad=12)
ax.set_ylabel("TPM", labelpad=10)

ax.minorticks_off()
ax.tick_params(axis="both", which="major",
               direction="out", length=5.0, width=1.2, color="black",
               bottom=True, top=False, left=True, right=False)
ax.tick_params(axis="x", pad=2)

for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.xaxis.grid(False)
ax.set_axisbelow(True)

ax.set_xlim(-0.6, n - 0.4)

vals = df["mean_TPM"].to_numpy(dtype=float)
ymax = np.nanmax(vals) if np.isfinite(vals).any() else 1.0
ax.set_ylim(0, ymax * 1.18)

LABEL_MIN_TPM = 1.0
LABEL_FS      = 9
for p in ax.patches:
    h = float(p.get_height())
    if np.isfinite(h) and h >= LABEL_MIN_TPM:
        ax.annotate(f"{h:.1f}",
                    (p.get_x() + p.get_width()/2.0, h),
                    ha="center", va="bottom",
                    fontsize=LABEL_FS, color="black",
                    xytext=(0, 2.5), textcoords="offset points",
                    clip_on=False)

handles, labels = ax.get_legend_handles_labels()
if ax.get_legend() is not None:
    ax.get_legend().remove()

ax_leg.legend(
    handles=handles[:len(family_order)],
    labels=labels[:len(family_order)],
    title="Gene family",
    loc="upper left",
    frameon=False,
    title_fontsize=14,
    fontsize=12
)

fig.subplots_adjust(bottom=0.43, left=0.08, right=0.98, top=0.92)

fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PDF}  (n={len(df)})")
