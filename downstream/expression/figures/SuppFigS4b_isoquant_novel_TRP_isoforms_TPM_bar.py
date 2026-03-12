# SuppFigS3b_isoquant_novel_TRP_isoforms_TPM_bar.py
# Plotting novel TRP isoforms TPM (Supplementary Figure)

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from pathlib import Path

BASE = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq")
IN_TSV  = BASE / "Mouse_Cortex_dRNA/TRP_novel/novel_TRP_isoforms_with_TPM_and_names.tsv"
OUT_PDF = BASE / "Mouse_Cortex_dRNA/TRP_novel/trp_novel_isoforms_TPM.pdf"
OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

# ---------- style ----------
available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 11,
    "axes.labelsize": 14,
    "xtick.labelsize": 9,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "pdf.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
})
sns.set_style("white")

# ---------- load ----------
df = pd.read_csv(IN_TSV, sep="\t")
if "TPM" not in df.columns and "mean_TPM" in df.columns:
    df = df.rename(columns={"mean_TPM": "TPM"})

df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")
df = df.dropna(subset=["TPM"]).copy()

df["label"] = (
    df["gene_name"].astype(str).str.strip()
    + " ("
    + df["transcript_id"].astype(str).str.strip()
    + ")"
)

df = df.sort_values("TPM", ascending=False).reset_index(drop=True)
labels_order = df["label"].tolist()
df["label"] = pd.Categorical(df["label"], categories=labels_order, ordered=True)

palette = {
    "novel_in_catalog": "#000000",
    "novel_not_in_catalog": "#9E9E9E"
}

# ---------- figure ----------
n = len(df)
fig_w = max(4.6, n * 0.50)
fig_h = 3.0

fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[4.8, 1.2], wspace=0.06)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    data=df,
    x="label", y="TPM",
    hue="structural_category",
    hue_order=["novel_in_catalog", "novel_not_in_catalog"],
    palette=palette,
    dodge=False,
    width=0.55,
    ax=ax
)

if ax.get_legend() is not None:
    ax.get_legend().remove()

ax.set_xlabel("Gene (Transcript ID)", labelpad=6)
ax.set_ylabel("TPM", labelpad=6)

ax.set_xticks(range(n))
ax.set_xticklabels(labels_order, rotation=60, ha="right")

ax.set_ylim(0, 8)
ax.set_yticks([0, 2, 4, 6, 8])

ax.minorticks_off()
ax.tick_params(
    axis="both",
    which="major",
    direction="out",
    length=5,
    width=1.1,
    bottom=True,
    left=True,
    top=False,
    right=False
)

for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.spines["left"].set_linewidth(1.1)
ax.spines["bottom"].set_linewidth(1.1)

ax.yaxis.grid(True, linestyle=":", linewidth=0.6, color="lightgrey")
ax.xaxis.grid(False)
ax.set_axisbelow(True)

for p in ax.patches:
    h = float(p.get_height())
    if h >= 0.8:
        ax.annotate(
            f"{h:.1f}",
            xy=(p.get_x() + p.get_width()/2, h),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=8,
            color="black"
        )

legend_handles = [
    Line2D([0],[0], lw=6, color=palette["novel_in_catalog"]),
    Line2D([0],[0], lw=6, color=palette["novel_not_in_catalog"]),
]
ax_leg.legend(
    legend_handles,
    ["novel_in_catalog", "novel_not_in_catalog"],
    title="Structural Category",
    frameon=False,
    loc="upper left",
    title_fontsize=11,
    fontsize=10
)

fig.subplots_adjust(bottom=0.40, left=0.12, right=0.98, top=0.92)
fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PDF}")
