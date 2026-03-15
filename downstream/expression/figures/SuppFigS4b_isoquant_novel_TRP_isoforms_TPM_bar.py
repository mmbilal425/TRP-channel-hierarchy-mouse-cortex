# === Novel TRP isoform expression (TPM) — panel style matched to known isoform figure ===

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm

# =================== paths ===================
BASE = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/Isoquant_results_from_fastq"
IN_TSV = f"{BASE}/Mouse_Cortex_dRNA/TRP_novel/novel_TRP_isoforms_with_TPM_and_names.tsv"
OUT_PDF = f"{BASE}/Mouse_Cortex_dRNA/TRP_novel/trp_novel_isoforms_TPM_panel_matched.pdf"
os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)

# =================== style ===================
available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "axes.labelsize": 18,
    "xtick.labelsize": 10,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
    "legend.title_fontsize": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
})
sns.set_style("white")

# =================== load ===================
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

# =================== plot ===================
n = len(df)
fig, ax = plt.subplots(figsize=(4.0, 5.2))

sns.barplot(
    data=df,
    x="label",
    y="TPM",
    hue="structural_category",
    hue_order=["novel_in_catalog", "novel_not_in_catalog"],
    dodge=False,
    palette=palette,
    width=0.82,
    ax=ax
)

# legend
leg = ax.legend(
    title="Structural Category",
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.36, 1.02),
    borderaxespad=0.0,
    handlelength=1.6,
    handletextpad=0.5,
    labelspacing=0.3,
    prop={"size": 10}
)
leg.get_title().set_fontsize(11)

# axes
ax.set_xticks(np.arange(n))
ax.set_xticklabels(labels_order, rotation=45, ha="right")
ax.set_xlabel("Gene (Transcript ID)", labelpad=12)
ax.set_ylabel("TPM", labelpad=10)

ax.set_ylim(0, 8)
ax.set_yticks([0, 2, 4, 6, 8])
ax.set_xlim(-0.6, n - 0.4)

ax.minorticks_off()
ax.tick_params(
    axis="both",
    which="major",
    direction="out",
    length=5.0,
    width=1.2,
    color="black",
    bottom=True,
    top=False,
    left=True,
    right=False
)
ax.tick_params(axis="x", pad=2)

# spines
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

# grid
ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.xaxis.grid(False)
ax.set_axisbelow(True)

# layout
fig.subplots_adjust(left=0.18, right=0.98, top=0.82, bottom=0.38)

fig.savefig(OUT_PDF, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PDF}")