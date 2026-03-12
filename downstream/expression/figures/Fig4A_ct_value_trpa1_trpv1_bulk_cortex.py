import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

data_file = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/Updated_Raw_Cq_Summary.csv"
out_dir   = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/facs_results"
out_pdf   = os.path.join(out_dir, "Mean_Ct_TRPA1_TRPV1_panel.pdf")
os.makedirs(out_dir, exist_ok=True)

available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else plt.rcParams.get("font.family", ["DejaVu Sans"])[0]
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
df["Gene"] = df["Gene"].astype(str).str.strip()

gene_map = {
    "Trpa1": "Trpa1", "TRPA1": "Trpa1",
    "Trpv1": "Trpv1", "TRPV1": "Trpv1",
    "Gapdh": "Gapdh", "GAPDH": "Gapdh"
}
df["GeneKey"] = df["Gene"].map(gene_map)
df = df[df["GeneKey"].isin(["Trpa1", "Trpv1", "Gapdh"])].copy()

df["Mean_Cq"] = pd.to_numeric(df["Mean_Cq"], errors="coerce")

if "Mean_Cq_SEM" in df.columns:
    df["Mean_Cq_SEM"] = pd.to_numeric(df["Mean_Cq_SEM"], errors="coerce")
else:
    rep_cols = [c for c in df.columns if c.startswith("Mean_Cq_N")]
    if rep_cols:
        for c in rep_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        n = df[rep_cols].notna().sum(axis=1).clip(lower=1)
        df["Mean_Cq_SEM"] = df[rep_cols].std(axis=1, ddof=1) / np.sqrt(n)
    else:
        df["Mean_Cq_SEM"] = np.nan

plot_df = df[df["GeneKey"].isin(["Trpa1", "Trpv1"])].copy()
plot_df = plot_df.groupby("GeneKey", as_index=False).agg(
    MeanCt=("Mean_Cq", "mean"),
    SEM=("Mean_Cq_SEM", "mean")
)

gene_order = ["Trpa1", "Trpv1"]
plot_df["GeneKey"] = pd.Categorical(plot_df["GeneKey"], categories=gene_order, ordered=True)
plot_df = plot_df.sort_values("GeneKey").reset_index(drop=True)

gapdh_vals = df.loc[df["GeneKey"].eq("Gapdh"), "Mean_Cq"].dropna().to_numpy()
gapdh_ct = float(np.nanmean(gapdh_vals)) if gapdh_vals.size else np.nan

fig_w, fig_h = 3.6, 3.8
fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1.4], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

xpos = np.array([0.0, 0.45])  
yvals = plot_df["MeanCt"].to_numpy(dtype=float)
yerrs = plot_df["SEM"].to_numpy(dtype=float)

ax.errorbar(
    x=xpos, y=yvals, yerr=yerrs,
    fmt="none", ecolor="#555555", alpha=0.7,
    capsize=3, lw=1, zorder=2
)

ax.scatter(
    xpos, yvals,
    s=80, color="black",
    edgecolor="none", linewidth=0, zorder=3
)

ax.set_xlabel("Gene", labelpad=4)
ax.set_ylabel("Mean Ct", labelpad=6)

ax.set_xticks(xpos)
ax.set_xticklabels(gene_order, rotation=60, ha="right")

ax.set_ylim(15, 45)
ax.set_yticks(np.arange(15, 46, 5))

ax.minorticks_off()
ax.tick_params(axis="both", which="major",
               direction="out", length=3.8, width=1.0, color="black",
               bottom=True, top=False, left=True, right=False)
ax.tick_params(axis="x", pad=1)

for side in ["top", "right"]:
    ax.spines[side].set_visible(False)
for side in ["left", "bottom"]:
    ax.spines[side].set_linewidth(1.2)

ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.xaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.set_axisbelow(True)

ax.set_xlim(-0.10, 0.55)

handles, labels = [], []
if not np.isnan(gapdh_ct):
    ax.axhline(y=gapdh_ct, color="red", linestyle="--", linewidth=1.2, zorder=1)
    handles.append(Line2D([0], [0], color="red", linestyle="--", linewidth=1.2))
    labels.append("Gapdh Ct")

if handles:
    ax_leg.legend(
        handles=handles,
        labels=labels,
        frameon=False,
        title_fontsize=14,
        fontsize=12,
        loc="upper left"
    )

fig.subplots_adjust(bottom=0.28, left=0.14, right=0.98, top=0.90)
fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
