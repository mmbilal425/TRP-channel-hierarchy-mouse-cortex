import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase

NOAMP_SYMBOL = "‡"
NOAMP_SIZE = 14  # plot symbol size (and legend will match)

class TextSymbolHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        t = plt.Text(
            x0 + width/2, y0 + height/2, NOAMP_SYMBOL,
            ha="center", va="center",
            fontsize=NOAMP_SIZE, color="black",
            transform=trans
        )
        return [t]

available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else plt.rcParams.get("font.family", ["sans-serif"])[0]
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

data_file = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/Updated_Raw_Cq_Summary.csv"
out_dir   = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/family_figures"
out_pdf   = os.path.join(out_dir, "Mean_Ct_SEM_Scatter_panel.pdf")
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(data_file)
for c in df.columns:
    if c.startswith("Mean_Cq"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

def get_family(g):
    n = str(g).lower().strip()
    if n.startswith("trpc"): return "TRPC"
    if n.startswith("trpv"): return "TRPV"
    if n.startswith("trpm"): return "TRPM"
    if n.startswith("trpa"): return "TRPA"
    if n.startswith("pkd"):  return "PKD"
    if n.startswith("mcoln"):return "MCOLN"
    return np.nan

df["family"] = df["Gene"].apply(get_family)

df_plot = df[~df["family"].isna() & (df["Gene"].str.strip() != "Gapdh")].copy()
df_plot["Gene"] = df_plot["Gene"].str.strip()

family_order = ["MCOLN", "PKD", "TRPA", "TRPC", "TRPM", "TRPV"]
order_blocks = {
    "MCOLN": ["Mcoln1", "Mcoln2", "Mcoln3"],
    "PKD"  : ["Pkd2", "Pkd2l1", "Pkd2l2"],
    "TRPA" : ["Trpa1"],
    "TRPC" : ["Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7"],
    "TRPM" : ["Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8"],
    "TRPV" : ["Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6"],
}
gene_order = [g for fam in family_order for g in order_blocks[fam] if g in df_plot["Gene"].values]
df_plot["Gene"] = pd.Categorical(df_plot["Gene"], categories=gene_order, ordered=True)
df_plot.sort_values("Gene", inplace=True, ignore_index=True)

rep_cols = [c for c in df_plot.columns if c.startswith("Mean_Cq_N")]
if "Mean_Cq_SEM" not in df_plot.columns:
    if rep_cols:
        n = df_plot[rep_cols].notna().sum(axis=1).clip(lower=1)
        df_plot["Mean_Cq_SEM"] = df_plot[rep_cols].std(axis=1, ddof=1) / np.sqrt(n)
    else:
        df_plot["Mean_Cq_SEM"] = np.nan
df_plot["Mean_Cq_SEM"] = pd.to_numeric(df_plot["Mean_Cq_SEM"], errors="coerce")

sns.set_style("white")

fig_w = max(7.6, len(gene_order) * 0.32)
fig_h = 4.3
fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[7.6, 1.0], wspace=0.04)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

palette = {
    "MCOLN": "#4682b4",
    "PKD"  : "#dda0dd",
    "TRPA" : "#696969",
    "TRPC" : "#fc8d62",
    "TRPM" : "#66a61e",
    "TRPV" : "#8da0cb",
}

sns.scatterplot(
    data=df_plot, x="Gene", y="Mean_Cq",
    hue="family", palette=palette, hue_order=family_order,
    s=80, edgecolor="none", linewidth=0, zorder=3, ax=ax
)

xpos  = np.arange(len(df_plot))
yvals = df_plot["Mean_Cq"].to_numpy(dtype=float)
yerrs = df_plot["Mean_Cq_SEM"].to_numpy(dtype=float)
ax.errorbar(x=xpos, y=yvals, yerr=yerrs,
            fmt="none", ecolor="#555555", alpha=0.7,
            capsize=3, lw=1, zorder=2)

gapdh_vals = df.loc[df["Gene"].str.strip().eq("Gapdh"), "Mean_Cq"].dropna().to_numpy()
gapdh_handle = None
if gapdh_vals.size:
    ax.axhline(y=float(gapdh_vals[0]), color="red", linestyle="--", linewidth=1.2)
    gapdh_handle = Line2D([0],[0], color="red", linestyle="--", linewidth=1.2, label="Gapdh Ct")

no_amp_genes = ["Mcoln3", "Trpv5"]
symbol_added = False
for g in no_amp_genes:
    if g in gene_order:
        vals = df_plot.loc[df_plot["Gene"]==g, "Mean_Cq"].to_numpy()
        if (vals.size == 0) or np.isnan(vals[0]):
            xi = gene_order.index(g)
            ax.text(xi, 44.4, NOAMP_SYMBOL,
                    ha="center", va="top",
                    fontsize=NOAMP_SIZE, color="black", zorder=5)
            symbol_added = True

ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("Mean Ct value", labelpad=10)
ax.set_xticks(np.arange(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=60, ha="right")
ax.set_yticks(np.arange(16, 46, 5))
ax.set_ylim(16, 45)

ax.minorticks_off()
ax.tick_params(axis="both", which="major",
               direction="out", length=3.8, width=1.0, color="black",
               bottom=True, left=True)

for side in ["top","right"]:
    ax.spines[side].set_visible(False)
for side in ["left","bottom"]:
    ax.spines[side].set_linewidth(1.2)

ax.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.set_axisbelow(True)
ax.set_xlim(-0.6, len(gene_order) - 0.4)

handles, labels = ax.get_legend_handles_labels()
if ax.get_legend() is not None:
    ax.get_legend().remove()

display_map = {"MCOLN": "TRPML", "PKD": "TRPP"}
fam_handles, fam_labels = [], []
for fam in family_order:
    if fam in labels:
        i = labels.index(fam)
        fam_handles.append(handles[i])
        fam_labels.append(display_map.get(fam, fam))

extra_handles, extra_labels = [], []
if gapdh_handle:
    extra_handles.append(gapdh_handle)
    extra_labels.append("Gapdh Ct")

# ---- KEY CHANGE: legend uses text handler, not marker ----
if symbol_added:
    symbol_proxy = "NOAMP"  # any string; handler_map will render ‡
    extra_handles.append(symbol_proxy)
    extra_labels.append("No amplification")

ax_leg.legend(
    handles=fam_handles + extra_handles,
    labels=fam_labels + extra_labels,
    title="Gene family",
    loc="center left",
    frameon=False,
    title_fontsize=14,
    fontsize=12,
    handletextpad=0.8,
    handler_map={str: TextSymbolHandler()}  # renders NOAMP as text at NOAMP_SIZE
)

fig.subplots_adjust(bottom=0.32, left=0.12, right=0.98, top=0.90)
fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
