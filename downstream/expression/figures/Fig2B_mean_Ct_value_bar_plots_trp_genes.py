import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
import matplotlib.ticker as mticker

# ---------------- symbol ----------------
NOAMP_SYMBOL = "‡"
NOAMP_FONTSIZE = 14

# ---------------- paths ----------------
data_file = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/Normalized_Expression_for_Plot.csv"
out_dir   = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/family_figures"
out_pdf   = os.path.join(out_dir, "RelativeExpression_panel.pdf")
os.makedirs(out_dir, exist_ok=True)

# ---------------- legend handler ----------------
class TextNoAmpHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        txt = plt.Text(
            x0 + width / 2, y0 + height / 2, NOAMP_SYMBOL,
            ha="center", va="center",
            fontsize=NOAMP_FONTSIZE, color="black",
            transform=trans
        )
        return [txt]

# ---------------- fonts / rcParams ----------------
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

sns.set_style("white")

# ---------------- load ----------------
df = pd.read_csv(data_file)
for c in ["Relative_Expression_Mean", "Relative_Expression_SEM"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

def get_family(g):
    n = str(g).lower().strip()
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpv"):  return "TRPV"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("pkd"):   return "PKD"
    if n.startswith("mcoln"): return "MCOLN"
    return np.nan

df["family"] = df["Gene"].apply(get_family)

# drop housekeeping for plot
df_plot = df[~df["family"].isna() & (df["Gene"].astype(str).str.strip() != "Gapdh")].copy()
df_plot["Gene"] = df_plot["Gene"].astype(str).str.strip()

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
df_plot = df_plot.sort_values("Gene").reset_index(drop=True)

palette = {
    "MCOLN": "#4682b4",
    "PKD"  : "#dda0dd",
    "TRPA" : "#696969",
    "TRPC" : "#fc8d62",
    "TRPM" : "#66a61e",
    "TRPV" : "#8da0cb",
}

# no amplification -> keep bars empty + add symbol
no_amp = df_plot["Relative_Expression_Mean"].isna() | (df_plot["Relative_Expression_Mean"] <= 0)
plot_df = df_plot.copy()
plot_df.loc[no_amp, ["Relative_Expression_Mean", "Relative_Expression_SEM"]] = np.nan

# ---------------- figure layout (panel + legend column) ----------------
fig_w = max(7.6, len(gene_order) * 0.32)
fig_h = 4.3
fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[7.6, 1.0], wspace=0.04)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    data=plot_df, x="Gene", y="Relative_Expression_Mean",
    hue="family", hue_order=family_order,
    dodge=False, palette=palette, edgecolor=None, ax=ax
)

# SEM error bars (grey)
xpos  = np.arange(len(plot_df))
yvals = plot_df["Relative_Expression_Mean"].to_numpy(dtype=float)
yerrs = plot_df["Relative_Expression_SEM"].to_numpy(dtype=float)

ax.errorbar(
    x=xpos, y=yvals, yerr=yerrs,
    fmt="none", ecolor="#555555", alpha=0.75,
    capsize=3, lw=1.0, zorder=3
)

# y-axis scaling + headroom
ax.set_yscale("linear")
finite = np.isfinite(yvals)
if finite.any():
    ymax = np.nanmax(yvals[finite] + np.nan_to_num(yerrs[finite], nan=0.0))
    ax.set_ylim(0, ymax * 1.14)
else:
    ax.set_ylim(0, 1)

# show ticks as ×10^-3 on axis
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x*1e3:.0f}"))

# multiplier moved UP a bit more (clearer)
ax.text(
    -0.075, 1.08, r"$\times 10^{-3}$",
    transform=ax.transAxes,
    ha="left", va="bottom",
    fontsize=14, color="black"
)

# labels
ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("Relative expression", labelpad=10)

ax.set_xticks(xpos)
ax.set_xticklabels(gene_order, rotation=60, ha="right")

# ---- FIX: ensure tick marks are visible ----
ax.minorticks_off()
ax.tick_params(
    axis="both", which="major",
    direction="out",
    length=6.0, width=1.2,
    color="black",
    bottom=True, top=False, left=True, right=False
)

# spines
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(1.2)

# horizontal grid (keep your original)
ax.yaxis.grid(True, linestyle=":", linewidth=0.7, color="lightgrey", alpha=0.7)
ax.set_axisbelow(True)

ax.set_xlim(-0.6, len(gene_order) - 0.4)

# ---- vertical dotted background guide lines (ADDED; does NOT change theme) ----
for i in range(len(gene_order)):
    ax.axvline(i, linestyle=":", linewidth=0.7,
               color="lightgrey", alpha=0.7, zorder=0)

# ---- symbol above empty bars ----
y_star = ax.get_ylim()[1] * 0.965
for i, flag in enumerate(no_amp.to_numpy()):
    if flag:
        ax.text(i, y_star, NOAMP_SYMBOL, ha="center", va="bottom",
                fontsize=NOAMP_FONTSIZE, color="black", zorder=5)

# legend (separate axis)
handles, labels = ax.get_legend_handles_labels()
if ax.get_legend() is not None:
    ax.get_legend().remove()

display_map = {"MCOLN": "TRPML", "PKD": "TRPP"}
fam_h, fam_l = [], []
for fam in family_order:
    if fam in labels:
        j = labels.index(fam)
        fam_h.append(handles[j])
        fam_l.append(display_map.get(fam, fam))

noamp_handle = "noamp"  # dummy string handle rendered by handler

ax_leg.legend(
    fam_h + [noamp_handle],
    fam_l + ["No amplification"],
    title="Gene family",
    loc="center left",
    frameon=False,
    title_fontsize=14,
    fontsize=14,
    handletextpad=0.8,
    handler_map={str: TextNoAmpHandler()}
)

fig.subplots_adjust(bottom=0.32, left=0.12, right=0.98, top=0.90)
fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
