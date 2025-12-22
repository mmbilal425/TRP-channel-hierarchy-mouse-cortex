import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from matplotlib import font_manager as fm

infile = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/TRP_TPM_CortexMean_NanoCount.tsv")
outdir = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results")
outdir.mkdir(parents=True, exist_ok=True)
out_pdf = outdir / "TRP_TPM_CortexMean_NanoCount_panel.pdf"
out_tsv = outdir / "TRP_TPM_CortexMean_NanoCount_ordered.tsv"
YTICKS_STEP = 10

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

trp = pd.read_csv(infile, sep="\t")

if "Gene" in trp.columns:
    trp["GeneName"] = trp["Gene"].astype(str)
elif "gene_name" in trp.columns:
    trp["GeneName"] = trp["gene_name"].astype(str)
else:
    raise ValueError("Neither 'Gene' nor 'gene_name' column found in input TSV.")

for c in ["TPM_rep1", "TPM_rep2", "TPM_rep3", "TPM_mean"]:
    if c not in trp.columns:
        trp[c] = 0.0

def get_family(name):
    n = str(name).lower()
    if n.startswith("trpc"):
        return "TRPC"
    if n.startswith("trpv"):
        return "TRPV"
    if n.startswith("trpm"):
        return "TRPM"
    if n.startswith("trpa"):
        return "TRPA"
    if n.startswith("pkd"):
        return "TRPP"
    if n.startswith("mcoln"):
        return "TRPML"
    return "Other"

def nice_case(s):
    s = str(s).lower()
    if s.startswith("trp") and len(s) >= 4:
        return s[:4].capitalize() + s[4:]
    if s.startswith("pkd"):
        return "Pkd" + s[3:]
    if s.startswith("mcoln"):
        return "Mcoln" + s[5:]
    return s.capitalize()

gene_order = [
    "mcoln1","mcoln2","mcoln3",
    "pkd2","pkd2l1","pkd2l2",
    "trpa1",
    "trpc1","trpc2","trpc3","trpc4","trpc5","trpc6","trpc7",
    "trpm1","trpm2","trpm3","trpm4","trpm5","trpm6","trpm7","trpm8",
    "trpv1","trpv2","trpv3","trpv4","trpv5","trpv6"
]

want = pd.DataFrame({
    "GeneName": [nice_case(g) for g in gene_order],
    "_key": [g.upper() for g in gene_order]
})

trp["family"] = trp["GeneName"].apply(get_family)
have = trp.copy()
have["_key"] = have["GeneName"].str.upper()

merged = want.merge(
    have[["GeneName","family","TPM_rep1","TPM_rep2","TPM_rep3","TPM_mean","_key"]],
    on="_key", how="left", suffixes=("", "_have")
)

for c in ["TPM_rep1","TPM_rep2","TPM_rep3","TPM_mean"]:
    merged[c] = merged[c].fillna(0.0)

merged["family"] = merged["family"].fillna(merged["GeneName"].map(get_family))
merged = merged[merged["family"] != "Other"].copy()

ordered_labels = [nice_case(g) for g in gene_order]
merged["label"] = pd.Categorical(merged["GeneName"], categories=ordered_labels, ordered=True)
merged = merged.sort_values("label", kind="stable")

merged.loc[:, ["GeneName","family","TPM_rep1","TPM_rep2","TPM_rep3","TPM_mean"]].to_csv(out_tsv, sep="\t", index=False)

palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}
hue_order = ["TRPML","TRPP","TRPA","TRPC","TRPM","TRPV"]

d = merged.copy()
max_val = float(d["TPM_mean"].max()) if len(d) else 0.0
top = float(np.ceil(max_val / YTICKS_STEP) * YTICKS_STEP) if max_val > 0 else YTICKS_STEP
yticks = np.arange(0, top + YTICKS_STEP, YTICKS_STEP)

fig = plt.figure(figsize=(8.5, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[7.2, 1.2], wspace=0.05)
xrot = 60
bottom = 0.32
left = 0.12
top_margin = 0.90

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    x="label", y="TPM_mean", hue="family",
    data=d, dodge=False, palette=palette,
    hue_order=hue_order, ax=ax
)
ax.margins(x=0.01)

ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("TPM (Cortex)", labelpad=10)

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

plt.setp(ax.get_xticklabels(), rotation=xrot, ha="right")

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

fig.subplots_adjust(bottom=bottom, left=left, right=0.98, top=top_margin)
fig.savefig(out_pdf, dpi=600, format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
print(f"Saved: {out_tsv}")
