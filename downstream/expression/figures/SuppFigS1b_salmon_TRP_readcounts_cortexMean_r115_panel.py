import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib import font_manager as fm

# -------------------- paths --------------------
RES = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
counts_file = RES / "Merged_Salmon_GeneSymbol_ReadCounts_r115.tsv"
trp_list    = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")
out_pdf     = RES / "TRP_ReadCounts_CortexMean_Salmon_r115_panel.pdf"

# -------------------- style --------------------
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

# -------------------- helpers --------------------
def load_trp_symbols(path: Path):
    """Robust TRP symbol loader (handles arrow '→' and whitespace/tab)."""
    syms = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or line.lower().startswith("geneid"):
                continue
            if "→" in line:
                sym = line.split("→", 1)[1].strip()
            else:
                sym = re.split(r"[,\t ]+", line)[-1].strip()
            sym = re.split(r"[,\t ]+", sym)[0].strip()
            if sym:
                syms.append(sym)
    return syms

def get_family(name: str) -> str:
    n = str(name).lower()
    if n.startswith("trpc"):  return "TRPC"
    if n.startswith("trpv"):  return "TRPV"
    if n.startswith("trpm"):  return "TRPM"
    if n.startswith("trpa"):  return "TRPA"
    if n.startswith("pkd"):   return "TRPP"
    if n.startswith("mcoln"): return "TRPML"
    return "Other"

def nice_case(s: str) -> str:
    s = str(s).strip().lower()
    if s.startswith("trp") and len(s) >= 4:
        return s[:4].capitalize() + s[4:]
    if s.startswith("pkd"):
        return "Pkd" + s[3:]
    if s.startswith("mcoln"):
        return "Mcoln" + s[5:]
    return s.capitalize()

# -------------------- load --------------------
if not counts_file.exists():
    raise FileNotFoundError(f"Missing: {counts_file}")

df = pd.read_csv(counts_file, sep="\t", dtype={"gene_name": str})
required = {"gene_name", "Counts_Cortex_mean"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {counts_file}:\n  {sorted(missing)}")

# TRP allowlist
trp_syms = load_trp_symbols(trp_list)
TRP_UP = {s.upper() for s in trp_syms}

# filter + clean
trp = df.copy()
trp["GeneName"] = trp["gene_name"].astype(str).str.strip()
trp = trp[trp["GeneName"].str.upper().isin(TRP_UP)].copy()

trp["Reads_mean"] = pd.to_numeric(trp["Counts_Cortex_mean"], errors="coerce").fillna(0.0)

trp["family"] = trp["GeneName"].apply(get_family)
trp = trp[trp["family"] != "Other"].copy()

# enforce your biological TRP order
gene_order = [
    "mcoln1","mcoln2","mcoln3",
    "pkd2","pkd2l1","pkd2l2",
    "trpa1",
    "trpc1","trpc2","trpc3","trpc4","trpc5","trpc6","trpc7",
    "trpm1","trpm2","trpm3","trpm4","trpm5","trpm6","trpm7","trpm8",
    "trpv1","trpv2","trpv3","trpv4","trpv5","trpv6"
]
order_map = {g.upper(): i for i, g in enumerate(gene_order)}
trp["_ord"] = trp["GeneName"].str.upper().map(order_map)
trp = trp.dropna(subset=["_ord"]).sort_values("_ord", kind="stable").copy()

# labels for x-axis
trp["label"] = trp["GeneName"].apply(nice_case)
ordered_labels = [nice_case(g) for g in gene_order]
trp["label"] = pd.Categorical(trp["label"], categories=ordered_labels, ordered=True)
trp = trp.sort_values("label", kind="stable").copy()

# palette + legend order
palette = {
    "TRPML": "#4682b4",
    "TRPP":  "#dda0dd",
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}
hue_order = ["TRPML","TRPP","TRPA","TRPC","TRPM","TRPV"]

# y-axis ticks (same logic as your old script)
vmax = float(trp["Reads_mean"].max()) if len(trp) else 0.0
yticks_step = 1000
top = float(np.ceil(vmax / yticks_step) * yticks_step) if vmax > 0 else yticks_step
yticks = np.arange(0, top + yticks_step, yticks_step)

# -------------------- plot (same layout) --------------------
fig = plt.figure(figsize=(8.5, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[7.2, 1.2], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis("off")

sns.barplot(
    x="label", y="Reads_mean", hue="family",
    data=trp, dodge=False, palette=palette,
    hue_order=hue_order, ax=ax
)

ax.margins(x=0.01)
ax.set_xlabel("Gene", labelpad=6)
ax.set_ylabel("Read count (Cortex mean)", labelpad=10)

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

# legend in the right panel
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

# annotate bar values (same as before)
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
