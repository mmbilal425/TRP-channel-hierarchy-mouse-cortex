# loading nanopore data and mean TPM
import pandas as pd
from pathlib import Path

# ONT NanoCount (gene-level) table
ONT_TSV = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/NanoCount_gene_TPM_all.tsv"

ont_raw = pd.read_csv(ONT_TSV, sep="\t")
ont_df = (ont_raw
          .rename(columns={"gene_name":"Gene"})
          .loc[:, ["Gene", "TPM_mean"]]
          .rename(columns={"TPM_mean":"ONT_TPM_mean"}))
ont_df.head(10)

#loading illumina data and mean TPM

# Illumina Salmon quant files (gene symbols)
REP1 = "/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output/cortex-rep1_genes_symbol/quant.genes.sf"
REP2 = "/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output/cortex-rep2_genes_symbol/quant.genes.sf"
REP3 = "/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output/cortex-rep3_genes_symbol/quant.genes.sf"

def load_salmon_gene_tpm(p):
    df = pd.read_csv(p, sep="\t", dtype={"Name":str})
    # in case of duplicates per gene, sum TPM
    return df[["Name","TPM"]].groupby("Name", as_index=False).sum()

df1 = load_salmon_gene_tpm(REP1).rename(columns={"TPM":"TPM_rep1"})
df2 = load_salmon_gene_tpm(REP2).rename(columns={"TPM":"TPM_rep2"})
df3 = load_salmon_gene_tpm(REP3).rename(columns={"TPM":"TPM_rep3"})

illumina_df = (df1.merge(df2, on="Name", how="outer")
                  .merge(df3, on="Name", how="outer"))
for c in ["TPM_rep1","TPM_rep2","TPM_rep3"]:
    illumina_df[c] = illumina_df[c].fillna(0.0)

illumina_df["ILLUMINA_TPM_mean"] = illumina_df[["TPM_rep1","TPM_rep2","TPM_rep3"]].mean(axis=1)
illumina_df = illumina_df.rename(columns={"Name":"Gene"})[["Gene","ILLUMINA_TPM_mean"]]
illumina_df.head(10)


# Load qPCR and compute relative expression (2^-ΔCt) vs GAPDH only

# qPCR summary
QPCR_CSV = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/pcr_data/Updated_Raw_Cq_Summary.csv"

qpcr_raw = pd.read_csv(QPCR_CSV)

# Use Mean_Cq if present; otherwise average any Mean_Cq_* columns
if "Mean_Cq" in qpcr_raw.columns:
    cq = qpcr_raw.rename(columns={"Mean_Cq":"Cq"})
else:
    cq_cols = [c for c in qpcr_raw.columns if c.lower().startswith("mean_cq")]
    if not cq_cols:
        raise ValueError("Could not find Mean_Cq or Mean_Cq_* columns in the qPCR file.")
    cq = qpcr_raw.assign(Cq=qpcr_raw[cq_cols].mean(axis=1))

cq["Gene"] = cq["Gene"].astype(str).str.strip()

# GAPDH reference Cq
ref_row = cq[cq["Gene"].str.upper()=="GAPDH"]
if ref_row.empty:
    raise ValueError("GAPDH not found in qPCR file.")
ref_cq = float(ref_row["Cq"].mean())

# Relative expression: 2^-(Ct_gene - Ct_ref)
qpcr_df = cq.copy()
qpcr_df["QPCR_relExp"] = (2.0 ** (-(qpcr_df["Cq"] - ref_cq))).astype(float)
qpcr_df = qpcr_df[["Gene","QPCR_relExp"]]
qpcr_df.head(10)


# Merge, force TRP order

# Output folder
OUT_DIR = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/corelation_plots/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "TRP_AllMethods_for_correlation.tsv"

# Strict TRP order
trp_order = [
    "Mcoln1","Mcoln2","Mcoln3",
    "Pkd2","Pkd2l1","Pkd2l2",
    "Trpa1",
    "Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7",
    "Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8",
    "Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6"
]

base = pd.DataFrame({"Gene": trp_order})
base["__key__"] = base["Gene"].str.upper()

for d in (ont_df, illumina_df, qpcr_df):
    d["__key__"] = d["Gene"].astype(str).str.strip().str.upper()

merged = (base[["__key__","Gene"]]
          .merge(ont_df[["__key__","ONT_TPM_mean"]], on="__key__", how="left")
          .merge(illumina_df[["__key__","ILLUMINA_TPM_mean"]], on="__key__", how="left")
          .merge(qpcr_df[["__key__","QPCR_relExp"]], on="__key__", how="left"))

for c in ["ONT_TPM_mean","ILLUMINA_TPM_mean","QPCR_relExp"]:
    merged[c] = merged[c].astype(float).fillna(0.0)

final = merged[["Gene","ONT_TPM_mean","ILLUMINA_TPM_mean","QPCR_relExp"]].reset_index(drop=True)
final.to_csv(OUT_TSV, sep="\t", index=False)
print("Saved:", OUT_TSV)
final.head(12)


# Log tranformation

import pandas as pd
import numpy as np
from pathlib import Path

OUT_DIR = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/corelation_plots/results")
in_tsv  = OUT_DIR / "TRP_AllMethods_for_correlation.tsv"
df      = pd.read_csv(in_tsv, sep="\t")

# TPM: log10(TPM + 1)
df["log10_ONT_TPMp1"]      = np.log10(df["ONT_TPM_mean"].astype(float) + 1.0)
df["log10_ILLUM_TPMp1"]    = np.log10(df["ILLUMINA_TPM_mean"].astype(float) + 1.0)

# qPCR: choose a small epsilon based on your data to keep zeros finite
pos = df["QPCR_relExp"][df["QPCR_relExp"] > 0]
eps = max(1e-6, float(pos.min())/10.0) if len(pos) else 1e-6
print(f"qPCR epsilon used: {eps:.2e}")

df["log10_QPCR_relExp"]    = np.log10(df["QPCR_relExp"].astype(float) + eps)

# Save log table
out_log = OUT_DIR / "TRP_AllMethods_for_correlation_LOG.tsv"
df.to_csv(out_log, sep="\t", index=False)
print("Saved:", out_log)
df.head(10)

# Correlation Plotting between Illumina, Nanopore, PCR

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from matplotlib.ticker import MaxNLocator
from matplotlib import font_manager as fm

RESULT_DIR = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/corelation_plots/results")
RESULT_DIR.mkdir(parents=True, exist_ok=True)

infile = RESULT_DIR / "TRP_AllMethods_for_correlation_LOG.tsv"
outfile_pdf = RESULT_DIR / "TRP_correlation_ALLfamilies_N26.pdf"

print("Reading:", infile)
df = pd.read_csv(infile, sep="\t", dtype={"Gene": str})
df["GeneName"] = df["Gene"].astype(str).str.strip()

# Use precomputed log columns
df["illum"] = pd.to_numeric(df["log10_ILLUM_TPMp1"], errors="coerce")
df["nano"]  = pd.to_numeric(df["log10_ONT_TPMp1"], errors="coerce")
df["qpcr"]  = pd.to_numeric(df["log10_QPCR_relExp"], errors="coerce")

# Exclude genes with no qPCR amplification
exclude_genes = {"Mcoln3", "Trpv5"}
df = df[~df["GeneName"].isin(exclude_genes)].copy()

def get_family(name: str) -> str:
    n = str(name).lower()
    if   n.startswith("mcoln"): return "MCOLN"
    elif n.startswith("pkd"):   return "PKD"
    elif n.startswith("trpa"):  return "TRPA"
    elif n.startswith("trpc"):  return "TRPC"
    elif n.startswith("trpm"):  return "TRPM"
    elif n.startswith("trpv"):  return "TRPV"
    else:                       return "Other"

df["family"] = df["GeneName"].apply(get_family)
df = df[df["family"] != "Other"].copy()

palette = {
    "MCOLN": "#4682b4",
    "PKD"  : "#dda0dd",
    "TRPA" : "#696969",
    "TRPC" : "#fc8d62",
    "TRPM" : "#66a61e",
    "TRPV" : "#8da0cb",
}

available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else plt.rcParams.get("font.family", ["DejaVu Sans"])[0]

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "axes.labelsize": 20,
    "xtick.labelsize": 14,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.6,
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.size": 14,
    "ytick.major.size": 14,
})

sns.set_style("white")

MARKER_SIZE    = 58
FIT_LW         = 2.0
UNITY_LW       = 1.6
GRID_LW        = 1.0
HIGHLIGHT_SIZE = 130
HIGHLIGHT_LW   = 1.3
LEADER_LW      = 1.2
LABEL_FONTSIZE = 14

HIGHLIGHT_GENES = ["Trpv2", "Trpm3", "Trpm7", "Trpc4", "Pkd2"]

LABEL_OFFSETS = {
    "panel1": {  # Illumina vs Nanopore
        "Trpv2": (-30, 20),
        "Trpm3": (30, -25),
        "Trpm7": (40, -10),
        "Trpc4": (30, -20),
        "Pkd2":  (-30, -20),
    },
    "panel2": {  # Illumina vs qPCR
        "Trpv2": (34, -10),
        "Trpm3": (34, 14),
        "Trpm7": (-30, 28),
        "Trpc4": (-34, -4),
        "Pkd2":  (-40, 15),
    },
    "panel3": {  # Nanopore vs qPCR
        "Trpv2": (34, -4),
        "Trpm3": (-34, 25),
        "Trpm7": (20, 30),
        "Trpc4": (-36, -10),
        "Pkd2":  (-15, 20),
    },
}

def square_limits(x, y, pad=0.06):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    lo = float(min(x[ok].min(), y[ok].min()))
    hi = float(max(x[ok].max(), y[ok].max()))
    if lo == hi:
        lo, hi = lo - 0.5, hi + 0.5
    rng = hi - lo
    return lo - pad * rng, hi + pad * rng

def add_fit(ax, x, y, same_scale=False, dashed=True):
    slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
    if same_scale:
        xx = np.linspace(min(x.min(), y.min()), max(x.max(), y.max()), 200)
    else:
        xx = np.linspace(x.min(), x.max(), 200)

    ax.plot(
        xx, slope * xx + intercept,
        linestyle="--" if dashed else "-",
        linewidth=FIT_LW,
        color="0.25",
        zorder=2
    )
    return slope, intercept, r_val, p_val

def annotate_stats(ax, x, y, loc="top-left"):
    pr_r, pr_p = stats.pearsonr(x, y)
    sr_r, sr_p = stats.spearmanr(x, y, nan_policy="omit")
    text = (
        f"Pearson r = {pr_r:.2f} (p={pr_p:.1e})\n"
        f"Spearman ρ = {sr_r:.2f} (p={sr_p:.1e})\n"
        f"N = {len(x)}"
    )

    x0, y0, ha, va = 0.03, 0.97, "left", "top"
    if loc == "bottom-right":
        x0, y0, ha, va = 0.97, 0.03, "right", "bottom"
    elif loc == "top-right":
        x0, y0, ha, va = 0.97, 0.97, "right", "top"
    elif loc == "bottom-left":
        x0, y0, ha, va = 0.03, 0.03, "left", "bottom"

    ax.text(
        x0, y0, text,
        transform=ax.transAxes,
        ha=ha, va=va,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.95),
        fontsize=12,
        zorder=10
    )

def enforce_ticks(ax, nbins_major=6):
    ax.minorticks_off()
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_major))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins_major))
    ax.tick_params(
        axis="both", which="major",
        length=8, width=1.3, direction="out",
        bottom=True, top=False, left=True, right=False
    )
    ax.tick_params(axis="both", which="minor", length=0, width=0)

def add_highlights(ax, xcol, ycol, panel_key):
    sub = df[df["GeneName"].isin(HIGHLIGHT_GENES)].copy()

    ax.scatter(
        sub[xcol], sub[ycol],
        s=HIGHLIGHT_SIZE,
        facecolors="none",
        edgecolors="black",
        linewidths=HIGHLIGHT_LW,
        zorder=5
    )

    for _, row in sub.iterrows():
        gene = row["GeneName"]
        x = row[xcol]
        y = row[ycol]
        dx, dy = LABEL_OFFSETS[panel_key].get(gene, (16, 10))

        ax.annotate(
            gene,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if dx >= 0 else "right",
            va="center",
            fontsize=LABEL_FONTSIZE,
            color="black",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9),
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                lw=LEADER_LW,
                shrinkA=0,
                shrinkB=5,
                connectionstyle="arc3,rad=0.0"
            ),
            zorder=6,
            clip_on=False
        )

def panel(ax, xcol, ycol, xlabel, ylabel, same_scale, show_unity, stats_loc="top-left", panel_key="panel1"):
    for fam, sub in df.groupby("family", sort=False):
        ax.scatter(
            sub[xcol], sub[ycol],
            s=MARKER_SIZE,
            color=palette.get(fam, "grey"),
            edgecolors="none",
            alpha=0.95,
            zorder=3
        )

    xv = df[xcol].to_numpy(float)
    yv = df[ycol].to_numpy(float)
    ok = np.isfinite(xv) & np.isfinite(yv)
    x = xv[ok]
    y = yv[ok]

    if same_scale:
        lo, hi = square_limits(x, y, pad=0.07)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    else:
        pad = 0.07
        x_lo, x_hi = x.min(), x.max()
        y_lo, y_hi = y.min(), y.max()
        ax.set_xlim(x_lo - pad * (x_hi - x_lo), x_hi + pad * (x_hi - x_lo))
        ax.set_ylim(y_lo - pad * (y_hi - y_lo), y_hi + pad * (y_hi - y_lo))

    if same_scale and show_unity:
        lo, hi = ax.get_xlim()
        ax.plot([lo, hi], [lo, hi], ls=":", c="0.75", lw=UNITY_LW, zorder=1)

    add_fit(ax, x, y, same_scale=same_scale, dashed=True)
    annotate_stats(ax, x, y, loc=stats_loc)
    add_highlights(ax, xcol, ycol, panel_key)

    ax.set_xlabel(xlabel, fontsize=22, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=22, labelpad=12)

    ax.xaxis.grid(True, linestyle=":", linewidth=GRID_LW, color="0.88")
    ax.yaxis.grid(True, linestyle=":", linewidth=GRID_LW, color="0.88")
    enforce_ticks(ax, nbins_major=6)

    ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax.tick_params(labelsize=13)

fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.6), constrained_layout=True)

panel(
    axes[0], "illum", "nano",
    xlabel=r"Illumina  $\log_{10}(\mathrm{TPM}+1)$",
    ylabel=r"Nanopore  $\log_{10}(\mathrm{TPM}+1)$",
    same_scale=True, show_unity=True, stats_loc="top-left", panel_key="panel1"
)

panel(
    axes[1], "illum", "qpcr",
    xlabel=r"Illumina  $\log_{10}(\mathrm{TPM}+1)$",
    ylabel=r"qPCR  $\log_{10}(2^{-\Delta C_t})$",
    same_scale=False, show_unity=False, stats_loc="bottom-right", panel_key="panel2"
)

panel(
    axes[2], "nano", "qpcr",
    xlabel=r"Nanopore  $\log_{10}(\mathrm{TPM}+1)$",
    ylabel=r"qPCR  $\log_{10}(2^{-\Delta C_t})$",
    same_scale=False, show_unity=False, stats_loc="bottom-right", panel_key="panel3"
)

fig.savefig(outfile_pdf, dpi=600, bbox_inches="tight", format="pdf")
plt.close(fig)

print("Saved:", outfile_pdf)
print("Final N used for correlations:", len(df))
print("Excluded genes:", ", ".join(sorted(exclude_genes)))
