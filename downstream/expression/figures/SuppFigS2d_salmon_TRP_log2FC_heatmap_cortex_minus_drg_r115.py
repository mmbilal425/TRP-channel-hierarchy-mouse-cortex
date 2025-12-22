import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager as fm
from pathlib import Path
import re

# -------------------- paths --------------------
in_file = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results/Merged_Salmon_GeneSymbol_TPMs_r115.tsv")
trp_list = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")

out_dir = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
out_dir.mkdir(parents=True, exist_ok=True)

out_pdf = out_dir / "TRP_log2FC_Cortex_minus_DRG_difference_heatmap_Salmon_r115.pdf"
out_tsv = out_dir / "TRP_TPM_Cortex_vs_DRG_Salmon_with_log2FC_r115.tsv"

# -------------------- font --------------------
available = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available else "DejaVu Sans"
plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "pdf.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
})

sns.set_style("white")

# -------------------- TRP list loader (robust) --------------------
def load_trp_symbols(path: Path):
    syms = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("geneid"):
                continue
            if "→" in line:
                sym = line.split("→", 1)[1].strip()
            else:
                sym = re.split(r"[,\t ]+", line)[-1].strip()
            sym = re.split(r"[,\t ]+", sym)[0]
            if sym:
                syms.append(sym)
    return syms

TRP_UP = {s.upper() for s in load_trp_symbols(trp_list)}
if not TRP_UP:
    raise ValueError(f"No TRP symbols loaded from: {trp_list}")

# -------------------- load merged TPM table --------------------
df = pd.read_csv(in_file, sep="\t")

required = {"gene_name", "TPM_Cortex_mean", "TPM_DRG"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {in_file}:\n  {sorted(missing)}")

df["GeneName"] = df["gene_name"].astype(str).str.strip()
df["TPM_Cortex_mean"] = pd.to_numeric(df["TPM_Cortex_mean"], errors="coerce").fillna(0.0)
df["TPM_DRG"]         = pd.to_numeric(df["TPM_DRG"],         errors="coerce").fillna(0.0)

# -------------------- FILTER TO TRP ONLY (fixes huge image bug) --------------------
df["GeneName_up"] = df["GeneName"].str.upper()
df = df[df["GeneName_up"].isin(TRP_UP)].copy()

# nice case for plotting
df["GeneName"] = df["GeneName"].str.capitalize()

# -------------------- compute log2 fold-change --------------------
df["log2FC_Cortex_minus_DRG"] = np.log2((df["TPM_Cortex_mean"] + 1.0) /
                                       (df["TPM_DRG"] + 1.0))

# save TRP-only combined table (same filename as before)
df_out = df[["GeneName", "TPM_Cortex_mean", "TPM_DRG", "log2FC_Cortex_minus_DRG"]].copy()
df_out.to_csv(out_tsv, sep="\t", index=False)

# -------------------- heatmap matrix --------------------
diff = (
    df_out.set_index("GeneName")["log2FC_Cortex_minus_DRG"]
          .to_frame(name="log2FC")
          .sort_values(by="log2FC", ascending=False)
)

abs_max = float(np.nanmax(np.abs(diff.values))) if diff.size else 1.0
if abs_max == 0:
    abs_max = 1.0

custom_cmap = LinearSegmentedColormap.from_list(
    "blue_white_red", ["blue", "white", "darkred"], N=256
)

fig_h = max(6, 0.30 * diff.shape[0])  # safe now because TRP-only
fig = plt.figure(figsize=(6.8, fig_h))
ax  = fig.add_axes([0.12, 0.10, 0.65, 0.80])
cax = fig.add_axes([0.82, 0.25, 0.03, 0.50])

sns.heatmap(
    diff,
    ax=ax,
    cmap=custom_cmap,
    center=0.0,
    vmin=-abs_max, vmax=abs_max,
    cbar=True, cbar_ax=cax,
    annot=True, fmt=".2f",
    xticklabels=False,
    yticklabels=False,
    linewidths=0.0,
    square=False
)

ax.set_xticks([0.5])
ax.set_xticklabels([r"Log$_2$FC (Cortex − DRG)"], rotation=0, ha="center", fontsize=18)

n_genes = diff.shape[0]
ax.set_yticks(np.arange(n_genes) + 0.5)
ax.set_yticklabels(diff.index, fontsize=12)

ax.set_xlabel("")
ax.set_ylabel("TRP genes", fontsize=18)

for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
plt.show()
