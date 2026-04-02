import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

# =========================
# Paths
# =========================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

INFILE = os.path.join(SCRIPT_DIR, "2026_01_30_TRP_rerum.xlsx")
OUTDIR = SCRIPT_DIR

OUTPDF = os.path.join(OUTDIR, "TRPV1_IP_MS_peptide_heatmap_ordered.pdf")
OUTPNG = os.path.join(OUTDIR, "TRPV1_IP_MS_peptide_heatmap_ordered.png")

# =========================
# Style
# =========================
available_fonts = {f.name for f in fm.fontManager.ttflist}
FONT_FAMILY = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 12,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# =========================
# Load data
# =========================
df = pd.read_excel(INFILE)

# Keep TRPV1 only in original Excel order
trpv1 = df[df["Protein.Names"].astype(str).str.contains("TRPV1", na=False)].copy()
trpa1 = df[df["Protein.Names"].astype(str).str.contains("TRPA1", na=False)].copy()

print(f"TRPV1 rows: {len(trpv1)}")
print(f"TRPA1 rows: {len(trpa1)}")

# =========================
# Sample columns
# =========================
bio_cols = [
    "Intensity Control1_pellet",
    "Intensity Control2_pellet",
    "Intensity Control3_pellet",
    "Intensity ControlPellet1",
    "Intensity ControlPellet2",
    "Intensity ControlPellet3",
    "Intensity DRG_pellet",
    "Intensity DRG8M",
    "Intensity DRGPellet",
    "Intensity DRGSample",
    "Intensity DRGSupernatant",
    "Intensity TRPV1_Cortex",
    "Intensity TRPV1_DRG",
]

ref_cols = [
    "Intensity TRPV1_gel",
    "Intensity TRPV1_Xeno",
]

bio_cols = [c for c in bio_cols if c in trpv1.columns]
ref_cols = [c for c in ref_cols if c in trpv1.columns]

# =========================
# Keep peptide order exactly as in Excel
# =========================
row_labels = trpv1["Stripped.Sequence"].astype(str).tolist()

bio_mat = trpv1[bio_cols].copy()
ref_mat = trpv1[ref_cols].copy()

rename_map = {
    "Intensity Control1_pellet": "Ctrl1_pellet",
    "Intensity Control2_pellet": "Ctrl2_pellet",
    "Intensity Control3_pellet": "Ctrl3_pellet",
    "Intensity ControlPellet1": "CtrlPellet1",
    "Intensity ControlPellet2": "CtrlPellet2",
    "Intensity ControlPellet3": "CtrlPellet3",
    "Intensity DRG_pellet": "DRG_pellet",
    "Intensity DRG8M": "DRG8M",
    "Intensity DRGPellet": "DRGPellet",
    "Intensity DRGSample": "DRGSample",
    "Intensity DRGSupernatant": "DRGSupernatant",
    "Intensity TRPV1_Cortex": "TRPV1_Cortex",
    "Intensity TRPV1_DRG": "TRPV1_DRG",
    "Intensity TRPV1_gel": "TRPV1_gel",
    "Intensity TRPV1_Xeno": "TRPV1_Xeno",
}
bio_mat = bio_mat.rename(columns=rename_map)
ref_mat = ref_mat.rename(columns=rename_map)

# =========================
# Log transform
# =========================
bio_plot = np.log10(bio_mat.astype(float) + 1)
ref_plot = np.log10(ref_mat.astype(float) + 1)

combined = pd.concat([bio_plot, ref_plot], axis=1)
vmin = 0
vmax = np.nanmax(combined.values)

# =========================
# Plot
# =========================
n_rows = len(row_labels)
fig_h = max(7.5, n_rows * 0.40)

fig = plt.figure(figsize=(12.5, fig_h), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[14, 3, 0.6], wspace=0.20)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
cax = fig.add_subplot(gs[0, 2])

cmap = plt.cm.Blues.copy()
cmap.set_bad("white")

im1 = ax1.imshow(
    bio_plot.values,
    aspect="auto",
    interpolation="none",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax
)
im2 = ax2.imshow(
    ref_plot.values,
    aspect="auto",
    interpolation="none",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax
)

# -------------------------
# Panel titles
# -------------------------
ax1.set_title("Biological / IP-related runs", pad=14, fontsize=14, weight="bold")
ax2.set_title("Reference runs", pad=14, fontsize=14, weight="bold")

# -------------------------
# X ticks
# -------------------------
ax1.set_xticks(range(bio_plot.shape[1]))
ax1.set_xticklabels(bio_plot.columns, rotation=60, ha="right", fontsize=11)

ax2.set_xticks(range(ref_plot.shape[1]))
ax2.set_xticklabels(ref_plot.columns, rotation=60, ha="right", fontsize=11)

# -------------------------
# Y ticks
# -------------------------
ax1.set_yticks(range(n_rows))
ax1.set_yticklabels(row_labels, fontsize=12)

ax2.set_yticks(range(n_rows))
ax2.set_yticklabels([])

# -------------------------
# Clean axes
# -------------------------
for ax in [ax1, ax2]:
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

# -------------------------
# Colorbar
# -------------------------
cb = fig.colorbar(im2, cax=cax)
cb.set_label(r"$\log_{10}(\mathrm{intensity} + 1)$", fontsize=13)

# -------------------------
# Highlight key peptide
# -------------------------
if "NFALVPLLR" in row_labels:
    idx = row_labels.index("NFALVPLLR")
    for ax, ncols in [(ax1, bio_plot.shape[1]), (ax2, ref_plot.shape[1])]:
        rect = plt.Rectangle(
            (-0.5, idx - 0.5),
            ncols,
            1,
            fill=False,
            linewidth=1.3,
            edgecolor="black"
        )
        ax.add_patch(rect)

# -------------------------
# Overall title
# -------------------------
fig.suptitle(
    "TRPV1 peptide detections in rerun LC–MS/MS analysis",
    fontsize=15,
    weight="bold",
    y=1.02
)

# =========================
# Save
# =========================
plt.savefig(OUTPDF, dpi=600, bbox_inches="tight")
plt.savefig(OUTPNG, dpi=600, bbox_inches="tight")
plt.show()