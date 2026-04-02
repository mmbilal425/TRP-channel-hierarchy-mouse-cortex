#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable

# ============================================================
# Font: Arial
# ============================================================
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# ============================================================
# Paths
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(SCRIPT_DIR, "ANU022_proteinGroupsTRP.csv")

OUT_FIG3E_PDF = os.path.join(SCRIPT_DIR, "TRP_heatmap_Fig3E_log10_final.pdf")
OUT_FIG3E_PNG = os.path.join(SCRIPT_DIR, "TRP_heatmap_Fig3E_log10_final.png")

OUT_S8D_PDF = os.path.join(SCRIPT_DIR, "TRP_heatmap_FigS8D_log10_final.pdf")
OUT_S8D_PNG = os.path.join(SCRIPT_DIR, "TRP_heatmap_FigS8D_log10_final.png")

# ============================================================
# Load data
# ============================================================
df = pd.read_csv(INPUT_FILE)
df = df[df["Genes"].notna()].copy()
df = df[df["Genes"] != "Genes"].copy()
df["Genes"] = df["Genes"].astype(str).str.upper().str.strip()
df = df.drop_duplicates()

# ============================================================
# Orders
# ============================================================
gene_order_full = [
    "MCOLN3",
    "PKD2",
    "TRPC4",
    "TRPM3",
    "TRPM7",
    "TRPV1",
    "TRPV2",
    "TRPV4",
]

gene_order_main = [
    "PKD2",
    "TRPC4",
    "TRPM3",
    "TRPM7",
    "TRPV1",
    "TRPV2",
    "TRPV4",
]

sample_order_3E = [
    "SN_1", "SN_2", "SN_3",
    "Control_1", "Control_2", "Control_3",
    "Pellet_1", "Pellet_2", "Pellet_3",
    "Pelletpellet_1", "Pelletpellet_2", "Pelletpellet_3",
    "FASP_8M_1", "FASP_8M_2", "FASP_8M_3",
    "FASP_SDS_1", "FASP_SDS_2", "FASP_SDS_3",
    "FASP_Pellet_1", "FASP_Pellet_2", "FASP_Pellet_3",
    "Gel_1", "Gel_2", "Gel_3",
]

sample_order_S8D = sample_order_3E + [
    "DRG_SN", "DRG_Sample", "DRG_Pellet", "DRG_Pelletpellet", "DRG8M",
    "Kidney_SN", "Kidney_Sample", "Kidney_Pellet", "Kidney_Pelletpellet",
    "Testes_SN", "Testes_Sample", "Testes_Pellet", "Testes_Pelletpellet",
    "Lung8M", "Skin8M",
    "TRPA1_gel", "TRPA1_Xeno", "TRPV1_gel", "TRPV1_Xeno",
]

# ============================================================
# Display labels
# ============================================================
display_labels = {
    "SN_1": "SN_1",
    "SN_2": "SN_2",
    "SN_3": "SN_3",
    "Control_1": "Control_1",
    "Control_2": "Control_2",
    "Control_3": "Control_3",
    "Pellet_1": "Pellet_1",
    "Pellet_2": "Pellet_2",
    "Pellet_3": "Pellet_3",
    "Pelletpellet_1": "Pelletpellet_1",
    "Pelletpellet_2": "Pelletpellet_2",
    "Pelletpellet_3": "Pelletpellet_3",
    "FASP_8M_1": "FASP_8M_1",
    "FASP_8M_2": "FASP_8M_2",
    "FASP_8M_3": "FASP_8M_3",
    "FASP_SDS_1": "FASP_SDS_1",
    "FASP_SDS_2": "FASP_SDS_2",
    "FASP_SDS_3": "FASP_SDS_3",
    "FASP_Pellet_1": "FASP_Pellet_1",
    "FASP_Pellet_2": "FASP_Pellet_2",
    "FASP_Pellet_3": "FASP_Pellet_3",
    "Gel_1": "Gel_1",
    "Gel_2": "Gel_2",
    "Gel_3": "Gel_3",

    "DRG_SN": "DRG_SN",
    "DRG_Sample": "DRG_Control",
    "DRG_Pellet": "DRG_Pellet",
    "DRG_Pelletpellet": "DRG_Pelletpellet",
    "DRG8M": "DRG_8M",

    "Kidney_SN": "Kidney_SN",
    "Kidney_Sample": "Kidney_Control",
    "Kidney_Pellet": "Kidney_Pellet",
    "Kidney_Pelletpellet": "Kidney_Pelletpellet",

    "Testes_SN": "Testes_SN",
    "Testes_Sample": "Testes_Control",
    "Testes_Pellet": "Testes_Pellet",
    "Testes_Pelletpellet": "Testes_Pelletpellet",

    "Lung8M": "Lung_8M",
    "Skin8M": "Skin_8M",

    "TRPA1_gel": "TRPA1_gel",
    "TRPA1_Xeno": "TRPA1_Xeno",
    "TRPV1_gel": "TRPV1_gel",
    "TRPV1_Xeno": "TRPV1_Xeno",
}

# ============================================================
# Colours
# ============================================================
cmap = LinearSegmentedColormap.from_list(
    "natasha_blue",
    ["#B6DBFF", "#1964B0"]
)

CELL_EDGE = "#c7c7c7"
EMPTY_FACE = "white"

# ============================================================
# Prepare matrix
# ============================================================
def prepare_matrix(df_in, genes, samples):
    tmp = df_in.copy()
    tmp = tmp[tmp["Genes"].isin(genes)].copy()
    tmp = tmp.drop_duplicates(subset=["Genes"], keep="first")
    tmp = tmp.set_index("Genes").reindex(genes)

    for s in samples:
        if s not in tmp.columns:
            tmp[s] = np.nan

    tmp = tmp[samples]

    for c in tmp.columns:
        tmp[c] = pd.to_numeric(tmp[c], errors="coerce")

    return tmp

# ============================================================
# Log10 transform
# ============================================================
def log_transform_matrix(mat):
    return np.log10(mat + 1)

# ============================================================
# Mathematical verification (not needed but good to verify)
# ============================================================
def verify_log_transform(raw_df, log_df, name="matrix"):
    raw_flat = raw_df.to_numpy(dtype=float).flatten()
    log_flat = log_df.to_numpy(dtype=float).flatten()

    # Keep only finite paired values
    mask = np.isfinite(raw_flat) & np.isfinite(log_flat)
    raw_valid = raw_flat[mask]
    log_valid = log_flat[mask]

    # Basic check: direct recomputation
    recomputed = np.log10(raw_valid + 1)
    transform_exact = np.allclose(log_valid, recomputed, rtol=0, atol=1e-12)

    # Monotonicity check: if x < y then log10(x+1) < log10(y+1)
    # Equivalent practical check: sorting order preserved
    raw_order = np.argsort(raw_valid, kind="mergesort")
    log_order = np.argsort(log_valid, kind="mergesort")
    rank_order_preserved = np.array_equal(raw_order, log_order)

    # Pairwise monotonicity check
    diffs_raw = np.diff(np.sort(raw_valid))
    diffs_log = np.diff(np.sort(log_valid))
    monotonic_non_decreasing = np.all(diffs_log >= -1e-12)

    # Summary
    print(f"\nVerification for {name}")
    print("-" * (18 + len(name)))
    print(f"Number of finite values checked: {len(raw_valid)}")
    print(f"Exact transform check [log10(raw + 1)]: {transform_exact}")
    print(f"Rank order preserved after transform: {rank_order_preserved}")
    print(f"Monotonic non-decreasing after sorting: {monotonic_non_decreasing}")

    if len(raw_valid) > 0:
        print(f"Raw range: {raw_valid.min():.4f} to {raw_valid.max():.4f}")
        print(f"Log10 range: {log_valid.min():.4f} to {log_valid.max():.4f}")

    return {
        "n_values": len(raw_valid),
        "transform_exact": transform_exact,
        "rank_order_preserved": rank_order_preserved,
        "monotonic_non_decreasing": monotonic_non_decreasing,
    }

# ============================================================
# Draw heatmap
# ============================================================
def draw_grid_heatmap(
    matrix_df,
    out_pdf,
    out_png,
    vmin,
    vmax,
    fig_width,
    fig_height,
    heat_left,
    heat_bottom,
    heat_width,
    heat_height,
    cbar_left,
    cbar_bottom,
    cbar_width,
    cbar_height,
    x_fontsize,
    y_fontsize,
    cbar_fontsize,
    display_labels=None,
    cbar_label=r"$\log_{10}$(intensity + 1)"
):
    genes = list(matrix_df.index)
    samples = list(matrix_df.columns)

    if display_labels is None:
        shown_samples = samples
    else:
        shown_samples = [display_labels.get(s, s) for s in samples]

    n_rows = len(genes)
    n_cols = len(samples)

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor="white")

    ax = fig.add_axes([heat_left, heat_bottom, heat_width, heat_height])
    ax.set_facecolor("white")

    norm = Normalize(vmin=vmin, vmax=vmax)

    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix_df.iloc[i, j]
            face = EMPTY_FACE if pd.isna(val) else cmap(norm(val))
            rect = Rectangle(
                (j, i), 1, 1,
                facecolor=face,
                edgecolor=CELL_EDGE,
                linewidth=0.8
            )
            ax.add_patch(rect)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)

    ax.set_xticklabels(
        shown_samples,
        rotation=90,
        fontsize=x_fontsize,
        fontweight="normal",
        ha="center",
        va="bottom"
    )
    ax.set_yticklabels(
        genes,
        fontsize=y_fontsize,
        fontweight="normal"
    )

    ax.xaxis.tick_top()
    ax.tick_params(axis="x", which="both", length=0, pad=2)
    ax.tick_params(axis="y", which="both", length=0, pad=3)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlabel("")
    ax.set_ylabel("")

    cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=cbar_fontsize, length=0, pad=1)
    cbar.set_label(cbar_label, fontsize=cbar_fontsize + 1, labelpad=3)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ============================================================
# Build matrices
# ============================================================
mat_3E = prepare_matrix(df, gene_order_main, sample_order_3E)
mat_S8D = prepare_matrix(df, gene_order_full, sample_order_S8D)

mat_3E_log = log_transform_matrix(mat_3E)
mat_S8D_log = log_transform_matrix(mat_S8D)

# ============================================================
# Verify transformation
# ============================================================
verify_3E = verify_log_transform(mat_3E, mat_3E_log, name="Figure 3E")
verify_S8D = verify_log_transform(mat_S8D, mat_S8D_log, name="Figure S8D")

# ============================================================
# Dynamic scale in log space
# ============================================================
vmin_3E = np.nanmin(mat_3E_log.values)
vmax_3E = np.nanmax(mat_3E_log.values)

vmin_S8D = np.nanmin(mat_S8D_log.values)
vmax_S8D = np.nanmax(mat_S8D_log.values)

# ============================================================
# Figure 3E
# ============================================================
draw_grid_heatmap(
    matrix_df=mat_3E_log,
    out_pdf=OUT_FIG3E_PDF,
    out_png=OUT_FIG3E_PNG,
    vmin=vmin_3E,
    vmax=vmax_3E,
    fig_width=8.2,
    fig_height=4.9,
    heat_left=0.12,
    heat_bottom=0.16,
    heat_width=0.72,
    heat_height=0.36,
    cbar_left=0.35,
    cbar_bottom=0.86,
    cbar_width=0.28,
    cbar_height=0.028,
    x_fontsize=12,
    y_fontsize=13,
    cbar_fontsize=10,
    display_labels=display_labels,
    cbar_label=r"$\log_{10}$(intensity + 1)"
)

# ============================================================
# Supplementary S8D
# ============================================================
draw_grid_heatmap(
    matrix_df=mat_S8D_log,
    out_pdf=OUT_S8D_PDF,
    out_png=OUT_S8D_PNG,
    vmin=vmin_S8D,
    vmax=vmax_S8D,
    fig_width=11.5,
    fig_height=5.0,
    heat_left=0.08,
    heat_bottom=0.16,
    heat_width=0.80,
    heat_height=0.26,
    cbar_left=0.14,
    cbar_bottom=0.87,
    cbar_width=0.20,
    cbar_height=0.022,
    x_fontsize=10.5,
    y_fontsize=13,
    cbar_fontsize=9,
    display_labels=display_labels,
    cbar_label=r"$\log_{10}$(intensity + 1)"
)

# ============================================================
# Final summary (print)
# ============================================================
print("\nDone: log10-transformed TRP heatmaps exported.")
print(f"Figure 3E scale: {vmin_3E:.3f} to {vmax_3E:.3f}")
print(f"Figure S8D scale: {vmin_S8D:.3f} to {vmax_S8D:.3f}")

all_checks_pass = (
    verify_3E["transform_exact"] and
    verify_3E["rank_order_preserved"] and
    verify_3E["monotonic_non_decreasing"] and
    verify_S8D["transform_exact"] and
    verify_S8D["rank_order_preserved"] and
    verify_S8D["monotonic_non_decreasing"]
)

print(f"\nOverall mathematical verification passed: {all_checks_pass}")