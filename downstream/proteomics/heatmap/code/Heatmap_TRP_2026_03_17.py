#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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
SUB_DIR = os.path.join(SCRIPT_DIR, "Heatmap_TRP")

OUT_CORTEX_PDF = os.path.join(SCRIPT_DIR, "TRP_intensity_heatmap_Cortex_UPDATED.pdf")
OUT_CORTEX_PNG = os.path.join(SCRIPT_DIR, "TRP_intensity_heatmap_Cortex_UPDATED.png")
OUT_ORGANS_PDF = os.path.join(SCRIPT_DIR, "TRP_intensity_heatmap_AllOrgans_UPDATED.pdf")
OUT_ORGANS_PNG = os.path.join(SCRIPT_DIR, "TRP_intensity_heatmap_AllOrgans_UPDATED.png")

# ============================================================
# Orders
# ============================================================
gene_order_cortex = [
    "Pkd2", "Trpc4", "Trpm3", "Trpm7", "Trpv1", "Trpv2", "Trpv4"
]

gene_order_organs = [
    "Mcoln3", "Pkd2", "Trpc4", "Trpm3", "Trpm7", "Trpv1", "Trpv2", "Trpv4"
]

sample_order_cortex = [
    "SN", "Control", "Pellet", "Pelletpellet",
    "FASP_8M", "FASP_SDS", "FASP_Pellet"
]

sample_order_organs = [
    "Cortex", "Cortex_Gel", "DRG", "Kidney",
    "Lung", "Skin", "Testes", "Xeno"
]

# ============================================================
# Labels
# ============================================================
display_labels_cortex = {
    "SN": "SN",
    "Control": "Control/membrane",
    "Pellet": "Pellet",
    "Pelletpellet": "Pellet-pellet",
    "FASP_8M": "FASP 8M",
    "FASP_SDS": "FASP SDS",
    "FASP_Pellet": "FASP pellet",
}

display_labels_organs = {
    "Cortex": "Cortex",
    "Cortex_Gel": "Cortex_Gel",
    "DRG": "DRG",
    "Kidney": "Kidney",
    "Lung": "Lung",
    "Skin": "Skin",
    "Testes": "Testes",
    "Xeno": "Xeno",
}

# ============================================================
# Colour map
# ============================================================
cmap = LinearSegmentedColormap.from_list(
    "natasha_blue",
    ["#B6DBFF", "#1964B0"]
)
cmap = cmap.copy()
cmap.set_bad("#f2f2f2")

# ============================================================
# Find files
# ============================================================
def find_existing_file(filename):
    candidates = [
        os.path.join(SCRIPT_DIR, filename),
        os.path.join(SUB_DIR, filename),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None

# ============================================================
# Load input
# ============================================================
def load_input_tables():
    xlsx_path = find_existing_file("TRP_heatmap_clear_inputs.xlsx")
    cortex_csv = find_existing_file("Fig_c_cortex_input.csv")
    organs_csv = find_existing_file("Fig_d_organs_input.csv")

    if xlsx_path is not None:
        print(f"Using workbook: {xlsx_path}")
        cortex_df = pd.read_excel(xlsx_path, sheet_name="Fig_c_cortex")
        organs_df = pd.read_excel(xlsx_path, sheet_name="Fig_d_organs")
        return clean_matrix(cortex_df), clean_matrix(organs_df)

    if cortex_csv is not None and organs_csv is not None:
        print(f"Using CSV files:\n{cortex_csv}\n{organs_csv}")
        cortex_df = pd.read_csv(cortex_csv)
        organs_df = pd.read_csv(organs_csv)
        return clean_matrix(cortex_df), clean_matrix(organs_df)

    raise FileNotFoundError(
        "Could not find input files.\n\n"
        "Place ONE of these in either:\n"
        f"  {SCRIPT_DIR}\n"
        f"  {SUB_DIR}\n\n"
        "Option 1:\n"
        "  TRP_heatmap_clear_inputs.xlsx\n\n"
        "Option 2:\n"
        "  Fig_c_cortex_input.csv\n"
        "  Fig_d_organs_input.csv"
    )

def clean_matrix(df):
    df = df.copy()

    first_col = df.columns[0]
    df[first_col] = df[first_col].astype(str).str.strip()

    if first_col != "TRP_Gene":
        df = df.rename(columns={first_col: "TRP_Gene"})

    df = df[df["TRP_Gene"].notna()]
    df = df[df["TRP_Gene"] != ""]
    df = df[df["TRP_Gene"].str.lower() != "trp_gene"]

    for col in df.columns[1:]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.drop_duplicates(subset=["TRP_Gene"], keep="first")
    return df.set_index("TRP_Gene")

# ============================================================
# Matrix prep
# ============================================================
def prepare_matrix(df, genes, samples):
    tmp = df.copy().reindex(genes)

    for s in samples:
        if s not in tmp.columns:
            tmp[s] = np.nan

    tmp = tmp[samples]

    for c in tmp.columns:
        tmp[c] = pd.to_numeric(tmp[c], errors="coerce")

    return tmp

def log_transform(mat):
    return np.log10(mat + 1)

# ============================================================
# Draw heatmap
# ============================================================
def draw_heatmap(matrix_df, out_pdf, out_png, title, labels, fig_width, fig_height):
    data = np.ma.masked_invalid(matrix_df.values.astype(float))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor="white")
    im = ax.imshow(data, aspect="auto", cmap=cmap, interpolation="none")

    ax.set_xticks(np.arange(matrix_df.shape[1]))
    ax.set_yticks(np.arange(matrix_df.shape[0]))

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
        fontsize=11
    )

    ax.set_yticklabels(
        [g.upper() for g in matrix_df.index],
        fontsize=13,
        fontweight="bold"
    )

    ax.set_title(title, fontsize=16, pad=10)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.06)
    cbar.set_label(r"$\log_{10}$(mean intensity + 1)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    fig.subplots_adjust(left=0.14, right=0.90, bottom=0.25, top=0.88)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ============================================================
# Main
# ============================================================
def main():
    cortex_df, organs_df = load_input_tables()

    mat_cortex = prepare_matrix(cortex_df, gene_order_cortex, sample_order_cortex)
    mat_cortex = log_transform(mat_cortex)
    labels_c = [display_labels_cortex.get(s, s) for s in sample_order_cortex]

    draw_heatmap(
        mat_cortex,
        OUT_CORTEX_PDF,
        OUT_CORTEX_PNG,
        "Cortex TRP intensity across extraction chemistries",
        labels_c,
        fig_width=9.2,
        fig_height=6.2
    )

    mat_organs = prepare_matrix(organs_df, gene_order_organs, sample_order_organs)
    mat_organs = log_transform(mat_organs)
    labels_d = [display_labels_organs.get(s, s) for s in sample_order_organs]

    draw_heatmap(
        mat_organs,
        OUT_ORGANS_PDF,
        OUT_ORGANS_PNG,
        "TRP intensity across organs (mean intensity)",
        labels_d,
        fig_width=9.8,
        fig_height=6.0
    )

    print("\nDone: heatmaps generated.")
    print(OUT_CORTEX_PDF)
    print(OUT_CORTEX_PNG)
    print(OUT_ORGANS_PDF)
    print(OUT_ORGANS_PNG)

if __name__ == "__main__":
    main()