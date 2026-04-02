#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ==========================================================
# PATHS
# ==========================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = SCRIPT_DIR

PROT_FILE    = os.path.join(BASE_DIR, "ANU022_proteinGroups.csv")
DESIGN2_FILE = os.path.join(BASE_DIR, "ANU022_DEP_Design2.csv")
DESIGN3_FILE = os.path.join(BASE_DIR, "ANU022_DEP_Design3.csv")  # kept only for reference

OUT_FIG3C_PDF = os.path.join(BASE_DIR, "Fig3C_Proteins_per_sample_REORDERED_filtered.pdf")
OUT_FIG3C_PNG = os.path.join(BASE_DIR, "Fig3C_Proteins_per_sample_REORDERED_filtered.png")

OUT_S8A_PDF   = os.path.join(BASE_DIR, "FigS8A_Proteins_per_sample_REORDERED_filtered.pdf")
OUT_S8A_PNG   = os.path.join(BASE_DIR, "FigS8A_Proteins_per_sample_REORDERED_filtered.png")

OUT_MAIN_PCA_LEGEND_PDF = os.path.join(BASE_DIR, "Fig3D_PCA_legend_main.pdf")
OUT_MAIN_PCA_LEGEND_PNG = os.path.join(BASE_DIR, "Fig3D_PCA_legend_main.png")

OUT_SUPP_PCA_LEGEND_PDF = os.path.join(BASE_DIR, "FigS8B_PCA_legend_supp.pdf")
OUT_SUPP_PCA_LEGEND_PNG = os.path.join(BASE_DIR, "FigS8B_PCA_legend_supp.png")

# ==========================================================
# STYLE
# ==========================================================

def set_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 18,
        "axes.titlesize": 22,
        "axes.labelsize": 20,
        "xtick.labelsize": 15,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

set_style()

# ==========================================================
# COLOURS
# ==========================================================

workflow_colors = {
    "Control": "#000000",
    "FASP_8M": "#E9DC6D",
    "FASP_Pellet": "#F4A637",
    "FASP_SDS": "#DB5829",
    "Gel": "#894B45",
    "OrganControl": "#00C992",
    "OrganFASP_8M": "#008A69",
    "OrganPellet": "#386350",
    "OrganPelletpellet": "#D2BBD7",
    "OrganSN": "#AE75A2",
    "Pellet": "#B6DBFF",
    "Pelletpellet": "#7BB0DF",
    "SN": "#1964B0",
    "TRPA1_gel": "#DEDEDE",
    "TRPA1_Xeno": "#C6C6C6",
    "TRPV1_gel": "#8C8C8C",
    "TRPV1_Xeno": "#4D4D4D",
}

# ==========================================================
# SAMPLE ORDERS
# ==========================================================

FIG3C_SAMPLE_ORDER = [
    "SN_1", "SN_2", "SN_3",
    "Control_1", "Control_2", "Control_3",
    "Pellet_1", "Pellet_2", "Pellet_3",
    "Pelletpellet_1", "Pelletpellet_2", "Pelletpellet_3",
    "FASP_8M_1", "FASP_8M_2", "FASP_8M_3",
    "FASP_SDS_1", "FASP_SDS_2", "FASP_SDS_3",
    "FASP_Pellet_1", "FASP_Pellet_2", "FASP_Pellet_3",
    "Gel_1", "Gel_2", "Gel_3",
]

S8A_SAMPLE_ORDER = [
    "SN_1", "SN_2", "SN_3",
    "Control_1", "Control_2", "Control_3",
    "Pellet_1", "Pellet_2", "Pellet_3",
    "Pelletpellet_1", "Pelletpellet_2", "Pelletpellet_3",
    "FASP_8M_1", "FASP_8M_2", "FASP_8M_3",
    "FASP_SDS_1", "FASP_SDS_2", "FASP_SDS_3",
    "FASP_Pellet_1", "FASP_Pellet_2", "FASP_Pellet_3",
    "Gel_1", "Gel_2", "Gel_3",
    "DRG_SN", "DRG_Sample", "DRG_Pellet", "DRG_Pelletpellet", "DRG8M",
    "Kidney_SN", "Kidney_Sample", "Kidney_Pellet", "Kidney_Pelletpellet",
    "Testes_SN", "Testes_Sample", "Testes_Pellet", "Testes_Pelletpellet",
    "Lung8M", "Skin8M",
    "TRPA1_gel", "TRPA1_Xeno", "TRPV1_gel", "TRPV1_Xeno",
]

MAIN_PCA_ORDER = [
    "SN",
    "Control",
    "Pellet",
    "Pelletpellet",
    "FASP_8M",
    "FASP_SDS",
    "FASP_Pellet",
    "Gel",
]

SUPP_PCA_ORDER = [
    "SN",
    "Control",
    "Pellet",
    "Pelletpellet",
    "FASP_8M",
    "FASP_SDS",
    "FASP_Pellet",
    "Gel",
    "OrganSN",
    "OrganControl",
    "OrganPellet",
    "OrganPelletpellet",
    "OrganFASP_8M",
    "TRPA1_gel",
    "TRPA1_Xeno",
    "TRPV1_gel",
    "TRPV1_Xeno",
]

# ==========================================================
# DISPLAY LABELS FOR X-AXIS
# Keep internal sample names unchanged; only improve figure labels
# ==========================================================

DISPLAY_LABELS = {
    # Cortex
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

    # Organ samples - changed from Sample to Control for display
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

    # Positive controls
    "TRPA1_gel": "TRPA1_gel",
    "TRPA1_Xeno": "TRPA1_Xeno",
    "TRPV1_gel": "TRPV1_gel",
    "TRPV1_Xeno": "TRPV1_Xeno",
}

# ==========================================================
# HELPER FUNCTIONS
# ==========================================================

def clean_label(x):
    x = str(x)
    x = re.sub(r"^LFQ Intensity\s+", "", x)
    x = re.sub(r"\.raw$", "", x)
    return x.strip()

def dep_filter_from_design(prot_df: pd.DataFrame, design_df: pd.DataFrame) -> pd.DataFrame:
    """
    Mimic DEP::filter_missval(data_se, thr = 0):
    keep proteins present in all samples of at least one condition.
    """
    keep_mask = pd.Series(False, index=prot_df.index)

    for cond in design_df["condition"].dropna().unique():
        samples = design_df.loc[design_df["condition"] == cond, "Sample"].tolist()
        samples = [s for s in samples if s in prot_df.columns]
        if len(samples) == 0:
            continue
        keep_mask = keep_mask | prot_df[samples].notna().all(axis=1)

    return prot_df.loc[keep_mask].copy()

def print_section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)

def plot_bar_from_order(
    counts_series: pd.Series,
    sample_order: list,
    condition_map: dict,
    out_pdf: str,
    out_png: str,
    display_labels: dict = None,
    ylim_top: int = 6000,
    fig_width: float = 10.0,
    fig_height: float = 6.0,
    bar_width: float = 0.82,
):
    samples_present = [s for s in sample_order if s in counts_series.index]

    missing_in_condition_map = [s for s in samples_present if s not in condition_map]
    if missing_in_condition_map:
        raise ValueError(
            f"These samples are present in counts_series but missing from condition_map: "
            f"{missing_in_condition_map}"
        )

    values = [counts_series[s] for s in samples_present]
    colors = [workflow_colors[condition_map[s]] for s in samples_present]

    if display_labels is None:
        xtick_labels = samples_present
    else:
        xtick_labels = [display_labels.get(s, s) for s in samples_present]

    x = np.arange(len(samples_present))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.bar(x, values, width=bar_width, color=colors, edgecolor="none")

    ax.set_title("Proteins per sample", fontweight="bold")
    ax.set_ylabel("Number of proteins")
    ax.set_xlabel("")
    ax.set_ylim(0, ylim_top)

    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, rotation=90)

    ax.set_xlim(-0.5, len(samples_present) - 0.5)
    ax.margins(x=0)

    ax.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close()

def make_pca_legend_only(
    out_pdf: str,
    out_png: str,
    legend_order: list,
    include_replicate: bool = True,
    figsize=(3.4, 5.2)
):
    missing_colors = [c for c in legend_order if c not in workflow_colors]
    if missing_colors:
        raise ValueError(f"Missing colours for legend conditions: {missing_colors}")

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    if include_replicate:
        rep_handles = [
            Line2D([0], [0], marker='o', linestyle='None',
                   markerfacecolor='black', markeredgecolor='black',
                   markersize=8, label='1'),
            Line2D([0], [0], marker='o', linestyle='None',
                   markerfacecolor='black', markeredgecolor='black',
                   markersize=8, label='2'),
            Line2D([0], [0], marker='o', linestyle='None',
                   markerfacecolor='black', markeredgecolor='black',
                   markersize=8, label='3'),
        ]

        leg1 = ax.legend(
            handles=rep_handles,
            title="replicate",
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            frameon=False,
            handlelength=1.0,
            handletextpad=0.8,
            borderaxespad=0.0,
            labelspacing=0.6
        )
        leg1.get_title().set_fontsize(18)
        for txt in leg1.get_texts():
            txt.set_fontsize(16)
        ax.add_artist(leg1)
        anchor_y = 0.56
    else:
        anchor_y = 1.0

    cond_handles = [
        Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=workflow_colors[c], markeredgecolor='none',
               markersize=8, label=c)
        for c in legend_order
    ]

    leg2 = ax.legend(
        handles=cond_handles,
        title="condition",
        loc="upper left",
        bbox_to_anchor=(0.0, anchor_y),
        frameon=False,
        handlelength=1.0,
        handletextpad=0.8,
        borderaxespad=0.0,
        labelspacing=0.6
    )
    leg2.get_title().set_fontsize(18)
    for txt in leg2.get_texts():
        txt.set_fontsize(16)

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight", transparent=True)
    plt.savefig(out_png, dpi=600, bbox_inches="tight", transparent=True)
    plt.close()

# ==========================================================
# LOAD DATA
# ==========================================================

design2 = pd.read_csv(DESIGN2_FILE)
prot = pd.read_csv(PROT_FILE)

design2["Sample"] = design2["label"].apply(clean_label)
design2["condition"] = design2["condition"].astype(str).str.strip()

# Remove contaminants
if "Protein.Group" in prot.columns:
    prot = prot[~prot["Protein.Group"].astype(str).str.contains("cRAP", na=False)].copy()

# LFQ matrix
intensity_cols_raw = [c for c in prot.columns if str(c).startswith("LFQ Intensity ")]
sample_map = {c: clean_label(c) for c in intensity_cols_raw}

prot_lfq = prot[intensity_cols_raw].copy()
prot_lfq.columns = [sample_map[c] for c in intensity_cols_raw]
prot_lfq = prot_lfq.apply(pd.to_numeric, errors="coerce")

# Keep only samples present in Design2
design2_samples = design2["Sample"].tolist()
present_samples = [s for s in design2_samples if s in prot_lfq.columns]
prot_lfq = prot_lfq[present_samples].copy()

design2_present = design2[design2["Sample"].isin(present_samples)].copy()

# DEP-like filter
prot_filt = dep_filter_from_design(prot_lfq, design2_present)

# Counts per sample
protein_counts = prot_filt.notna().sum(axis=0)

# Sample -> condition map
condition_map = dict(zip(design2_present["Sample"], design2_present["condition"]))

# ==========================================================
# DIAGNOSTIC PRINTS (just to double-check)
# ==========================================================

print_section("Samples in Design2 and protein table")
for s in present_samples:
    print(f"{s}\t{condition_map[s]}")

print_section("Requested Fig 3C samples missing from filtered input")
missing_fig3c = [s for s in FIG3C_SAMPLE_ORDER if s not in protein_counts.index]
if missing_fig3c:
    for s in missing_fig3c:
        print(s)
else:
    print("None")

print_section("Requested S8A samples missing from filtered input")
missing_s8a = [s for s in S8A_SAMPLE_ORDER if s not in protein_counts.index]
if missing_s8a:
    for s in missing_s8a:
        print(s)
else:
    print("None")

print_section("Conditions present in filtered data")
present_conditions = sorted(set(condition_map.values()))
for c in present_conditions:
    print(c)

print_section("Conditions missing from MAIN_PCA_ORDER")
missing_main_pca = sorted(set(present_conditions) - set(MAIN_PCA_ORDER))
if missing_main_pca:
    for c in missing_main_pca:
        print(c)
else:
    print("None")

print_section("Conditions missing from SUPP_PCA_ORDER")
missing_supp_pca = sorted(set(present_conditions) - set(SUPP_PCA_ORDER))
if missing_supp_pca:
    for c in missing_supp_pca:
        print(c)
else:
    print("None")

print_section("Conditions defined in workflow_colors but not present in filtered data")
unused_colors = sorted(set(workflow_colors.keys()) - set(present_conditions))
if unused_colors:
    for c in unused_colors:
        print(c)
else:
    print("None")

print_section("Conditions present in filtered data but missing from workflow_colors")
missing_colour_conditions = sorted(set(present_conditions) - set(workflow_colors.keys()))
if missing_colour_conditions:
    for c in missing_colour_conditions:
        print(c)
    sys.exit("Stop: add missing colours before plotting.")
else:
    print("None")

print_section("Display labels used for organ controls")
for raw_name in ["DRG_Sample", "Kidney_Sample", "Testes_Sample"]:
    if raw_name in DISPLAY_LABELS:
        print(f"{raw_name} -> {DISPLAY_LABELS[raw_name]}")

# ==========================================================
# FIG 3C
# ==========================================================

fig3c_condition_map = {s: condition_map[s] for s in FIG3C_SAMPLE_ORDER if s in condition_map}

plot_bar_from_order(
    counts_series=protein_counts,
    sample_order=FIG3C_SAMPLE_ORDER,
    condition_map=fig3c_condition_map,
    out_pdf=OUT_FIG3C_PDF,
    out_png=OUT_FIG3C_PNG,
    display_labels=DISPLAY_LABELS,
    ylim_top=6000,
    fig_width=9.0,
    fig_height=6.0,
    bar_width=0.82
)

# ==========================================================
# FIG S8A
# ==========================================================

s8a_condition_map = {s: condition_map[s] for s in S8A_SAMPLE_ORDER if s in condition_map}

plot_bar_from_order(
    counts_series=protein_counts,
    sample_order=S8A_SAMPLE_ORDER,
    condition_map=s8a_condition_map,
    out_pdf=OUT_S8A_PDF,
    out_png=OUT_S8A_PNG,
    display_labels=DISPLAY_LABELS,
    ylim_top=6000,
    fig_width=10.8,
    fig_height=6.0,
    bar_width=0.82
)

# ==========================================================
# PCA LEGENDS
# ==========================================================

make_pca_legend_only(
    OUT_MAIN_PCA_LEGEND_PDF,
    OUT_MAIN_PCA_LEGEND_PNG,
    legend_order=MAIN_PCA_ORDER,
    include_replicate=True,
    figsize=(3.2, 4.6)
)

make_pca_legend_only(
    OUT_SUPP_PCA_LEGEND_PDF,
    OUT_SUPP_PCA_LEGEND_PNG,
    legend_order=SUPP_PCA_ORDER,
    include_replicate=True,
    figsize=(3.6, 7.2)
)

print_section("Done")
print(OUT_FIG3C_PDF)
print(OUT_FIG3C_PNG)
print(OUT_S8A_PDF)
print(OUT_S8A_PNG)
print(OUT_MAIN_PCA_LEGEND_PDF)
print(OUT_MAIN_PCA_LEGEND_PNG)
print(OUT_SUPP_PCA_LEGEND_PDF)
print(OUT_SUPP_PCA_LEGEND_PNG)