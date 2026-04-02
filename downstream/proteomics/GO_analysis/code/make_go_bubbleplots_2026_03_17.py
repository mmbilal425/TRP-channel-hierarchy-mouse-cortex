#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# =========================
# User defaults
# =========================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BASE_DIR = SCRIPT_DIR

# Reordered
DEFAULT_PROTOCOLS = [
    "SN_GO",
    "Control_GO",
    "Pellet_GO",
    "Pelletpellet_GO",
    "FASP_Pellet_GO",
    "FASP_8M_GO",
    "FASP_SDS_GO",
    "Gel_GO",
]

# Natasha's GO colours
GO_LOW = "#7BB0DF"
GO_HIGH = "#DB5829"

# Natasha's blue heatmap colours (stored here for consistency)
HEATMAP_LOW = "#B6DBFF"
HEATMAP_HIGH = "#1964B0"


def wrap(s: str, width: int = 52) -> str:
    return "\n".join(textwrap.wrap(str(s), width=width))


def protocol_name_from_key(key: str) -> str:
    name = re.sub(r"_GO.*$", "", key)
    name = name.replace("Control", "Control/membrane")
    name = name.replace("Pelletpellet", "Pellet-pellet")
    name = name.replace("FASP_Pellet", "FASP-Pellet")
    name = name.replace("FASP_8M_PNG", "FASP 8M+PNGaseF")
    name = name.replace("FASP_SDS_PNG", "FASP SDS+PNGaseF")
    name = name.replace("FASP_8M", "FASP 8M urea")
    name = name.replace("FASP_SDS", "FASP SDS")
    return name


def read_go_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]

    if "term_name" not in df.columns and "Description" in df.columns:
        df["term_name"] = df["Description"]

    if "negative_log10_of_adjusted_p_value" not in df.columns and "adjusted_p_value" in df.columns:
        adj = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
        adj = adj.clip(lower=np.finfo(float).tiny)
        df["negative_log10_of_adjusted_p_value"] = -np.log10(adj)

    if "intersection_size" not in df.columns and "Count" in df.columns:
        df["intersection_size"] = pd.to_numeric(df["Count"], errors="coerce")

    for col in [
        "adjusted_p_value",
        "negative_log10_of_adjusted_p_value",
        "intersection_size",
        "term_size",
        "query_size",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def set_plot_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 20,
        "axes.titlesize": 24,
        "axes.labelsize": 22,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 15,
        "legend.title_fontsize": 16,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
    })


def bubbleplot_gridspec(
    df: pd.DataFrame,
    protocols: list[str],
    terms: list[str],
    out_pdf: str,
    title: str,
    xlabel: str,
    ylabel: str,
    size_scale: float,
    size_legend_values,
    fig_width: float,
    fig_height: float,
    y_wrap: int = 56,
) -> None:
    set_plot_style()

    x_map = {p: i for i, p in enumerate(protocols)}
    y_map = {t: i for i, t in enumerate(terms)}

    plot_df = df.copy()
    plot_df["x"] = plot_df["protocol"].map(x_map)
    plot_df["y"] = plot_df["term_name"].map(y_map)
    plot_df = plot_df.dropna(subset=["x", "y"]).copy()

    plot_df["size"] = plot_df["intersection_size"].astype(float)
    plot_df["color"] = plot_df["negative_log10_of_adjusted_p_value"].astype(float)

    custom_cmap = LinearSegmentedColormap.from_list(
        "natasha_go_gradient",
        [GO_LOW, GO_HIGH]
    )

    fig = plt.figure(figsize=(fig_width, fig_height))

    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[1.0, 0.08],
        wspace=0.06
    )

    ax = fig.add_subplot(gs[0, 0])

    gs_right = gs[0, 1].subgridspec(
        nrows=5,
        ncols=1,
        height_ratios=[0.14, 0.18, 0.08, 0.48, 0.12],
        hspace=0.02
    )
    top_spacer = fig.add_subplot(gs_right[0, 0])
    cax = fig.add_subplot(gs_right[1, 0])
    mid_spacer = fig.add_subplot(gs_right[2, 0])
    lax = fig.add_subplot(gs_right[3, 0])
    bottom_spacer = fig.add_subplot(gs_right[4, 0])

    for extra_ax in [top_spacer, mid_spacer, bottom_spacer]:
        extra_ax.axis("off")

    sc = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["size"] * size_scale,
        c=plot_df["color"],
        cmap=custom_cmap,
        alpha=0.95,
        edgecolors="black",
        linewidths=0.28
    )

    ax.set_xticks(range(len(protocols)))
    ax.set_xticklabels(protocols, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([wrap(t, y_wrap) for t in terms])

    ax.set_xlabel(xlabel, labelpad=16)
    ax.set_ylabel(ylabel, labelpad=18)
    ax.set_title(title, pad=18)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    fig.subplots_adjust(
        left=0.38,
        right=0.935,
        top=0.92,
        bottom=0.19
    )

    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label(r"$-\log_{10}$(adjusted p-value)", fontsize=13, labelpad=5)
    cbar.ax.tick_params(labelsize=11, width=0.8, length=3)

    lax.axis("off")
    handles, labels = [], []
    for s in size_legend_values:
        h = ax.scatter(
            [],
            [],
            s=s * size_scale,
            c="gray",
            alpha=0.75,
            edgecolors="black",
            linewidths=0.28
        )
        handles.append(h)
        labels.append(str(s))

    lax.legend(
        handles,
        labels,
        title="Intersection size",
        loc="upper left",
        frameon=False,
        borderaxespad=0.0,
        labelspacing=0.8,
        handletextpad=0.8,
        scatterpoints=1,
        markerscale=0.85
    )

    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Make GO:CC bubble plots from precomputed GO CSV results."
    )
    parser.add_argument(
        "--in_dir",
        default=DEFAULT_BASE_DIR,
        help="Folder containing *_GO*.csv files"
    )
    parser.add_argument(
        "--out_dir",
        default=DEFAULT_BASE_DIR,
        help="Overwrite outputs in the same folder by default"
    )
    parser.add_argument(
        "--protocols",
        nargs="*",
        default=DEFAULT_PROTOCOLS,
        help="Basenames (no .csv) of protocol GO files"
    )
    parser.add_argument(
        "--topn_per_protocol",
        type=int,
        default=12,
        help="Number of top membrane-associated GO:CC terms retained per protocol"
    )
    args = parser.parse_args()

    in_dir = args.in_dir
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    csv_paths = glob.glob(os.path.join(in_dir, "*_GO*.csv"))
    dfs = {}
    for p in csv_paths:
        key = os.path.splitext(os.path.basename(p))[0]
        dfs[key] = read_go_csv(p)

    missing = [k for k in args.protocols if k not in dfs]
    if missing:
        raise SystemExit(
            f"Missing protocol CSV(s): {missing}\nAvailable files: {sorted(dfs.keys())}"
        )

    prot_dfs = []
    for k in args.protocols:
        df = dfs[k].copy()
        if "source" in df.columns:
            df = df[df["source"].astype(str).str.contains("GO:CC", na=False)].copy()
        df["protocol"] = protocol_name_from_key(k)
        prot_dfs.append(df)

    all_cc = pd.concat(prot_dfs, ignore_index=True)

    mem_mask = all_cc["term_name"].astype(str).str.contains(
        r"membrane|vesicle|synap|endoplasmic|mitochond|axon|dendrit|cell projection|raft|microdomain|\bER\b",
        flags=re.IGNORECASE,
        regex=True,
        na=False,
    )
    mem_cc = all_cc[mem_mask].copy()

    top_terms = set()
    for prot, dfp in mem_cc.groupby("protocol"):
        dfp = dfp.sort_values("negative_log10_of_adjusted_p_value", ascending=False).head(args.topn_per_protocol)
        top_terms |= set(dfp["term_name"].astype(str))

    broad_df = mem_cc[mem_cc["term_name"].isin(top_terms)].copy()
    broad_df = broad_df.sort_values("adjusted_p_value").drop_duplicates(
        ["protocol", "term_name"], keep="first"
    )

    term_order = (
        broad_df.groupby("term_name")["negative_log10_of_adjusted_p_value"]
        .max()
        .sort_values(ascending=True)
        .index.tolist()
    )
    protocol_order = [protocol_name_from_key(k) for k in args.protocols]

    broad_pdf = os.path.join(out_dir, "GO_CC_membrane_bubbleplot.pdf")
    bubbleplot_gridspec(
        broad_df,
        protocol_order,
        term_order,
        broad_pdf,
        "Membrane-associated GO:CC enrichment across extraction workflows",
        "Extraction / digestion protocol",
        "GO:CC (membrane / vesicle-associated terms)",
        size_scale=22.0,
        size_legend_values=(4, 8, 16, 24),
        fig_width=19.5,
        fig_height=max(11.0, len(term_order) * 0.42),
        y_wrap=58,
    )

    broad_tidy = broad_df[
        [
            "protocol",
            "term_name",
            "term_id",
            "adjusted_p_value",
            "negative_log10_of_adjusted_p_value",
            "intersection_size",
            "term_size",
            "query_size",
        ]
    ].copy()
    broad_tidy.to_csv(
        os.path.join(out_dir, "GO_CC_membrane_bubbleplot_tidy.csv"),
        index=False
    )

    focused_key = "SN_GO_CC_Membrane"
    if focused_key in dfs:
        focus_src = dfs[focused_key].copy()

        if "Description" not in focus_src.columns and "term_name" in focus_src.columns:
            focus_src["Description"] = focus_src["term_name"]

        if "term_id" in focus_src.columns and "Description" in focus_src.columns:
            target_term_ids = set(focus_src["term_id"].astype(str))
            term_map = (
                focus_src[["term_id", "Description"]]
                .drop_duplicates()
                .set_index("term_id")["Description"]
                .to_dict()
            )

            rows = []
            for k in args.protocols:
                df = dfs[k].copy()
                if "source" in df.columns:
                    df = df[df["source"].astype(str).str.contains("GO:CC", na=False)].copy()
                df = df[df["term_id"].astype(str).isin(target_term_ids)].copy()
                df["protocol"] = protocol_name_from_key(k)
                df["term_name"] = df["term_id"].astype(str).map(term_map)
                rows.append(df)

            focused_df = pd.concat(rows, ignore_index=True)
            focused_df = focused_df.sort_values("adjusted_p_value").drop_duplicates(
                ["protocol", "term_name"], keep="first"
            )

            term_order2 = (
                focused_df.groupby("term_name")["negative_log10_of_adjusted_p_value"]
                .max()
                .sort_values(ascending=True)
                .index.tolist()
            )

            focused_pdf = os.path.join(out_dir, "GO_CC_synaptic_membrane_terms_bubbleplot.pdf")
            bubbleplot_gridspec(
                focused_df,
                protocol_order,
                term_order2,
                focused_pdf,
                "Synaptic / vesicular membrane GO:CC terms across workflows",
                "Extraction / digestion protocol",
                "GO:CC membrane-associated terms",
                size_scale=25.0,
                size_legend_values=(2, 4, 8, 12),
                fig_width=19.0,
                fig_height=max(10.5, len(term_order2) * 0.55),
                y_wrap=58,
            )

            focused_tidy = focused_df[
                [
                    "protocol",
                    "term_name",
                    "term_id",
                    "adjusted_p_value",
                    "negative_log10_of_adjusted_p_value",
                    "intersection_size",
                    "term_size",
                    "query_size",
                ]
            ].copy()
            focused_tidy.to_csv(
                os.path.join(out_dir, "GO_CC_synaptic_membrane_terms_tidy.csv"),
                index=False
            )

    report_path = os.path.join(out_dir, "run_summary.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("GO bubble plot run summary\n")
        f.write("==========================\n")
        f.write(f"Input directory: {in_dir}\n")
        f.write(f"Output directory: {out_dir}\n")
        f.write(f"Protocols: {', '.join(args.protocols)}\n")
        f.write(f"Top N per protocol: {args.topn_per_protocol}\n")
        f.write(f"GO colours: low={GO_LOW}, high={GO_HIGH}\n")
        f.write(f"Heatmap colours stored for reuse: low={HEATMAP_LOW}, high={HEATMAP_HIGH}\n")

    print("Done.")
    print(f"Output directory: {out_dir}")
    print(f"Wrote: {broad_pdf}")
    print(f"Wrote: {os.path.join(out_dir, 'GO_CC_membrane_bubbleplot_tidy.csv')}")
    if focused_key in dfs:
        print(f"Wrote: {os.path.join(out_dir, 'GO_CC_synaptic_membrane_terms_bubbleplot.pdf')}")
        print(f"Wrote: {os.path.join(out_dir, 'GO_CC_synaptic_membrane_terms_tidy.csv')}")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    main()