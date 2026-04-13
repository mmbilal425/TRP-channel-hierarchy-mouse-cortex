#!/usr/bin/env python3
"""
Supplementary Figure S5a
m6A stoichiometry distribution across TRP-associated mRNA sites
(Nanopore direct RNA; modkit)

Input:
- result_df (loaded upstream from modkit pileup parquet)

Output:
- stoichiometry_histogram.tsv
- stoichiometry_histogram.pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import logging
from matplotlib.ticker import FuncFormatter

# Silence fontTools warnings
logging.getLogger("fontTools").setLevel(logging.ERROR)

# ------------------------------------------------------------------
# Sanity checks (explicit for reproducibility)
# ------------------------------------------------------------------
assert "result_df" in globals(), (
    "result_df not found. Load modkit pileup table before running."
)
assert "plot_dir" in globals(), (
    "plot_dir not defined. Set output directory before running."
)

# ------------------------------------------------------------------
# Filter: m6A ('a') + modified sites only
# ------------------------------------------------------------------
hist_df = result_df[
    (result_df["mod_type"] == "a") &
    (result_df["is_modified"] == 1)
].copy()

print(f"[INFO] Number of m6A-modified sites: {len(hist_df)}")

# Save filtered table (for traceability)
summary_file = os.path.join(plot_dir, "stoichiometry_histogram.tsv")
hist_df.to_csv(summary_file, sep="\t", index=False)
print(f"[INFO] Filtered data written to: {summary_file}")

# ------------------------------------------------------------------
# Plot configuration (Nature-style)
# ------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "pdf.fonttype": 42,
    "xtick.major.pad": 1.4,
    "ytick.major.pad": 1.4,
})

fig = plt.figure(figsize=(2.3 / 2.54, 2.7 / 2.54))
fig.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.92)
ax = fig.add_subplot(111)

# Histogram
n, bins, _ = ax.hist(
    hist_df["merged_bam_stoich"].to_numpy(dtype=float),
    bins=20,
    range=(0, 100),
    color="blue",
    edgecolor="none",
    align="left"
)

# Axis labels
ax.set_xlabel("m6A stoichiometry (%)", fontsize=6, labelpad=1.4)
ax.set_ylabel("Number of sites (×10)", fontsize=6, labelpad=1.4)

# Spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)

# Ticks
tick_length = 0.06 * 28.35  # mm → points
ax.tick_params(
    axis="both",
    which="major",
    length=tick_length,
    width=0.5,
    direction="out",
    pad=1.4
)

ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")

# Y-axis formatter (scaled to ×10)
def y_formatter(x, pos):
    return f"{int(x / 1e4)}"

ax.yaxis.set_major_formatter(FuncFormatter(y_formatter))
ax.set_yticks(np.arange(0, 5) * 1e4)
ax.set_ylim(0, 4 * 1e4)

# X-axis (fixed, no padding)
ax.set_xticks([0, 50, 100])
ax.set_xticklabels(["0", "50", "100"], fontsize=6)
ax.set_xlim(-0.01, 100)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
out_pdf = os.path.join(plot_dir, "stoichiometry_histogram.pdf")
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"[OK] Figure saved to: {out_pdf}")
