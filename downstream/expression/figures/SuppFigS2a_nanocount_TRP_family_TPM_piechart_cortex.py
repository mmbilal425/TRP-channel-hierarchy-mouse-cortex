import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

# -------- font settings (same as Illumina plot) --------
available_fonts = {f.name for f in fm.fontManager.ttflist}
font_choice = "Arial" if "Arial" in available_fonts else "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_choice,
    "font.size": 14,
    "pdf.fonttype": 42,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black"
})

# -------- paths --------
tpm_file = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results/TRP_TPM_CortexMean_NanoCount_ordered.tsv"
out_dir  = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/nanocount/results"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "TRP_TPM_CortexMean_NanoCount_family_TPM_piechart.pdf")
out_tsv = os.path.join(out_dir, "TRP_TPM_CortexMean_NanoCount_family_TPM_totals.tsv")

# -------- load TPM table --------
df = pd.read_csv(tpm_file, sep="\t")

# Expecting columns: GeneName, family, TPM_rep1/2/3, TPM_mean
df = df[["GeneName", "family", "TPM_mean"]].copy()
df["GeneName"] = df["GeneName"].astype(str)
df["family"]   = df["family"].astype(str)
df["TPM_mean"] = pd.to_numeric(df["TPM_mean"], errors="coerce").fillna(0.0)

# -------- aggregate TPM by family --------
family_tpm = (
    df.groupby("family", as_index=False)["TPM_mean"]
      .sum()
      .rename(columns={"TPM_mean": "TPM"})
      .sort_values("TPM", ascending=False)
)

# save table (optional but useful for supplement)
family_tpm.to_csv(out_tsv, sep="\t", index=False)

# -------- colour palette (same as Illumina pie) --------
palette = {
    "TRPML": "#4682b4",  # Mcoln
    "TRPP":  "#dda0dd",  # Pkd
    "TRPA":  "#696969",
    "TRPC":  "#fc8d62",
    "TRPM":  "#66a61e",
    "TRPV":  "#8da0cb",
}
colors = [palette.get(f, "#999999") for f in family_tpm["family"]]

# -------- TPM-based pie chart --------
plt.figure(figsize=(7, 7))
wedges, texts, autotexts = plt.pie(
    family_tpm["TPM"].values,
    labels=family_tpm["family"].values,
    autopct="%1.1f%%",
    startangle=90,
    counterclock=False,
    colors=colors,
    textprops={"fontsize": 12},
)

# Optional: nudge TRPA label up if present (tiny slice)
if "TRPA" in list(family_tpm["family"]):
    idx = list(family_tpm["family"]).index("TRPA")
    x, y = autotexts[idx].get_position()
    autotexts[idx].set_position((x, y + 0.12))

plt.title("TPM distribution by TRP gene family (Nanopore, cortex)", fontsize=18)
plt.tight_layout()

plt.savefig(out_pdf, dpi=600, format="pdf", bbox_inches="tight")
plt.show()

print("Saved NanoCount TPM-based family pie chart to:", out_pdf)
print("Saved NanoCount family TPM table to:", out_tsv)
