import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from pathlib import Path
import re

# ---------------- FONT SETTINGS ----------------
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

# ---------------- PATHS ----------------
outdir = Path("/g/data/lf10/mb1232/illumina_data/mouse_brain_illumina/salmon_output_r115/results")
outdir.mkdir(exist_ok=True, parents=True)

merged_tpm_file = outdir / "Merged_Salmon_GeneSymbol_TPMs_r115.tsv"
trp_list_file   = Path("/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data/trp_gene_ids.txt")

TPM_CUTOFF = 0.5

# ---------------- TRP LIST LOADER (robust) ----------------
def load_trp_symbols(path: Path):
    syms = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("geneid"):
                continue
            # handles "ENSMUSG...→Trpv1" or whitespace/tab separated
            if "→" in line:
                sym = line.split("→", 1)[1].strip()
            else:
                sym = re.split(r"[,\t ]+", line)[-1].strip()
            sym = re.split(r"[,\t ]+", sym)[0]
            if sym:
                syms.append(sym)
    return syms

trp_syms = load_trp_symbols(trp_list_file)
TRP_UP = {s.upper() for s in trp_syms}
if not TRP_UP:
    raise ValueError("TRP list loaded empty. Check trp_gene_ids.txt formatting/path.")

# ---------------- LOAD MERGED TPM TABLE ----------------
df = pd.read_csv(merged_tpm_file, sep="\t")

required = {"gene_name", "TPM_Cortex_mean", "TPM_DRG"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {merged_tpm_file}:\n  {sorted(missing)}")

df["GeneName"] = df["gene_name"].astype(str)
df["TPM_Cortex_mean"] = pd.to_numeric(df["TPM_Cortex_mean"], errors="coerce").fillna(0.0)
df["TPM_DRG"]         = pd.to_numeric(df["TPM_DRG"],         errors="coerce").fillna(0.0)

# ---------------- FILTER TO TRP ONLY (THIS WAS THE MISSING STEP) ----------------
df = df[df["GeneName"].str.upper().isin(TRP_UP)].copy()

# ---------------- DETERMINE EXPRESSED TRP GENES ----------------
cortex_expr = set(df.loc[df["TPM_Cortex_mean"] >= TPM_CUTOFF, "GeneName"])
drg_expr    = set(df.loc[df["TPM_DRG"]         >= TPM_CUTOFF, "GeneName"])

only_cortex = sorted(list(cortex_expr - drg_expr))
only_drg    = sorted(list(drg_expr - cortex_expr))
both        = sorted(list(cortex_expr & drg_expr))

# ---------------- SAVE TSV OUTPUTS ----------------
cut_str = str(TPM_CUTOFF).replace(".", "p")
prefix = f"TRP_VENN_TPMge{cut_str}"

(outdir / f"{prefix}_CortexOnly.tsv").write_text(
    "GeneName\n" + "\n".join(only_cortex) + "\n"
)
(outdir / f"{prefix}_DRGOnly.tsv").write_text(
    "GeneName\n" + "\n".join(only_drg) + "\n"
)
(outdir / f"{prefix}_Both.tsv").write_text(
    "GeneName\n" + "\n".join(both) + "\n"
)

print("Counts (TRP only):")
print(" Cortex only:", len(only_cortex))
print(" DRG only:   ", len(only_drg))
print(" Both:       ", len(both))

# ---------------- PLOT (Venn-style circles, no external libs) ----------------
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.set_aspect("equal")
ax.axis("off")

teal   = "#1b9e77"   # Cortex
purple = "#7570b3"   # DRG

r = 1.6
x1, y1 = -1.0, 0.0
x2, y2 =  1.0, 0.0

c1 = plt.Circle((x1, y1), r, color=teal,   alpha=0.35, lw=1.5, ec=teal)
c2 = plt.Circle((x2, y2), r, color=purple, alpha=0.35, lw=1.5, ec=purple)
ax.add_patch(c1)
ax.add_patch(c2)

ax.text(x1, y1 + r + 0.35, "Cortex", color=teal,   ha="center", va="bottom", fontsize=13)
ax.text(x2, y2 + r + 0.35, "DRG",    color=purple, ha="center", va="bottom", fontsize=13)

ax.text(x1 - 0.55, 0.0, str(len(only_cortex)), fontsize=13, ha="center", va="center")
ax.text(0.0,       0.0, str(len(both)),        fontsize=13, ha="center", va="center")
ax.text(x2 + 0.55, 0.0, str(len(only_drg)),    fontsize=13, ha="center", va="center")

ax.set_xlim(-3.2, 3.2)
ax.set_ylim(-2.1, 2.4)

plt.tight_layout()

out_pdf = outdir / f"{prefix}_Venn_TRP.pdf"
plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
plt.show()

print("\nSaved files:")
print(" ", outdir / f"{prefix}_CortexOnly.tsv")
print(" ", outdir / f"{prefix}_DRGOnly.tsv")
print(" ", outdir / f"{prefix}_Both.tsv")
print(" ", out_pdf)
