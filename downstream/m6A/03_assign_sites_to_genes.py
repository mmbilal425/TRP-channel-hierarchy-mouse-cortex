#!/usr/bin/env python3
"""
03_assign_sites_to_genes.py

Assigns each site to a gene (gene_id + biotype) by overlap with gene intervals.

Input: modified_sites_table.parquet (from 02)  OR a user-specified parquet/tsv
Output: modified_sites_with_genes.parquet (+ tsv)
"""

import os
import logging
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Tuple
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import tqdm

# ----------------------------
# CONFIG (your defaults)
# ----------------------------
PLOTS_DIR = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/plots"

IN_PARQUET = "modified_sites_table.parquet"
OUT_PARQUET = "modified_sites_with_genes.parquet"
OUT_TSV = "modified_sites_with_genes.tsv"

GENE_BED = "/g/data/qq78/as7425/rna_biogenesis_maps/dFORCE_prod_dec24/mouse/second_pass_index/annotation_with_clusters/updated_gene.bed"
GTF = "/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.115.gtf"

N_CORES = 32
CHUNK_SIZE = 4_590_000

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")


@dataclass
class ChromosomeData:
    starts: np.ndarray
    ends: np.ndarray
    gene_ids: np.ndarray
    biotypes: np.ndarray


def parse_gtf_gene_biotypes(gtf_path: str) -> pd.DataFrame:
    gtf_df = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=["chrom","source","feature","start","end","score","strand","frame","attribute"],
        dtype=str
    )
    gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()

    def parse_attr(s):
        gene_id = None
        biotype = None
        for attr in s.strip().split(";"):
            kv = attr.strip().split()
            if len(kv) >= 2:
                key = kv[0]
                val = kv[1].replace('"', "")
                if key == "gene_id":
                    gene_id = val
                elif key == "gene_biotype":
                    biotype = val
        return gene_id, biotype

    parsed = gtf_df["attribute"].apply(parse_attr)
    gtf_df["gene_id"] = parsed.apply(lambda x: x[0])
    gtf_df["biotype"] = parsed.apply(lambda x: x[1] if x[1] else "unknown")

    return gtf_df[["gene_id", "biotype"]].drop_duplicates()


def load_gene_intervals(gene_bed_path: str, gtf_path: str) -> Dict[Tuple[str, str], ChromosomeData]:
    gene_df = pd.read_csv(
        gene_bed_path,
        sep="\t",
        header=None,
        names=["chr","start","end","gene_id","source","strand"],
        dtype={"chr": str, "start": np.int64, "end": np.int64, "gene_id": str, "strand": str, "source": str},
    )

    biotype_df = parse_gtf_gene_biotypes(gtf_path)
    gene_df = gene_df.merge(biotype_df, on="gene_id", how="left")
    gene_df["biotype"] = gene_df["biotype"].fillna("unknown")

    gene_df["length"] = gene_df["end"] - gene_df["start"]
    gene_df = gene_df.sort_values(["chr","strand","start","length"], ascending=[True,True,True,False])

    out = {}
    for (chrom, strand), g in gene_df.groupby(["chr","strand"]):
        starts = g["start"].values
        ends = g["end"].values
        gene_ids = g["gene_id"].values
        biotypes = g["biotype"].values

        # remove overlaps (keep longest-first)
        keep = np.ones(len(starts), dtype=bool)
        for i in range(1, len(starts)):
            if starts[i] < ends[i-1] and keep[i-1]:
                keep[i] = False

        out[(chrom, strand)] = ChromosomeData(starts[keep], ends[keep], gene_ids[keep], biotypes[keep])

    return out


def find_gene(pos: int, chr_data: ChromosomeData):
    idx = np.searchsorted(chr_data.starts, pos, side="right") - 1
    if idx >= 0 and chr_data.starts[idx] <= pos < chr_data.ends[idx]:
        return chr_data.gene_ids[idx], chr_data.biotypes[idx]
    return None, "intergenic"


def process_chunk(args):
    chunk, chr_data_dict, chrom, strand = args
    chr_data = chr_data_dict[(chrom, strand)]

    positions = chunk["start"].values
    gene_ids = np.empty(len(positions), dtype=object)
    biotypes = np.empty(len(positions), dtype=object)

    for i, p in enumerate(positions):
        gene_ids[i], biotypes[i] = find_gene(int(p), chr_data)

    out = chunk.copy()
    out["gene_id"] = gene_ids
    out["biotype"] = biotypes
    return out


def main():
    base = pd.read_parquet(os.path.join(PLOTS_DIR, IN_PARQUET))
    base = base.sort_values(["chr","strand","start"]).reset_index(drop=True)

    chr_data = load_gene_intervals(GENE_BED, GTF)

    # build chunk list
    tasks = []
    for (chrom, strand), g in base.groupby(["chr","strand"]):
        if (chrom, strand) not in chr_data:
            continue
        for i in range(0, len(g), CHUNK_SIZE):
            tasks.append((g.iloc[i:i+CHUNK_SIZE].copy(), chr_data, chrom, strand))

    results = []
    with ProcessPoolExecutor(max_workers=N_CORES) as ex:
        futures = [ex.submit(process_chunk, t) for t in tasks]
        for fut in tqdm.tqdm(futures, total=len(futures), desc="Assign genes"):
            results.append(fut.result())

    out_df = pd.concat(results, ignore_index=True)
    out_df = out_df.sort_values(["chr","strand","start"]).reset_index(drop=True)

    out_pq = os.path.join(PLOTS_DIR, OUT_PARQUET)
    out_tsv = os.path.join(PLOTS_DIR, OUT_TSV)
    out_df.to_parquet(out_pq, index=False, compression="snappy")
    out_df.to_csv(out_tsv, sep="\t", index=False)

    assigned = out_df["gene_id"].notna().sum()
    print(f"Assigned {assigned:,} of {len(out_df):,} sites ({assigned/len(out_df)*100:.1f}%)")
    print(f"Saved: {out_pq}")
    print(f"Saved: {out_tsv}")


if __name__ == "__main__":
    main()
