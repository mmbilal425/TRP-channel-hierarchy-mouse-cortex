#!/usr/bin/env python3
"""
01_attach_motifs.py

Adds 5-mer motif centered at the modified base (start position),
and flags DRACH / RAC motifs.

Input: merged_m6a_pileups.parquet
Output: merged_m6a_pileups_with_motifs.parquet
"""

import re
import time
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
from tqdm import tqdm

# ----------------------------
# CONFIG (your defaults)
# ----------------------------
DATA_DIR = "/g/data/lf10/mb1232/nanopore_data/2024_neurotranscriptomics/data"
IN_PARQUET = "merged_m6a_pileups.parquet"
OUT_PARQUET = "merged_m6a_pileups_with_motifs.parquet"

REF_FASTA = "/g/data/lf10/mb1232/reference_genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa"

NWORKERS = 32
CHUNKSIZE = 4_590_000

# patterns (same logic as you used)
DRACH_PATTERN = re.compile(r'^[AGT][AG]A[C][ATC]$')
RAC_PATTERN   = re.compile(r'^[ACGT][AG]A[C][ACGT]$')

genome = None


def init_worker():
    global genome
    genome = Fasta(REF_FASTA, as_raw=True, build_index=False)


def extract_one(chr_, start, strand):
    """
    Extract 5-mer centered on start:
      motif_start = start-2, motif_end = start+3
    """
    try:
        chr_ = str(chr_)
        start = int(start)
        if start < 2:
            return "NNNNN", 0, 0

        motif_start = start - 2
        motif_end = start + 3

        seq = str(genome[chr_][motif_start:motif_end]).upper()
        if len(seq) != 5:
            return "NNNNN", 0, 0

        if strand == "-":
            seq = str(Seq(seq).reverse_complement())

        if len(seq) != 5:
            return "NNNNN", 0, 0

        is_drach = int(bool(DRACH_PATTERN.match(seq)))
        is_rac = int(bool(RAC_PATTERN.match(seq)))

        return seq, is_drach, is_rac

    except Exception:
        return "NNNNN", 0, 0


def process_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    chrs = chunk["chr"].astype(str).values
    starts = chunk["start"].values
    strands = chunk["strand"].astype(str).values

    motifs = []
    drach = []
    rac = []

    for c, s, st in zip(chrs, starts, strands):
        m, d, r = extract_one(c, s, st)
        motifs.append(m)
        drach.append(d)
        rac.append(r)

    out = chunk.copy()
    out["motif"] = motifs
    out["is_drach"] = np.array(drach, dtype=np.int8)
    out["is_rac"] = np.array(rac, dtype=np.int8)
    return out


def main():
    t0 = time.time()
    data_dir = Path(DATA_DIR)

    df = pd.read_parquet(data_dir / IN_PARQUET)
    df = df.sort_values(["chr", "start", "strand"]).reset_index(drop=True)

    # chunking
    chunks = [df[i:i + CHUNKSIZE] for i in range(0, len(df), CHUNKSIZE)]

    out_chunks = []
    with Pool(processes=min(NWORKERS, cpu_count()), initializer=init_worker) as pool:
        for out_chunk in tqdm(pool.imap(process_chunk, chunks), total=len(chunks), desc="Motifs"):
            out_chunks.append(out_chunk)

    out_df = pd.concat(out_chunks, ignore_index=True)
    out_df = out_df.sort_values(["chr", "start", "strand"]).reset_index(drop=True)

    out_path = data_dir / OUT_PARQUET
    out_df.to_parquet(out_path, index=False, compression="snappy")

    print(f"Saved: {out_path}")
    print(f"Rows={len(out_df):,} time={time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
