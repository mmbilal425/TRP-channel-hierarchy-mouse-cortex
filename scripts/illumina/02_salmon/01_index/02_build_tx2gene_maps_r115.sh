#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$OUTDIR" && cd "$OUTDIR"

catz() { case "$1" in *.gz) zcat "$1" ;; *) cat "$1" ;; esac; }

# (a) transcript_id -> Ensembl gene_id
catz "$TXOME" | awk '
  /^>/{
    tid=""; gid="";
    split($1,hdr,">"); tid=hdr[2];
    if (match($0,/gene:([A-Za-z0-9\.]+)/,m)) gid=m[1];
    if (tid!="" && gid!="") print tid "\t" gid;
  }' | sort -u > tx2gene_ensid.r115.from_fasta.tsv

# (b) transcript_id -> gene_symbol
catz "$TXOME" | awk '
  /^>/{
    tid=""; sym="";
    split($1,hdr,">"); tid=hdr[2];
    if (match($0,/gene_symbol:([^ ]+)/,m)) sym=m[1];
    if (tid!="" && sym!="") print tid "\t" sym;
  }' | sort -u > tx2gene_symbol.r115.from_fasta.tsv

echo "[ensid]"; head tx2gene_ensid.r115.from_fasta.tsv
echo "[symbol]"; head tx2gene_symbol.r115.from_fasta.tsv
echo "[counts]"; wc -l tx2gene_*.r115.from_fasta.tsv
