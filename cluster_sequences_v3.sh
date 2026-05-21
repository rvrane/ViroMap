#!/bin/bash
# cluster_sequences_v3.sh
# Single-threshold (96 %) CD-HIT clustering for viral probe design
set -euo pipefail

# ------- auto-detect reference FASTA in current directory -------
FASTA=$(ls *.fa *.fasta 2>/dev/null | head -n 1)
[ -z "$FASTA" ] && { echo "No FASTA found."; exit 1; }
BASE=$(basename "${FASTA%.*}")

# ------- parameters -------
SIM=0.96        # 96 % similarity
WORD=10         # word size appropriate for ≥90 %
THREADS=8
OUTDIR="data/clustered"
mkdir -p "$OUTDIR"

echo "Clustering $FASTA at 96 % identity (-n $WORD)…"
cd-hit-est -i "$FASTA" -o "$OUTDIR/${BASE}_clustered_96.fasta" \
           -c $SIM -n $WORD -l 1000 -d 0 -M 16000 -T $THREADS
echo "Finished clustering"

# ------- summary -------
echo -e "\n========== CLUSTERING SUMMARY =========="
echo "Input sequences : $(grep -c '^>' "$FASTA")"
echo "96 % clusters    : $(grep -c '^>' "$OUTDIR/${BASE}_clustered_96.fasta")"
echo "Output file      : $OUTDIR/${BASE}_clustered_96.fasta"
echo "========================================"
