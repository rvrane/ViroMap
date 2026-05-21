#!/bin/bash

############################################
# filter_probes_by_human_hits.sh
#
# Remove probe sequences that have any BLAST hits
# against the human genome (above given thresholds).
############################################

set -euo pipefail

# -------- PARAMETERS --------

# Directory containing probe FASTA files
PROBE_DIR="/path/to/probe_results"
# Directory containing BLAST outputs vs. human genome
BLAST_DIR="${PROBE_DIR}/blast_results"
# Output directory for filtered probes
FILTERED_DIR="${PROBE_DIR}/filtered_probes"
# Minimum percent identity to consider a “hit”
MIN_IDENTITY=80
# Minimum alignment length (nt) to consider a “hit”
MIN_ALIGN_LEN=30

THREADS=4

############################################
# PREPARE
############################################

mkdir -p "$FILTERED_DIR"

# For each probe FASTA + corresponding BLAST result
for PROBE_FASTA in "$PROBE_DIR"/*.fasta; do
  BASENAME=$(basename "$PROBE_FASTA" .fasta)
  BLAST_OUT="${BLAST_DIR}/${BASENAME}_vs_human.txt"
  FILTERED_FASTA="${FILTERED_DIR}/${BASENAME}_filtered.fasta"

  # If no BLAST output, just copy original
  if [ ! -s "$BLAST_OUT" ]; then
    echo "No human hits for $BASENAME – copying all probes"
    cp "$PROBE_FASTA" "$FILTERED_FASTA"
    continue
  fi

  # Identify probe IDs to remove:
  # BLAST outfmt 6 columns: qseqid sseqid pident length ...
  awk -v pid="$MIN_IDENTITY" -v plen="$MIN_ALIGN_LEN" '
    $3 >= pid && $4 >= plen { print $1 }
  ' "$BLAST_OUT" | sort -u > "${FILTERED_DIR}/${BASENAME}_to_remove.txt"

  # Filter FASTA: exclude any header matching to_remove
  awk '
    BEGIN { while (getline id < "'"${FILTERED_DIR}/${BASENAME}_to_remove.txt"'" ) rem[id]=1 }
    /^>/ {
      header = substr($0,2)
      split(header, a, /[ \t]/)
      seqid = a[1]
      keep = !(seqid in rem)
    }
    {
      if (keep) print
    }
  ' "$PROBE_FASTA" > "$FILTERED_FASTA"

  echo "Filtered $BASENAME: removed $(wc -l < "${FILTERED_DIR}/${BASENAME}_to_remove.txt") probes, kept $(grep -c '^>' "$FILTERED_FASTA") probes."
done

echo "Filtering complete. Filtered FASTAs in: $FILTERED_DIR"
