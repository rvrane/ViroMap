#!/bin/bash

############################################
# Configuration
############################################

# NCBI reference genome fasta (use correct path/file!)
REF_FASTA="/media/user/36d81ee3-f383-4bc4-b619-20abd785be33/Hepatitis_probe/hep_C/HEP_C_sequences.fasta"

# Directory containing probe FASTA files
PROBE_DIR="/media/user/36d81ee3-f383-4bc4-b619-20abd785be33/Hepatitis_probe/hep_C/data/clustered/probe_results"

# BLAST results will be saved here
BLAST_DIR="${PROBE_DIR}/blast_results"

# Threads for BLAST
THREADS=4

############################################
# Validate input directory
############################################

if [ ! -d "$PROBE_DIR" ]; then
    echo "Error: Probe directory not found: $PROBE_DIR"
    exit 1
fi

# Count FASTA files in probe directory
FASTA_COUNT=$(find "$PROBE_DIR" -maxdepth 1 -name "*.fasta" | wc -l)
if [ "$FASTA_COUNT" -eq 0 ]; then
    echo "Error: No FASTA files found in $PROBE_DIR"
    exit 1
fi
echo "Found $FASTA_COUNT FASTA files in $PROBE_DIR"

############################################
# Make BLAST database (only if not already created)
############################################

if [ ! -f "${REF_FASTA}.nin" ]; then
    echo "Building BLAST database from $REF_FASTA"
    if ! makeblastdb -in "$REF_FASTA" -dbtype nucl -out "${REF_FASTA%.fasta}_db"; then
        echo "Error: Failed to create BLAST database"
        exit 1
    fi
else
    echo "BLAST database already exists for $REF_FASTA"
fi

DB="${REF_FASTA%.fasta}_db"

############################################
# BLAST each probe FASTA file
############################################

mkdir -p "$BLAST_DIR"

for PROBE_FASTA in "$PROBE_DIR"/*.fasta; do
    [ -f "$PROBE_FASTA" ] || continue
    
    # Skip if file is empty
    if [ ! -s "$PROBE_FASTA" ]; then
        echo "Warning: Empty file, skipping: $PROBE_FASTA"
        continue
    fi
    
    # Create output filename based on input filename
    BASENAME=$(basename "$PROBE_FASTA" .fasta)
    OUTFILE="${BLAST_DIR}/${BASENAME}_blast_results.txt"
    
    echo "Running BLAST for: $BASENAME"
    echo "  Input: $PROBE_FASTA"
    echo "  Output: $OUTFILE"
    
    # Run BLAST
    if ! blastn -query "$PROBE_FASTA" -db "$DB" -out "$OUTFILE" -outfmt 6 -num_threads $THREADS; then
        echo "Error: BLAST failed for $PROBE_FASTA"
        continue
    fi
    
    # Check if BLAST produced results
    if [ ! -s "$OUTFILE" ]; then
        echo "Warning: No BLAST hits found for $PROBE_FASTA"
    else
        HIT_COUNT=$(wc -l < "$OUTFILE")
        echo "  Found $HIT_COUNT hits"
    fi
done

# Run BLAST against human genome
HUMAN_DB="/media/user/36d81ee3-f383-4bc4-b619-20abd785be33/Hg38_blast/GRCh38_db"
OUT_HUMAN="${BLAST_DIR}/${BASENAME}_vs_human.txt"
blastn -query "$PROBE_FASTA" -db "$HUMAN_DB" -out "$OUT_HUMAN" -outfmt 6 -num_threads $THREADS


echo "\nBLAST analysis complete. Results saved in: $BLAST_DIR"
