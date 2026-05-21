#!/bin/bash
# cluster_sequences.sh
# Script to cluster HCV sequences using CD-HIT at different similarity thresholds

# Set parameters
INPUT_FILE="/media/user/36d81ee3-f383-4bc4-b619-20abd785be33/Hepatitis_probe/hep_C/HEP_C_sequences.fasta"
OUTPUT_DIR="data/clustered"
THREADS=8

# Create output directory
mkdir -p $OUTPUT_DIR

# Function to run CD-HIT with appropriate word size based on similarity threshold
cluster_sequences() {
    local similarity=$1
    local word_size=$2
    local output_suffix=$3

    echo "Clustering at ${similarity}% similarity..."
    cd-hit-est -i $INPUT_FILE \
        -o $OUTPUT_DIR/hcv_clustered_${output_suffix}.fasta \
        -c $similarity \
        -n $word_size \
        -l 1000 \
        -d 0 \
        -M 16000 \
        -T $THREADS

    echo "  Completed ${similarity}% clustering"
    echo ""
}

# Cluster at different similarity thresholds
# Note: Word size (-n) should be adjusted based on similarity threshold
# For >= 90% similarity: use -n 10
# For 85-89% similarity: use -n 8  
# For 80-84% similarity: use -n 7

#cluster_sequences 0.98 10 "98"
cluster_sequences 0.96 10 "96" 
#cluster_sequences 0.94 10 "94"
#cluster_sequences 0.92 10 "92"
#cluster_sequences 0.90 10 "90"
#cluster_sequences 0.88 8 "88"

# Generate clustering statistics
echo "========================================="
echo "CLUSTERING STATISTICS SUMMARY"
echo "========================================="
echo "Original sequences: $(grep -c '^>' $INPUT_FILE)"
echo ""
echo "Clustered sequences by similarity threshold:"
#echo "98% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_98.fasta)"
echo "96% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_96.fasta)"
#echo "94% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_94.fasta)"
#echo "92% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_92.fasta)"
#echo "90% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_90.fasta)"
#echo "88% clusters: $(grep -c '^>' $OUTPUT_DIR/hcv_clustered_88.fasta)"
echo ""
echo "Note: Lower similarity thresholds result in fewer, more diverse clusters"
echo "Higher similarity thresholds result in more, more specific clusters"
