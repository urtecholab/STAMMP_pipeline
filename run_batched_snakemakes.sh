#!/bin/bash

SNAKEFILE="capscan_batch_snakemake_v2"
CORES=32
BATCH_DIR="/home/maholton/archive/raw_data/2026/capscan_STAMMP_pipeline/inputs/agora_inputs"
AGORA_DIR="outputs/agora_outputs/AGORA_mapped"
COUNTS_TABLE="${AGORA_DIR}/agora_mapped_counts.csv"

# Find all batch directories in order
batches=($(ls -d ${BATCH_DIR}/batch_*/  | sort))

echo "Found ${#batches[@]} batches: ${batches[@]}"

# Run each batch sequentially
for batch in "${batches[@]}"; do
    batch="${batch%/}"  # strip trailing slash
    echo "=========================================="
    echo "Running batch: $batch"
    echo "=========================================="

    # Temporarily override RAW_DIR for this batch by passing it as a config value
    snakemake -s "$SNAKEFILE" \
        --cores $CORES \
        --rerun-incomplete \
        --config raw_dir="$batch" \
        --until minimap_agora_count

    if [ $? -ne 0 ]; then
        echo "ERROR: Batch $batch failed. Stopping."
        exit 1
    fi

    echo "Batch $batch complete."
done

echo "=========================================="
echo "All batches done. Running final merge..."
echo "=========================================="

snakemake -s "$SNAKEFILE" --cores $CORES --config merge_only=True