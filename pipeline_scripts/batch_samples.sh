#!/bin/bash
# Usage: batch_samples.sh <dir> <r1_pattern> [batch_size]
#
# Moves paired-end fastq files in <dir> into numbered batch subdirectories.
# Useful for splitting a large sample set into manageable chunks for a pipeline.
#
# Arguments:
#   <dir>         Directory containing the fastq files
#   <r1_pattern>  Glob pattern matching R1 files (e.g. "*_1.fastq.gz")
#                 R2 files are inferred by replacing _1 with _2 in each filename
#   [batch_size]  Number of samples per batch (default: 5)
#
# Output: batch_01/, batch_02/, ... subdirectories created inside <dir>,
#         each containing up to <batch_size> sample pairs (R1 + R2).

# --- Argument validation ---

# Require at least 2 arguments (dir and pattern); batch size is optional
if [ $# -lt 2 ]; then
    echo "Usage: $(basename "$0") <dir> <r1_pattern> [batch_size]"
    echo "  Example: $(basename "$0") /data/fastqs '*_1.fastq.gz' 5"
    exit 1
fi

DIR="$1"
R1_PATTERN="$2"
# Use $3 if provided, otherwise default to 5
BATCH_SIZE="${3:-5}"

# Confirm the target directory exists before proceeding
if [ ! -d "$DIR" ]; then
    echo "Error: directory not found: $DIR"
    exit 1
fi

# --- Collect R1 files ---

# Expand the glob inside DIR, suppress errors if nothing matches (2>/dev/null),
# then sort so batches are assigned in a consistent, reproducible order
r1_files=($(ls "${DIR}"/${R1_PATTERN} 2>/dev/null | sort))

# Exit early if the pattern matched nothing — likely a typo or wrong directory
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "Error: no files found matching '${R1_PATTERN}' in ${DIR}"
    exit 1
fi

echo "Found ${#r1_files[@]} R1 files. Batching into groups of ${BATCH_SIZE}."

# --- Batch counter state ---

batch_num=1   # which batch we are currently filling (used to name the directory)
count=0       # how many samples have been added to the current batch

# --- Main loop: assign each R1 (and its R2 pair) to a batch directory ---

for r1 in "${r1_files[@]}"; do

    # Derive the R2 filepath from the R1 filepath.
    # sed replaces the first occurrence of _1 immediately before .fastq with _2,
    # which handles naming conventions like sample_1.fastq.gz -> sample_2.fastq.gz.
    # The \(\.fastq\) capture group ensures we only match _1 at the read-identifier
    # position and not elsewhere in the sample name.
    r2=$(echo "$r1" | sed 's/_1\(\.fastq\)/_2\1/')

    # Skip this sample if the expected R2 file does not exist
    if [ ! -f "$r2" ]; then
        echo "Warning: R2 file not found for $(basename "$r1"), skipping."
        continue
    fi

    # When count is 0 we are starting a fresh batch: create its directory.
    # printf '%02d' zero-pads the number so directories sort correctly (01, 02, ..., 10).
    if [ $count -eq 0 ]; then
        batch_dir="${DIR}/batch_$(printf '%02d' $batch_num)"
        mkdir -p "$batch_dir"
    fi

    # Move both reads for this sample into the current batch directory
    mv "$r1" "$batch_dir/"
    mv "$r2" "$batch_dir/"

    count=$((count + 1))

    # Once the batch is full, reset the counter and advance to the next batch number
    if [ $count -eq $BATCH_SIZE ]; then
        count=0
        batch_num=$((batch_num + 1))
    fi
done

# --- Final report ---

# If count > 0, the last batch was partial and batch_num was never incremented past it,
# so total = batch_num. If count == 0, the last increment already moved batch_num one
# ahead of the last real batch, so subtract 1.
if [ $count -gt 0 ]; then
    total_batches=$batch_num
else
    total_batches=$((batch_num - 1))
fi

echo "Done. Created ${total_batches} batch(es) in ${DIR}."