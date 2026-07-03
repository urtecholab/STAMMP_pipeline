#!/bin/bash

BASE="/home/gurtecho/raw_data/2026/20260417_SDDRC_IBD_Biobank/"
OUTPUT="/home/maholton/archive/raw_data/2026/Marvic_IBD_biobank_STAMMP_pipeline/inputs/agora_inputs"

# Get sorted unique sample names
samples=($(ls ${BASE}/*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' | sort))

batch_size=5
batch_num=1
count=0

# format for BASE is SRR15540846_1.fastq.gz

for sample in "${samples[@]}"; do
    # Start new batch directory if needed
    if [ $count -eq 0 ]; then
        dir="${BASE}/batch_$(printf '%02d' $batch_num)"
        mkdir -p "$dir"
    fi

    mv "${BASE}/${sample}_1.fastq.gz" "$dir/"
    mv "${BASE}/${sample}_2.fastq.gz" "$dir/"

    count=$((count + 1))

    if [ $count -eq $batch_size ]; then
        count=0
        batch_num=$((batch_num + 1))
    fi
done

echo "Done. Created $batch_num batches."