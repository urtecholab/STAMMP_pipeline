#!/bin/bash
# run_batched_snakemakes.sh — STAMMP batched AGORA + MICOM pipeline
#
# Usage:
#   run_batched_snakemakes.sh --input <raw_data_dir> --output <project_dir> [OPTIONS]
#
# Required:
#   --input   <dir>   Directory containing batch_01/, batch_02/, ... subdirectories
#                     (produced by batch_samples.sh)
#   --output  <dir>   Project output root
#                     (e.g. /home/maholton/scratch/projects/2026_HIV_MASLD_STAMMP)
#
# Optional:
#   --cores   <N>     Snakemake threads per batch (default: 32)
#   --skip-micom      Do not run MICOM_pipeline.py after the merge step
#
# STAMMP install directory is auto-detected from this script's own location.
# To share STAMMP, zip the entire STAMMP directory; collaborators unzip and call
# this script directly — no install step needed.
#
# Cluster-specific database paths (KRAKEN_DB, BOWTIE_IDX, AGORA_MINIMAP_DB) are
# set inside capscan_batch_snakemake_v2. Override them per-run with:
#   snakemake ... --config kraken_db=<path> bowtie_idx=<path> agora_db=<path>

set -uo pipefail

# ── Auto-detect STAMMP install directory ──────────────────────────────────────
# realpath resolves $0 against CWD, so this works whether you call the script
# as "bash run_batched_snakemakes.sh" (from pipeline_scripts/) or by full path.
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
STAMMP_DIR="$(realpath "${SCRIPT_DIR}/..")"
SNAKEFILE="${STAMMP_DIR}/pipeline_scripts/capscan_batch_snakemake_v2"

# ── Defaults ──────────────────────────────────────────────────────────────────
CORES=32              # total cores available (passed to snakemake --cores)
THREADS_PER_SAMPLE=8  # threads allocated per rule; controls how many samples run
                      # in parallel: floor(CORES / THREADS_PER_SAMPLE) at once.
                      # e.g. 32 cores / 8 threads = 4 samples simultaneously.
                      # Setting this equal to CORES serializes everything.
BATCH_DIR=""
PROJECT_DIR=""
RUN_MICOM=true

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)      BATCH_DIR="$2";    shift 2 ;;
        --output)     PROJECT_DIR="$2";  shift 2 ;;
        --cores)      CORES="$2";        shift 2 ;;
        --skip-micom) RUN_MICOM=false;   shift   ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $(basename "$0") --input <raw_data_dir> --output <project_dir> [--cores N] [--skip-micom]"
            exit 1
            ;;
    esac
done

if [[ -z "$BATCH_DIR" || -z "$PROJECT_DIR" ]]; then
    echo "ERROR: --input and --output are required."
    echo "Usage: $(basename "$0") --input <raw_data_dir> --output <project_dir> [--cores N] [--skip-micom]"
    exit 1
fi

if [[ ! -d "$BATCH_DIR" ]]; then
    echo "ERROR: input directory not found: $BATCH_DIR"
    exit 1
fi

# ── Convenience paths ─────────────────────────────────────────────────────────
AGORA_OUT="${PROJECT_DIR}/outputs/agora_outputs"
AGORA_MAPPED="${AGORA_OUT}/AGORA_mapped"
SAMPLE_LIST="${AGORA_OUT}/sample_list.txt"
COUNTS_TABLE="${AGORA_MAPPED}/agora_mapped_counts.csv"

# ── Create output directory structure ─────────────────────────────────────────
mkdir -p "${AGORA_OUT}/cleaned"
mkdir -p "${AGORA_OUT}/kraken_filtered"
mkdir -p "${AGORA_OUT}/bowtie_filtered"
mkdir -p "${AGORA_MAPPED}"
mkdir -p "${PROJECT_DIR}/outputs/humann_outputs"
mkdir -p "${PROJECT_DIR}/outputs/micom_outputs"

echo "STAMMP directory : ${STAMMP_DIR}"
echo "Input (batches)  : ${BATCH_DIR}"
echo "Output (project) : ${PROJECT_DIR}"
echo ""

# ── Auto-generate master sample list from all batches ─────────────────────────
find "${BATCH_DIR}/batch_"* -maxdepth 1 -name "*_1.fastq.gz" \
    | xargs -I{} basename {} _1.fastq.gz \
    | sort \
    > "${SAMPLE_LIST}"

N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "Sample list: ${N_SAMPLES} samples -> ${SAMPLE_LIST}"
echo ""

# ── Find and sort batch directories ───────────────────────────────────────────
batches=($(ls -d "${BATCH_DIR}"/batch_*/ 2>/dev/null | sort))

if [[ ${#batches[@]} -eq 0 ]]; then
    echo "ERROR: no batch_*/ directories found in ${BATCH_DIR}"
    exit 1
fi

echo "Found ${#batches[@]} batch(es): ${batches[*]}"
echo ""

# ── Step 1: run each batch through the AGORA rules ────────────────────────────
# (produces per-sample .counts files).
# Batches run sequentially; the script stops on the first failure.
for batch in "${batches[@]}"; do
    batch="${batch%/}"   # strip trailing slash
    echo "=========================================="
    echo "Running batch: ${batch}"
    echo "=========================================="

    if ! snakemake -s "$SNAKEFILE" \
        --cores "$CORES" \
        --rerun-incomplete \
        --config \
            raw_dir="$batch" \
            output_dir="$PROJECT_DIR" \
            stammp_dir="$STAMMP_DIR" \
            sample_list="$SAMPLE_LIST" \
            threads="$THREADS_PER_SAMPLE"; then
        echo "ERROR: batch ${batch} failed. Stopping."
        exit 1
    fi

    echo "Batch ${batch} complete."
    echo ""
done

# ── Step 2: merge all per-sample count files into a single table ──────────────
# All batches are done; collect every .counts file in AGORA_MAPPED and combine.
echo "=========================================="
echo "Merging count files..."
echo "=========================================="

count_files=("${AGORA_MAPPED}"/*.counts)

if [[ ${#count_files[@]} -eq 0 ]]; then
    echo "ERROR: no .counts files found in ${AGORA_MAPPED}"
    exit 1
fi

echo "Merging ${#count_files[@]} count file(s) -> ${COUNTS_TABLE}"

if ! python "${STAMMP_DIR}/agora_filter_scripts/count_files_to_table.py" \
    "${count_files[@]}" \
    > "${COUNTS_TABLE}"; then
    echo "ERROR: merge step failed."
    exit 1
fi

echo "Merge complete: ${COUNTS_TABLE}"
echo ""

# ── Step 3: MICOM (requires the merged counts table from step 2) ──────────────
# micom_references lives inside the STAMMP install — no copies needed.
if [[ "$RUN_MICOM" == true ]]; then
    echo "=========================================="
    echo "Running MICOM pipeline..."
    echo "=========================================="

    if ! python "${STAMMP_DIR}/pipeline_scripts/MICOM_pipeline.py" \
        --agora     "${COUNTS_TABLE}" \
        --output    "${PROJECT_DIR}/outputs" \
        --reference "${STAMMP_DIR}/micom_references"; then
        echo "ERROR: MICOM pipeline failed."
        exit 1
    fi

    echo "MICOM complete."
    echo ""
fi

echo "=========================================="
echo "Pipeline finished."
echo "All outputs in: ${PROJECT_DIR}/outputs"
echo "  agora_outputs/  -> AGORA mapped counts and filtered reads"
echo "  humann_outputs/ -> (HUMAnN outputs, if run separately)"
echo "  micom_outputs/  -> MICOM models and growth results"
echo "=========================================="
