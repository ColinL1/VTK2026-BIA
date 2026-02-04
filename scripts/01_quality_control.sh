#!/bin/bash
# Week 1 - Day 1-2: Quality Control and Preprocessing
# Raw data assessment, adapter trimming, and read filtering

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load configuration
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/." && pwd)"

# Parse command line arguments
SAMPLE_PATH=${1:-""}
THREADS=${2:-8}

if [ -z "$SAMPLE_PATH" ]; then
    echo "Usage: $0 <sample_path> [threads]"
    echo "Example: $0 F003_M_enclense.fastq.gz 8"
    exit 1
fi

# Get paths from sample info
RAW_DATA_DIR="${BASE_DIR}/."
# RAW_DATA_DIR="${BASE_DIR}/data/mock_data/fastq_out"
QC_DIR="${BASE_DIR}/results/qc"
TRIMMED_DIR="${QC_DIR}/trimmed"
FILTERED_DIR="${QC_DIR}/filtered"

# Create output directories
mkdir -p "${QC_DIR}/nanoplot_raw"
mkdir -p "${QC_DIR}/nanoplot_filtered"
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${FILTERED_DIR}"

echo "========================================="
echo "Quality Control Pipeline for ${SAMPLE_PATH}"
echo "========================================="

# Find input fastq file
INPUT_FASTQ=$(find "${RAW_DATA_DIR}" -name "${SAMPLE_PATH}" -o -name "${SAMPLE_PATH}*.fastq.gz" | head -n 1)

if [ -z "$INPUT_FASTQ" ]; then
    echo "Error: No input file found for ${SAMPLE_PATH} in ${RAW_DATA_DIR}"
    exit 1
fi

echo "Input file: ${INPUT_FASTQ}"

# Extract sample name without extension for directory naming
SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .fastq.gz)
SAMPLE_NAME=$(basename "${SAMPLE_NAME}" .fastq)

# Step 1: Raw data assessment with NanoPlot
echo ""
echo "Step 1: Raw data assessment with NanoPlot..."
NanoPlot -t ${THREADS} \
    --fastq "${INPUT_FASTQ}" \
    --loglength \
    --plots dot \
    -o "${QC_DIR}/nanoplot_raw/${SAMPLE_NAME}" \
    --N50 \
    --title "${SAMPLE_NAME} - Raw Reads"

# Generate summary statistics with NanoStat
NanoStat --fastq "${INPUT_FASTQ}" \
    -t ${THREADS} \
    --name ""${QC_DIR}/${SAMPLE_PATH}_raw_stats.txt"

echo "✓ Raw data assessment complete"

# Step 2: Adapter trimming with Porechop_ABI
echo ""
echo "Step 2: Adapter trimming with Porechop_ABI..."
TRIMMED_FASTQ="${TRIMMED_DIR}/${SAMPLE_PATH}_trimmed.fastq.gz"

porechop_abi \
    -i "${INPUT_FASTQ}" \
    -o "${TRIMMED_FASTQ}" \
    --threads ${THREADS} \
    -abi \
    --discard_middle

echo "✓ Adapter trimming complete"

# # Step 3: Calculate genome size for coverage estimation
# # Microbacterium: 3-5 MB, Alteromonas: 4-5 MB
# # We'll use 4 MB as average for 100x target
GENOME_SIZE=4500000
TARGET_COV=150
TARGET_BASES=$((GENOME_SIZE * TARGET_COV))

echo ""
echo "Step 3: Read filtering with Filtlong..."
echo "Target coverage: ${TARGET_COV}x"
echo "Estimated genome size: ${GENOME_SIZE} bp"
echo "Target bases: ${TARGET_BASES} bp"

FILTERED_FASTQ="${FILTERED_DIR}/${SAMPLE_PATH}_filtered.fastq.gz"

filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    --target_bases ${TARGET_BASES} \
    --length_weight 10 \
    --mean_q_weight 10 \
    "${TRIMMED_FASTQ}" | gzip > "${FILTERED_FASTQ}"

echo "✓ Read filtering complete"

# Step 4: Post-filtering assessment
echo ""
echo "Step 4: Post-filtering assessment with NanoPlot..."
NanoPlot -t ${THREADS} \
    --fastq "${FILTERED_FASTQ}" \
    --loglength \
    --plots dot \
    -o "${QC_DIR}/nanoplot_filtered/${SAMPLE_PATH}" \
    --N50 \
    --title "${SAMPLE_PATH} - Filtered Reads"

NanoStat --fastq "${FILTERED_FASTQ}" \
    --name "${QC_DIR}/${SAMPLE_PATH}_filtered_stats.txt" \
    -t ${THREADS} 

echo "✓ Post-filtering assessment complete"

# Generate summary report
echo ""
echo "========================================="
echo "QC Summary for ${SAMPLE_PATH}"
echo "========================================="
echo ""
echo "Raw reads statistics:"
grep -E "Number of reads|Total bases|Mean read length|Mean read quality|Read length N50" \
    "${QC_DIR}/${SAMPLE_PATH}_raw_stats.txt" || true
echo ""
echo "Filtered reads statistics:"
grep -E "Number of reads|Total bases|Mean read length|Mean read quality|Read length N50" \
    "${QC_DIR}/${SAMPLE_PATH}_filtered_stats.txt" || true
echo ""
echo "Output files:"
echo "  - Trimmed reads: ${TRIMMED_FASTQ}"
echo "  - Filtered reads: ${FILTERED_FASTQ}"
echo "  - QC reports: ${QC_DIR}/nanoplot_filtered/${SAMPLE_PATH}"
echo ""
echo "✓ Quality control pipeline complete for ${SAMPLE_PATH}"
echo "========================================="
