#!/bin/bash
# Week 1 - Day 3-4: Genome Assembly
# De novo assembly with Flye, polishing with Medaka, assembly QC

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load configuration
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/." && pwd)"

# Parse command line arguments
SAMPLE_PATH=${1:-""}
THREADS=${2:-8}
GENOME_SIZE=${3:-"4m"}  # 4m for 4 Mb, adjust based on expected size

if [ -z "$SAMPLE_PATH" ]; then
    echo "Usage: $0 <sample_path> [threads] [genome_size]"
    echo "Example: $0 F003_M_enclense_filtered.fastq.gz 8 4m"
    echo "  genome_size: 3m for Microbacterium, 4.5m for Alteromonas"
    exit 1
fi

# Get paths
FILTERED_DIR="${BASE_DIR}/results/qc/filtered"
ASSEMBLY_DIR="${BASE_DIR}/results/assembly"

echo "========================================="
echo "Assembly Pipeline for ${SAMPLE_PATH}"
echo "========================================="

# Find filtered reads
FILTERED_FASTQ=$(find "${FILTERED_DIR}" -name "${SAMPLE_PATH}" -o -name "${SAMPLE_PATH}*.fastq.gz" | head -n 1)

if [ -z "$FILTERED_FASTQ" ]; then
    echo "Error: No filtered reads found for ${SAMPLE_PATH} in ${FILTERED_DIR}"
    echo "Please run 01_quality_control.sh first"
    exit 1
fi

echo "Input file: ${FILTERED_FASTQ}"

# Extract sample name without extension for directory naming
SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .fastq.gz)
SAMPLE_NAME=$(basename "${SAMPLE_NAME}" .fastq)
SAMPLE_NAME=$(basename "${SAMPLE_NAME}" _filtered)

# Create output directories
mkdir -p "${ASSEMBLY_DIR}/${SAMPLE_NAME}/flye"
mkdir -p "${ASSEMBLY_DIR}/${SAMPLE_NAME}/medaka"
mkdir -p "${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc"

echo "Sample name: ${SAMPLE_NAME}"
echo "Expected genome size: ${GENOME_SIZE}"

# Step 1: De novo assembly with Flye
echo ""
echo "Step 1: De novo assembly with Flye..."
echo "This may take 1-3 hours depending on coverage and genome size"

flye \
    --nano-hq "${FILTERED_FASTQ}" \
    --genome-size ${GENOME_SIZE} \
    --threads ${THREADS} \
    --out-dir "${ASSEMBLY_DIR}/${SAMPLE_NAME}/flye" \
    --iterations 4 \
    --meta

echo "✓ Flye assembly complete"

# Copy assembly for polishing
cp "${ASSEMBLY_DIR}/${SAMPLE_NAME}/flye/assembly.fasta" \
   "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_flye.fasta"

# Step 2: Polish with Medaka
echo ""
echo "Step 2: Polishing with Medaka..."
echo "This may take 30-60 minutes"

# Determine Medaka model based on basecaller (adjust as needed)
# For newer basecallers, use r1041_e82_400bps_sup_v5.2.0
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.2.0"

medaka_consensus \
    -i "${FILTERED_FASTQ}" \
    -d "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_flye.fasta" \
    -o "${ASSEMBLY_DIR}/${SAMPLE_NAME}/medaka" \
    -t ${THREADS} \
    -m ${MEDAKA_MODEL}

# Copy polished assembly
cp "${ASSEMBLY_DIR}/${SAMPLE_NAME}/medaka/consensus.fasta" \
   "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_polished.fasta"

echo "✓ Medaka polishing complete"

# Step 3: Assembly QC with QUAST
echo ""
echo "Step 3: Assembly quality assessment with QUAST..."

quast.py \
    "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_flye.fasta" \
    "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_polished.fasta" \
    -o "${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc/quast" \
    --threads ${THREADS} \
    --min-contig 500 \
    --labels "Flye,Flye+Medaka"

echo "✓ QUAST complete"

# Step 4: BUSCO completeness assessment
echo ""
echo "Step 4: BUSCO completeness assessment..."
echo "Using Bacteria_odb10 lineage dataset"

# Download BUSCO datasets if not present
if [ ! -d "${HOME}/.busco_downloads/lineages/bacteria_odb10" ]; then
    echo "Downloading BUSCO bacteria dataset..."
    busco --download bacteria_odb10
fi

busco \
    -i "${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_polished.fasta" \
    -o "${SAMPLE_NAME}_busco" \
    --out_path "${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc" \
    -m genome \
    -l bacteria_odb10 \
    --cpu ${THREADS} \
    --offline

echo "✓ BUSCO complete"

# Generate assembly summary
echo ""
echo "========================================="
echo "Assembly Summary for ${SAMPLE_NAME}"
echo "========================================="
echo ""
echo "QUAST Results (Polished Assembly):"
grep -A 20 "Flye+Medaka" "${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc/quast/report.txt" | head -n 15 || true
echo ""
echo "BUSCO Results:"
cat "${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc/${SAMPLE_NAME}_busco/short_summary*.txt" | \
    grep -E "C:|S:|D:|F:|M:" || true
echo ""
echo "Assembly files:"
echo "  - Flye assembly: ${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_flye.fasta"
echo "  - Polished assembly: ${ASSEMBLY_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_polished.fasta"
echo "  - Assembly graph: ${ASSEMBLY_DIR}/${SAMPLE_NAME}/flye/assembly_graph.gfa"
echo "  - QUAST report: ${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc/quast/report.html"
echo "  - BUSCO report: ${ASSEMBLY_DIR}/${SAMPLE_NAME}/qc/${SAMPLE_NAME}_busco/"
echo ""
echo "✓ Assembly pipeline complete for ${SAMPLE_NAME}"
echo "========================================="
