#!/bin/bash
# Week 1 - Day 5: Taxon-Independent Functional Annotation
# Reference-free annotation with Bakta

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load configuration
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/." && pwd)"

# Parse command line arguments
SAMPLE_PATH=${1:-""}
THREADS=${2:-8}

if [ -z "$SAMPLE_PATH" ]; then
    echo "Usage: $0 <sample_path> [threads]"
    echo "Example: $0 F003_M_enclense_polished.fasta 8"
    exit 1
fi

# Get paths
ASSEMBLY_DIR="${BASE_DIR}/results/assembly"
ANNOTATION_DIR="${BASE_DIR}/results/annotation"

echo "========================================="
echo "Annotation Pipeline for ${SAMPLE_PATH}"
echo "========================================="

# Find polished assembly
POLISHED_ASSEMBLY=$(find "${ASSEMBLY_DIR}" -name "${SAMPLE_PATH}" -o -name "${SAMPLE_PATH}*.fasta" | head -n 1)

if [ -z "$POLISHED_ASSEMBLY" ] || [ ! -f "$POLISHED_ASSEMBLY" ]; then
    echo "Error: Polished assembly not found for ${SAMPLE_PATH} in ${ASSEMBLY_DIR}"
    echo "Please run 02_assembly.sh first"
    exit 1
fi

echo "Input file: ${POLISHED_ASSEMBLY}"

# Extract sample name without extension for directory naming
SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .fasta)
SAMPLE_NAME=$(basename "${SAMPLE_NAME}" .fa)
SAMPLE_NAME=$(basename "${SAMPLE_NAME}" _polished)

# Create output directory
mkdir -p "${ANNOTATION_DIR}/${SAMPLE_NAME}"

echo "Sample name: ${SAMPLE_NAME}"

# Check if Bakta database is set up
if [ -z "${BAKTA_DB:-}" ]; then
    echo ""
    echo "WARNING: BAKTA_DB environment variable not set"
    echo "Please download the Bakta database:"
    echo "  bakta_db download --output <db_path> --type full"
    echo "Then set: export BAKTA_DB=<db_path>"
    echo ""
    read -p "Enter Bakta database path (or press Enter to skip): " BAKTA_DB_INPUT
    BAKTA_DB=${BAKTA_DB_INPUT:-""}
fi

if [ -z "$BAKTA_DB" ]; then
    echo "Error: Bakta database path required"
    exit 1
fi

# Step 1: Annotation with Bakta
echo ""
echo "Step 1: Annotation with Bakta..."
echo "This may take 30-60 minutes depending on genome size"

bakta \
    --db "${BAKTA_DB}" \
    --output "${ANNOTATION_DIR}/${SAMPLE_NAME}" \
    --prefix "${SAMPLE_NAME}" \
    --threads ${THREADS} \
    --verbose \
    --keep-contig-headers \
    --compliant \
    "${POLISHED_ASSEMBLY}"

echo "✓ Bakta annotation complete"

# Step 2: Generate annotation statistics
echo ""
echo "Step 2: Generating annotation statistics..."

STATS_FILE="${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_annotation_stats.txt"

cat > "${STATS_FILE}" <<EOF
=========================================
Annotation Statistics for ${SAMPLE_NAME}
=========================================

EOF

# Extract statistics from Bakta JSON output
if command -v python3 &> /dev/null; then
    python3 - << PYTHON_EOF >> "${STATS_FILE}"
import json
import sys

try:
    with open("${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.json", 'r') as f:
        data = json.load(f)
    
    stats = data.get('stats', {})
    
    print(f"Genome size: {stats.get('size', 'N/A'):,} bp")
    print(f"GC content: {stats.get('gc', 'N/A'):.2f}%")
    print(f"Number of contigs: {stats.get('no_sequences', 'N/A')}")
    print(f"N50: {stats.get('n50', 'N/A'):,} bp")
    print()
    print("Gene Counts:")
    print(f"  Total genes: {stats.get('no_genes', 'N/A')}")
    print(f"  CDS: {stats.get('no_cds', 'N/A')}")
    print(f"  tRNA: {stats.get('no_t_rna', 'N/A')}")
    print(f"  rRNA: {stats.get('no_r_rna', 'N/A')}")
    print(f"  ncRNA: {stats.get('no_nc_rna', 'N/A')}")
    print(f"  CRISPR arrays: {stats.get('no_crispr', 'N/A')}")
    
except Exception as e:
    print(f"Could not parse JSON: {e}")
PYTHON_EOF
fi

# Count hypothetical proteins
HYPOTHETICAL_COUNT=$(grep -c "hypothetical protein" "${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.tsv" || echo "0")
TOTAL_CDS=$(grep -c "CDS" "${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.gff3" || echo "1")
HYPOTHETICAL_PERCENT=$(echo "scale=2; ${HYPOTHETICAL_COUNT} * 100 / ${TOTAL_CDS}" | bc || echo "N/A")

cat >> "${STATS_FILE}" <<EOF

Functional Annotation:
  Hypothetical proteins: ${HYPOTHETICAL_COUNT} (${HYPOTHETICAL_PERCENT}%)
  Functionally annotated: $((TOTAL_CDS - HYPOTHETICAL_COUNT))

EOF

echo "✓ Statistics generation complete"

# Display summary
echo ""
echo "========================================="
echo "Annotation Summary for ${SAMPLE_NAME}"
echo "========================================="
cat "${STATS_FILE}"
echo ""
echo "Output files:"
echo "  - GenBank: ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.gbff"
echo "  - GFF3: ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.gff3"
echo "  - FASTA (proteins): ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.faa"
echo "  - FASTA (nucleotides): ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.fna"
echo "  - TSV summary: ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.tsv"
echo "  - JSON: ${ANNOTATION_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.json"
echo "  - Statistics: ${STATS_FILE}"
echo ""
echo "✓ Annotation pipeline complete for ${SAMPLE_NAME}"
echo "========================================="
