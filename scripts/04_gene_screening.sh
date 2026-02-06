#!/bin/bash
# Gene Screening and Trait Detection
# Secondary metabolite biosynthesis, antimicrobial resistance genes, and coral-associated trait screening

set -e  # Exit on error
set -u  # Exit on undefined variable

# Parse command line arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <assembly_path> <gff_path> <protein_fasta_path> [threads] [output_dir]"
    echo ""
    echo "Example: $0 ./assembly.fasta ./annotation.gff ./proteins.faa 8 ./results"
    echo ""
    echo "Arguments:"
    echo "  assembly_path      Path to polished assembly FASTA file"
    echo "  gff_path           Path to GFF3 annotation file"
    echo "  protein_fasta_path Path to protein FASTA file"
    echo "  threads            Number of CPU threads (default: 8)"
    echo "  output_dir         Output directory (default: ./results)"
    exit 1
fi

SAMPLE_ID="$(basename "${2}" | sed 's/\.[^.]*$//')"
ASSEMBLY="${1}"
GFF="${2}"
PROTEIN_FASTA="${3}"
THREADS="${4:-8}"
OUTPUT_BASE_DIR="${5:-.}/results"

# Convert to absolute paths
ASSEMBLY="$(cd "$(dirname "${ASSEMBLY}")" && pwd)/$(basename "${ASSEMBLY}")"
GFF="$(cd "$(dirname "${GFF}")" && pwd)/$(basename "${GFF}")"
PROTEIN_FASTA="$(cd "$(dirname "${PROTEIN_FASTA}")" && pwd)/$(basename "${PROTEIN_FASTA}")"

# Create output directories
SCREENING_DIR="${OUTPUT_BASE_DIR}/screening/${SAMPLE_ID}"
mkdir -p "${SCREENING_DIR}/antismash"
mkdir -p "${SCREENING_DIR}/amr"

echo "========================================="
echo "Gene Screening Pipeline for ${SAMPLE_ID}"
echo "========================================="
echo "Threads: ${THREADS}"
echo ""

# Verify input files exist
echo "Checking input files..."

if [ ! -f "$ASSEMBLY" ]; then
    echo "Error: Assembly not found: ${ASSEMBLY}"
    exit 1
fi

if [ ! -f "$GFF" ]; then
    echo "Error: GFF file not found: ${GFF}"
    exit 1
fi

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "Error: Protein FASTA not found: ${PROTEIN_FASTA}"
    exit 1
fi

echo "Input assembly: ${ASSEMBLY}"
echo "Input GFF: ${GFF}"
echo "Input proteins: ${PROTEIN_FASTA}"
echo "Output directory: ${SCREENING_DIR}"
echo ""

# Step 1: Secondary metabolite biosynthesis with antiSMASH
echo "Step 1: Detecting biosynthetic gene clusters with antiSMASH..."
echo "This may take long..."
antismash \
    --genefinding-tool none \
    --genefinding-gff3 "${GFF}" \
    --output-dir "${SCREENING_DIR}/antismash" \
    --output-basename "${SAMPLE_ID}" \
    --cpus ${THREADS} \
    --taxon bacteria \
    --cb-general \
    --clusterhmmer \
    --cb-subclusters \
    --cb-knownclusters \
    --asf \
    --pfam2go \
    "${ASSEMBLY}"

echo "✓ antiSMASH complete"

# Step 2: AMR gene screening with ABRicate
echo ""
echo "Step 2: Screening for antimicrobial resistance genes..."

# Update ABRicate databases
echo "Updating ABRicate databases..."
abricate --setupdb

# Screen against multiple databases
for DB in card ncbi resfinder; do
    echo "  Screening against ${DB} database..."
    abricate \
        --db ${DB} \
        --minid 75 \
        --mincov 50 \
        "${ASSEMBLY}" \
        > "${SCREENING_DIR}/amr/${SAMPLE_ID}_${DB}.tab"
done

# Combine results
echo "FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION" > "${SCREENING_DIR}/amr/${SAMPLE_ID}_amr_all.tab"
cat "${SCREENING_DIR}/amr/"*.tab | grep -v "^#" >> "${SCREENING_DIR}/amr/${SAMPLE_ID}_amr_all.tab" || true

echo "✓ AMR screening complete"

# Step 3: AMRFinderPlus for comprehensive AMR analysis
echo ""
echo "Step 3: Running AMRFinderPlus..."

amrfinder \
    --plus \
    --protein "${PROTEIN_FASTA}" \
    --threads ${THREADS} \
    --output "${SCREENING_DIR}/amr/${SAMPLE_ID}_amrfinder.tsv"

echo "✓ AMRFinderPlus complete"

# Step 5: Generate screening summary
echo ""
echo "Step 5: Generating screening summary..."

SUMMARY_FILE="${SCREENING_DIR}/${SAMPLE_ID}_screening_summary.txt"

cat > "${SUMMARY_FILE}" <<EOF
========================================
Gene Screening Summary for ${SAMPLE_ID}
========================================

EOF

# Count BGCs from antiSMASH
if [ -f "${SCREENING_DIR}/antismash/${SAMPLE_ID}.json" ]; then
    echo "Secondary Metabolite Biosynthetic Gene Clusters:" >> "${SUMMARY_FILE}"
    python3 - << PYTHON_EOF >> "${SUMMARY_FILE}" 2>/dev/null || echo "  Could not parse antiSMASH results" >> "${SUMMARY_FILE}"
import json
try:
    with open("${SCREENING_DIR}/antismash/${SAMPLE_ID}.json", 'r') as f:
        data = json.load(f)
    records = data.get('records', [])
    total_clusters = sum(len(r.get('areas', [])) for r in records)
    print(f"  Total BGCs detected: {total_clusters}")
    
    # Count by type
    cluster_types = {}
    for record in records:
        for area in record.get('areas', []):
            for product in area.get('products', []):
                cluster_types[product] = cluster_types.get(product, 0) + 1
    
    if cluster_types:
        print("  BGC types:")
        for ctype, count in sorted(cluster_types.items()):
            print(f"    - {ctype}: {count}")
except Exception as e:
    print(f"  Error: {e}")
PYTHON_EOF
fi

echo "" >> "${SUMMARY_FILE}"

# Count AMR genes
AMR_COUNT=$(tail -n +2 "${SCREENING_DIR}/amr/${SAMPLE_ID}_amrfinder.tsv" 2>/dev/null | wc -l || echo "0")
cat >> "${SUMMARY_FILE}" <<EOF
Antimicrobial Resistance Genes:
  Total AMR genes detected: ${AMR_COUNT}

EOF

echo "✓ Summary generation complete"

# Display summary
echo ""
echo "========================================="
cat "${SUMMARY_FILE}"
echo ""
echo "Output files:"
echo "  - antiSMASH results: ${SCREENING_DIR}/antismash/index.html"
echo "  - AMR results: ${SCREENING_DIR}/amr/"
echo "  - Summary: ${SUMMARY_FILE}"
echo ""
echo "✓ Gene screening pipeline complete for ${SAMPLE_ID}"
echo "========================================="