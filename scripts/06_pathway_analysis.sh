#!/bin/bash
# Week 2 - Day 10: Pathway and Metabolic Analysis
# KEGG pathway reconstruction, metabolic profiling, functional enrichment, data visualization

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get paths
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# KofamScan environment variables (set these before running if using KofamScan)
KOFAMSCAN_PROFILES="${KOFAMSCAN_PROFILES:-}"
KOFAMSCAN_KO_LIST="${KOFAMSCAN_KO_LIST:-}"

# Parse command line arguments
SAMPLE_ID=${1:-""}
THREADS=${2:-8}

if [ -z "$SAMPLE_ID" ]; then
    echo "Usage: $0 <sample_id> [threads]"
    echo "       $0 comparative [threads]  # For comparative pathway analysis"
    echo "Example: $0 F003_M_enclense 8"
    exit 1
fi

ANNOTATION_DIR="${BASE_DIR}/results/annotation"
SCREENING_DIR="${BASE_DIR}/results/screening"
PATHWAY_DIR="${BASE_DIR}/results/pathways"

# Create output directories
mkdir -p "${PATHWAY_DIR}/${SAMPLE_ID}"
mkdir -p "${PATHWAY_DIR}/comparative"

echo "========================================="
echo "Pathway Analysis Pipeline for ${SAMPLE_ID}"
echo "========================================="
echo "Threads: ${THREADS}"
echo ""

if [ "${SAMPLE_ID}" == "comparative" ]; then
    echo "Running comparative pathway analysis for all samples..."
    bash "${BASE_DIR}/scripts/helpers/comparative_pathways.sh" "${PATHWAY_DIR}" "${ANNOTATION_DIR}" "${SCREENING_DIR}"
    exit 0
fi

# Validate input files
echo "Checking input files..."

# Find annotation files
TSV_FILE="${ANNOTATION_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.tsv"
PROTEIN_FILE="${ANNOTATION_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.faa"

if [ ! -f "$TSV_FILE" ]; then
    echo "Error: Annotation TSV not found: ${TSV_FILE}"
    echo "Please run 03_annotation.sh first"
    exit 1
fi

echo "Input annotations: ${TSV_FILE}"

# Step 1: Extract EC numbers from Bakta annotation
echo ""
echo "Step 1: Extracting EC numbers and functional annotations..."

# Extract EC numbers
grep -v "^#" "${TSV_FILE}" | awk -F'\t' '$13 != "" {print $13}' | \
    tr ',' '\n' | sort -u > "${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt" || \
    touch "${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt"

EC_COUNT=$(wc -l < "${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt")
echo "  Found ${EC_COUNT} unique EC numbers"

# Extract GO terms
grep -v "^#" "${TSV_FILE}" | awk -F'\t' '$14 != "" {print $14}' | \
    tr ',' '\n' | sort -u > "${PATHWAY_DIR}/${SAMPLE_ID}/go_terms.txt" || \
    touch "${PATHWAY_DIR}/${SAMPLE_ID}/go_terms.txt"

GO_COUNT=$(wc -l < "${PATHWAY_DIR}/${SAMPLE_ID}/go_terms.txt")
echo "  Found ${GO_COUNT} unique GO terms"

# Extract COG categories
grep -v "^#" "${TSV_FILE}" | awk -F'\t' '$15 != "" {print $15}' | \
    tr ',' '\n' | sort -u > "${PATHWAY_DIR}/${SAMPLE_ID}/cog_categories.txt" || \
    touch "${PATHWAY_DIR}/${SAMPLE_ID}/cog_categories.txt"

COG_COUNT=$(wc -l < "${PATHWAY_DIR}/${SAMPLE_ID}/cog_categories.txt")
echo "  Found ${COG_COUNT} unique COG categories"

echo "✓ Functional annotation extraction complete"

# Step 2: KEGG pathway mapping
echo ""
echo "Step 2: Mapping genes to KEGG pathways..."

bash "${BASE_DIR}/scripts/helpers/map_kegg_pathways.sh" \
    "${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt" \
    "${PATHWAY_DIR}/${SAMPLE_ID}" \
    "${SAMPLE_ID}"

echo "✓ KEGG pathway mapping complete"

# Step 3: Analyze coral-relevant pathways
echo ""
echo "Step 3: Analyzing coral-relevant metabolic pathways..."

python3 - <<PYTHON_EOF > "${PATHWAY_DIR}/${SAMPLE_ID}/coral_pathways_summary.txt"
import sys
from pathlib import Path

# Define coral-relevant pathways
CORAL_PATHWAYS = {
    'Nitrogen metabolism': {
        'pathways': ['ko00910'],
        'genes': ['nif', 'nir', 'nor', 'nos', 'nar', 'nap'],
        'importance': 'Nitrogen provision to coral host'
    },
    'Sulfur metabolism': {
        'pathways': ['ko00920'],
        'genes': ['dmd', 'ddd', 'sox', 'cys'],
        'importance': 'DMSP cycling, sulfur compound transformation'
    },
    'Thiamine (Vitamin B1) biosynthesis': {
        'pathways': ['ko00730'],
        'genes': ['thiC', 'thiE', 'thiG', 'thiS'],
        'importance': 'Vitamin provision to coral'
    },
    'Cobalamin (Vitamin B12) biosynthesis': {
        'pathways': ['ko00860'],
        'genes': ['cob', 'cbi', 'hem'],
        'importance': 'Essential cofactor provision'
    },
    'Oxidative phosphorylation': {
        'pathways': ['ko00190'],
        'genes': ['atp', 'ndh', 'cox', 'cyc'],
        'importance': 'Energy metabolism'
    },
    'Oxidative stress response': {
        'pathways': ['ko00480', 'ko00250'],
        'genes': ['cat', 'sod', 'gpx', 'ahp', 'trx'],
        'importance': 'Protection against oxidative damage'
    }
}

print("=" * 60)
print(f"Coral-Relevant Pathway Analysis: ${SAMPLE_ID}")
print("=" * 60)
print()

# Read EC numbers
ec_file = Path("${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt")
ec_numbers = set()
if ec_file.exists():
    with open(ec_file) as f:
        ec_numbers = set(line.strip() for line in f if line.strip())
    print(f"Total EC numbers detected: {len(ec_numbers)}")
    print()

# Analyze each pathway category
for pathway_name, details in CORAL_PATHWAYS.items():
    print(f"{pathway_name}:")
    print(f"  Importance: {details['importance']}")
    
    # Check for gene presence (basic matching on gene name substrings)
    genes_found = [g for g in details['genes'] if any(g.lower() in ec.lower() for ec in ec_numbers)]
    if genes_found:
        print(f"  Detected genes: {', '.join(genes_found)}")
    else:
        print(f"  Key genes to look for: {', '.join(details['genes'])}")
    print()

print()
print("Note: Detailed pathway completeness requires KEGG Mapper analysis")
print("Visit: https://www.genome.jp/kegg/mapper/")
PYTHON_EOF

cat "${PATHWAY_DIR}/${SAMPLE_ID}/coral_pathways_summary.txt"

echo "✓ Coral pathway analysis complete"

# Step 4: Create pathway visualization
echo ""
echo "Step 4: Creating pathway visualizations..."

python3 - <<PYTHON_EOF
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

# Read annotations
try:
    df = pd.read_csv("${TSV_FILE}", sep='\t', comment='#', header=0)
    
    # COG functional categories
    cog_data = []
    for idx, row in df.iterrows():
        cog = row.get('COG', '')
        if pd.notna(cog) and cog:
            cog_data.extend(cog.split(','))
    
    if cog_data:
        cog_counts = Counter(cog_data)
        
        # COG category names (simplified)
        cog_names = {
            'J': 'Translation', 'K': 'Transcription', 'L': 'Replication',
            'D': 'Cell division', 'V': 'Defense', 'T': 'Signal transduction',
            'M': 'Cell wall', 'N': 'Motility', 'U': 'Secretion',
            'O': 'Chaperones', 'C': 'Energy', 'G': 'Carbohydrate',
            'E': 'Amino acid', 'F': 'Nucleotide', 'H': 'Coenzyme',
            'I': 'Lipid', 'P': 'Inorganic ion', 'Q': 'Secondary metabolites',
            'R': 'General function', 'S': 'Unknown function'
        }
        
        # Plot COG distribution
        fig, ax = plt.subplots(figsize=(12, 6))
        
        sorted_cogs = sorted(cog_counts.items(), key=lambda x: x[1], reverse=True)[:15]
        categories = [cog_names.get(c[0], c[0]) for c in sorted_cogs]
        counts = [c[1] for c in sorted_cogs]
        
        ax.barh(categories, counts, color='steelblue')
        ax.set_xlabel('Number of Genes', fontsize=12)
        ax.set_ylabel('COG Category', fontsize=12)
        ax.set_title('COG Functional Category Distribution', fontsize=14, pad=20)
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        plt.savefig("${PATHWAY_DIR}/${SAMPLE_ID}/cog_distribution.png", 
                   dpi=300, bbox_inches='tight')
        print("✓ Saved COG distribution plot")
    
except Exception as e:
    print(f"Warning: Could not create COG visualization: {e}")

# Create pathway presence/absence summary
try:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Read gene screening results if available
    screening_file = "${SCREENING_DIR}/${SAMPLE_ID}/custom_genes/summary.txt"
    import os
    if os.path.exists(screening_file):
        with open(screening_file) as f:
            content = f.read()
        
        ax.text(0.1, 0.9, 'Coral-Associated Traits:', 
               fontsize=16, weight='bold', transform=ax.transAxes)
        ax.text(0.1, 0.1, content, 
               fontsize=10, family='monospace', 
               transform=ax.transAxes, verticalalignment='bottom')
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig("${PATHWAY_DIR}/${SAMPLE_ID}/traits_summary.png", 
                   dpi=300, bbox_inches='tight')
        print("✓ Saved traits summary plot")
    
except Exception as e:
    print(f"Warning: Could not create traits summary: {e}")

PYTHON_EOF

echo "✓ Pathway visualizations complete"

# Step 5: Generate comprehensive summary
echo ""
echo "Step 5: Generating comprehensive pathway summary..."

SUMMARY_FILE="${PATHWAY_DIR}/${SAMPLE_ID}/pathway_summary.txt"

cat > "${SUMMARY_FILE}" <<EOF
========================================
Pathway Analysis Summary for ${SAMPLE_ID}
========================================

Functional Annotation Statistics:
  EC numbers: ${EC_COUNT}
  GO terms: ${GO_COUNT}
  COG categories: ${COG_COUNT}

EOF

# Add coral pathway summary
cat "${PATHWAY_DIR}/${SAMPLE_ID}/coral_pathways_summary.txt" >> "${SUMMARY_FILE}"

echo "✓ Summary generation complete"

# Display summary
echo ""
echo "========================================="
cat "${SUMMARY_FILE}"
echo ""
echo "Output files:"
echo "  - EC numbers: ${PATHWAY_DIR}/${SAMPLE_ID}/ec_numbers.txt"
echo "  - GO terms: ${PATHWAY_DIR}/${SAMPLE_ID}/go_terms.txt"
echo "  - COG categories: ${PATHWAY_DIR}/${SAMPLE_ID}/cog_categories.txt"
echo "  - KEGG pathways: ${PATHWAY_DIR}/${SAMPLE_ID}/kegg_pathways.txt"
echo "  - Coral pathway summary: ${PATHWAY_DIR}/${SAMPLE_ID}/coral_pathways_summary.txt"
echo "  - Visualizations: ${PATHWAY_DIR}/${SAMPLE_ID}/*.png"
echo "  - Summary: ${SUMMARY_FILE}"
echo ""
echo "Next steps:"
echo "  - Upload EC numbers to KEGG Mapper: https://www.genome.jp/kegg/mapper/"
echo "  - Run comparative pathway analysis: $0 comparative"
echo ""

# Step 6: KofamScan analysis for KEGG Orthology assignment ##### untested
echo ""
echo "Step 6: Running KofamScan for KEGG Orthology assignment..."

KOFAMSCAN_DIR="${PATHWAY_DIR}/${SAMPLE_ID}/kofamscan"
mkdir -p "${KOFAMSCAN_DIR}"

if ! command -v exec_annotation &> /dev/null; then
    echo "Warning: KofamScan not found in PATH. Skipping KofamScan analysis."
    echo "  Install from: https://github.com/takaram/kofam_scan"
else
    if [ -z "${KOFAMSCAN_PROFILES}" ] || [ -z "${KOFAMSCAN_KO_LIST}" ]; then
        echo "Warning: KOFAMSCAN_PROFILES or KOFAMSCAN_KO_LIST not set."
        echo "  Set environment variables before running:"
        echo "    export KOFAMSCAN_PROFILES=/path/to/profiles"
        echo "    export KOFAMSCAN_KO_LIST=/path/to/ko_list"
        echo "  Skipping KofamScan analysis"
    elif [ ! -f "$PROTEIN_FILE" ]; then
        echo "Warning: Protein file not found: ${PROTEIN_FILE}"
        echo "  Skipping KofamScan analysis"
    else
        # Run KofamScan
        echo "Running KofamScan on ${PROTEIN_FILE}..."
        
        exec_annotation \
            -p "${KOFAMSCAN_PROFILES}" \
            -k "${KOFAMSCAN_KO_LIST}" \
            -o "${KOFAMSCAN_DIR}/kofamscan_results.txt" \
            --cpu "${THREADS}" \
            --format detail-tsv \
            "${PROTEIN_FILE}" 2>&1 | tee "${KOFAMSCAN_DIR}/kofamscan.log"
        
        # Parse KofamScan results
        if [ -f "${KOFAMSCAN_DIR}/kofamscan_results.txt" ]; then
            # Extract KO assignments (threshold met)
            grep -v "^#" "${KOFAMSCAN_DIR}/kofamscan_results.txt" | \
                awk '$1 == "*" {print $3}' | sort -u > "${KOFAMSCAN_DIR}/ko_list.txt" || \
                touch "${KOFAMSCAN_DIR}/ko_list.txt"
            
            KO_COUNT=$(wc -l < "${KOFAMSCAN_DIR}/ko_list.txt")
            echo "  Found ${KO_COUNT} KEGG Orthology assignments"
            
            # Create KO summary
            cat > "${KOFAMSCAN_DIR}/ko_summary.txt" <<KOEOF

KofamScan Analysis Summary
==========================

Total KO assignments (threshold met): ${KO_COUNT}

Top 20 assigned KOs:
$(grep -v "^#" "${KOFAMSCAN_DIR}/kofamscan_results.txt" | \
  awk '$1 == "*" {print $3"\t"$7}' | \
  cut -f1 | sort | uniq -c | sort -rn | head -20)

KOEOF
            
            cat "${KOFAMSCAN_DIR}/ko_summary.txt"
            
            echo "✓ KofamScan analysis complete"
            echo "  Results: ${KOFAMSCAN_DIR}/kofamscan_results.txt"
            echo "  KO list: ${KOFAMSCAN_DIR}/ko_list.txt"
        else
            echo "Warning: KofamScan did not produce output file"
        fi
    fi
fi

echo ""
echo "========================================="
echo "✓ Pathway analysis pipeline complete for ${SAMPLE_ID}"
echo "========================================="