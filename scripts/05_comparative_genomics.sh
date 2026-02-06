#!/bin/bash
# Week 2 - Day 8-9: Comparative Genomics and Pangenome Analysis
# Average Nucleotide Identity (ANI), OrthoFinder pangenome analysis,
# core/accessory genome identification, and phylogenomic tree construction

set -e  # Exit on error
set -u  # Exit on undefined variable

# Parse command line arguments
SAMPLE_PATTERN=${1:-""}
THREADS=${2:-8}

if [ -z "$SAMPLE_PATTERN" ]; then
    echo "Usage: $0 <sample_pattern|all> [threads]"
    echo ""
    echo "Examples:"
    echo "  $0 all 8                    # Compare all samples"
    echo "  $0 'F003_*' 8              # Compare all F003 samples"
    echo "  $0 'F003_124 F003_177' 8   # Compare specific samples"
    echo ""
    exit 1
fi

# Get paths
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOTATION_DIR="${BASE_DIR}/results/annotation"
ASSEMBLY_DIR="${BASE_DIR}/results/assembly"
COMPARATIVE_DIR="${BASE_DIR}/results/comparative"

# Create output directories
mkdir -p "${COMPARATIVE_DIR}/ani"
mkdir -p "${COMPARATIVE_DIR}/orthofinder"
mkdir -p "${COMPARATIVE_DIR}/pangenome"

echo "========================================="
echo "Comparative Genomics Pipeline"
echo "========================================="
echo "Sample pattern: ${SAMPLE_PATTERN}"
echo "Threads: ${THREADS}"

# Step 1: Prepare input files for ANI analysis
echo ""
echo "Step 1: Preparing genome assemblies for ANI analysis..."

cd "${COMPARATIVE_DIR}/ani"

# Create a list of polished assemblies based on sample pattern
echo "Creating genome list..."
> genome_list.txt  # Clear the file

if [ "$SAMPLE_PATTERN" = "all" ]; then
    # Find all polished assemblies
    find "${ASSEMBLY_DIR}" -name "*_polished.fasta" >> genome_list.txt
else
    # Find specific samples matching the pattern
    for PATTERN in ${SAMPLE_PATTERN}; do
        # Handle wildcards
        if [[ "$PATTERN" == *"*"* ]]; then
            find "${ASSEMBLY_DIR}" -name "${PATTERN}" -name "*_polished.fasta" >> genome_list.txt
        else
            # Look for specific sample directory
            if [ -d "${ASSEMBLY_DIR}/${PATTERN}" ]; then
                ASSEMBLY="${ASSEMBLY_DIR}/${PATTERN}/${PATTERN}_polished.fasta"
                if [ -f "$ASSEMBLY" ]; then
                    echo "$ASSEMBLY" >> genome_list.txt
                else
                    echo "Warning: Assembly not found for ${PATTERN}"
                fi
            else
                # Try finding it in any directory
                FOUND=$(find "${ASSEMBLY_DIR}" -type f -name "${PATTERN}*_polished.fasta" | head -1)
                if [ -n "$FOUND" ]; then
                    echo "$FOUND" >> genome_list.txt
                else
                    echo "Warning: No assembly found matching ${PATTERN}"
                fi
            fi
        fi
    done
fi

GENOME_COUNT=$(wc -l < genome_list.txt)

if [ "$GENOME_COUNT" -lt 2 ]; then
    echo "Error: Need at least 2 genomes for comparative analysis"
    echo "Found only ${GENOME_COUNT} genome(s)"
    echo ""
    echo "Available assemblies:"
    find "${ASSEMBLY_DIR}" -name "*_polished.fasta" -exec basename {} \; | sed 's/_polished.fasta//'
    exit 1
fi

echo "Found ${GENOME_COUNT} polished genomes for comparison"
echo ""
echo "Genomes to compare:"
cat genome_list.txt | sed 's|.*/||; s|_polished.fasta||' | awk '{print "  - " $0}'
echo ""
echo "✓ Genomes prepared"

# Step 2: Run FastANI (All-vs-All)
echo ""
echo "Step 2: Calculating Average Nucleotide Identity (ANI)..."
echo "ANI ≥ 95% indicates same species"
echo ""

# Run FastANI to calculate pairwise ANI
echo "Running FastANI analysis..."
fastANI \
    --ql genome_list.txt \
    --rl genome_list.txt \
    -o ani_matrix.txt \
    -t ${THREADS} \
    --matrix

echo "✓ FastANI analysis complete"

# Step 3: Create ANI matrix with proper formatting
echo ""
echo "Step 3: Creating formatted ANI matrix..."

# Extract genome names and create header
cut -f1 ani_matrix.txt | sed 's|.*/||; s|_polished.fasta||' | sort -u > genome_names.txt

# Process ANI results into matrix format
awk -F'\t' '
BEGIN {OFS="\t"}
{
    # Extract genome names from paths
    gsub(".*/", "", $1); gsub("_polished.fasta", "", $1);
    gsub(".*/", "", $2); gsub("_polished.fasta", "", $2);
    # Store ANI value
    ani[$1,$2] = $3;
    # Track genome names
    genomes[$1] = 1;
    genomes[$2] = 1;
}
END {
    # Sort genome names
    n = asorti(genomes, sorted);

    # Print header
    printf "Genome";
    for (i=1; i<=n; i++) printf "\t%s", sorted[i];
    printf "\n";

    # Print matrix rows
    for (i=1; i<=n; i++) {
        printf "%s", sorted[i];
        for (j=1; j<=n; j++) {
            if (sorted[i] == sorted[j]) {
                printf "\t100.00";
            } else if ((sorted[i],sorted[j]) in ani) {
                printf "\t%.2f", ani[sorted[i],sorted[j]];
            } else {
                printf "\tNA";
            }
        }
        printf "\n";
    }
}' ani_matrix.txt > ani_matrix.tsv

echo "✓ ANI matrix created"

# Display the matrix
echo ""
echo "ANI Matrix:"
column -t ani_matrix.tsv

# Identify species groups (ANI ≥ 95%)
echo ""
echo "=========================================================="
echo "Species Groups (ANI ≥ 95%):"
echo "=========================================================="
awk -F'\t' '
NR==1 {
    for (i=2; i<=NF; i++) header[i] = $i;
    next;
}
{
    query = $1;
    for (i=2; i<=NF; i++) {
        if ($i != "100.00" && $i != "NA" && $i >= 95) {
            print "  " query " <-> " header[i] ": " $i "%";
        }
    }
}' ani_matrix.tsv | sort -u

# Step 4: Prepare protein sequences for OrthoFinder
echo ""
echo "Step 4: Preparing protein sequences for OrthoFinder..."

cd "${COMPARATIVE_DIR}"
PROTEINS_DIR="${COMPARATIVE_DIR}/orthofinder/proteins"
mkdir -p "${PROTEINS_DIR}"

echo "Copying protein FAA files for selected samples..."
# Extract sample names from genome list
cat "${COMPARATIVE_DIR}/ani/genome_list.txt" | while read GENOME_PATH; do
    SAMPLE=$(basename "$GENOME_PATH" | sed 's/_polished.fasta//')
    FAA_FILE="${ANNOTATION_DIR}/${SAMPLE}/${SAMPLE}.faa"
    
    if [ -f "$FAA_FILE" ]; then
        cp "$FAA_FILE" "${PROTEINS_DIR}/"
        echo "  Copied: ${SAMPLE}.faa"
    else
        echo "  Warning: No protein file found for ${SAMPLE}"
    fi
done

PROTEIN_COUNT=$(ls -1 "${PROTEINS_DIR}"/*.faa 2>/dev/null | wc -l)
if [ "$PROTEIN_COUNT" -lt 2 ]; then
    echo "Error: Need at least 2 protein files for OrthoFinder analysis"
    echo "Found only ${PROTEIN_COUNT} protein file(s)"
    exit 1
fi

echo "✓ Protein sequences prepared (${PROTEIN_COUNT} samples)"

# Step 5: Run OrthoFinder
echo ""
echo "Step 5: Running OrthoFinder for ortholog clustering..."
echo "This may take 30 minutes to 1 hour"
echo ""

orthofinder \
    -f "${PROTEINS_DIR}" \
    -t ${THREADS} \
    -a 4 \
    -M msa \
    -S diamond \
    -o "${COMPARATIVE_DIR}/orthofinder"

echo "✓ OrthoFinder analysis complete"

# Step 6: Navigate to OrthoFinder results
echo ""
echo "Step 6: Processing OrthoFinder results..."

RESULTS_DIR=$(find "${COMPARATIVE_DIR}/orthofinder" -name "Results_*" -type d | head -1)

if [ -z "$RESULTS_DIR" ] || [ ! -d "$RESULTS_DIR" ]; then
    echo "Warning: OrthoFinder results directory not found"
    echo "Expected: ${COMPARATIVE_DIR}/orthofinder/Results_*"
else
    echo "Found results in: $RESULTS_DIR"

    # Copy species tree if it exists
    if [ -f "${RESULTS_DIR}/Species_Tree/SpeciesTree_rooted.txt" ]; then
        cp "${RESULTS_DIR}/Species_Tree/SpeciesTree_rooted.txt" "${COMPARATIVE_DIR}/pangenome/species_tree.nwk"
        echo "✓ Species tree extracted"
    fi
fi

# Step 7: Core and Accessory Genome Analysis
echo ""
echo "Step 7: Analyzing core and accessory genome..."

if [ -z "$RESULTS_DIR" ] || [ ! -d "$RESULTS_DIR" ]; then
    echo "Skipping pangenome analysis - OrthoFinder results not found"
else
    OG_FILE="${RESULTS_DIR}/Orthogroups/Orthogroups.GeneCount.tsv"

    if [ ! -f "$OG_FILE" ]; then
        echo "Warning: Orthogroups file not found at ${OG_FILE}"
    else
        cd "${COMPARATIVE_DIR}/pangenome"

        # Count number of genomes
        N_GENOMES=$(head -1 "$OG_FILE" | awk -F'\t' '{print NF-2}')
        echo "Number of genomes: ${N_GENOMES}"

        # Find core genes (present in ALL genomes)
        echo ""
        echo "Identifying core genes..."
        awk -F'\t' -v n="$N_GENOMES" '
        NR==1 {print; next}
        NR>1 {
            count = 0;
            for(i=2; i<NF; i++) {
                if($i > 0) count++;
            }
            if(count == n) print;
        }' "$OG_FILE" > core_genes.tsv

        CORE_COUNT=$(($(wc -l < core_genes.tsv) - 1))

        # Find accessory genes (in some but not all genomes)
        echo "Identifying accessory genes..."
        awk -F'\t' -v n="$N_GENOMES" '
        NR==1 {print $0 "\tgenome_count"; next}
        NR>1 {
            count = 0;
            for(i=2; i<NF; i++) {
                if($i > 0) count++;
            }
            if(count > 0 && count < n) print $0 "\t" count;
        }' "$OG_FILE" > accessory_genes.tsv

        ACCESSORY_COUNT=$(($(wc -l < accessory_genes.tsv) - 1))

        # Find unique genes (in only ONE genome)
        echo "Identifying unique genes..."
        awk -F'\t' '
        NR==1 {print; next}
        NR>1 {
            count = 0;
            for(i=2; i<NF; i++) {
                if($i > 0) count++;
            }
            if(count == 1) print;
        }' "$OG_FILE" > unique_genes.tsv

        UNIQUE_COUNT=$(($(wc -l < unique_genes.tsv) - 1))
        TOTAL_OG=$(($(wc -l < "$OG_FILE") - 1))

        # Print pangenome summary
        echo ""
        echo "=========================================================="
        echo "PANGENOME SUMMARY"
        echo "=========================================================="
        echo "Total orthogroups: ${TOTAL_OG}"
        echo "Core genes (100%):      ${CORE_COUNT} ($(awk "BEGIN {printf \"%.1f\", 100*${CORE_COUNT}/${TOTAL_OG}")%)"
        echo "Accessory genes:        ${ACCESSORY_COUNT} ($(awk "BEGIN {printf \"%.1f\", 100*${ACCESSORY_COUNT}/${TOTAL_OG}")%)"
        echo "Unique genes (1 genome): ${UNIQUE_COUNT} ($(awk "BEGIN {printf \"%.1f\", 100*${UNIQUE_COUNT}/${TOTAL_OG}")%)"
        echo "=========================================================="

        echo ""
        echo "Saved:"
        echo "  - core_genes.tsv (${CORE_COUNT} orthogroups)"
        echo "  - accessory_genes.tsv (${ACCESSORY_COUNT} orthogroups)"
        echo "  - unique_genes.tsv (${UNIQUE_COUNT} orthogroups)"

        echo "✓ Pangenome analysis complete"
    fi
fi

# Final summary
echo ""
echo "========================================="
echo "Comparative Genomics Summary"
echo "========================================="
echo ""
echo "Output files:"
echo "  - ANI matrix: ${COMPARATIVE_DIR}/ani/ani_matrix.tsv"
echo "  - OrthoFinder results: ${COMPARATIVE_DIR}/orthofinder/"
echo "  - Core genes: ${COMPARATIVE_DIR}/pangenome/core_genes.tsv"
echo "  - Accessory genes: ${COMPARATIVE_DIR}/pangenome/accessory_genes.tsv"
echo "  - Unique genes: ${COMPARATIVE_DIR}/pangenome/unique_genes.tsv"
echo "  - Species tree: ${COMPARATIVE_DIR}/pangenome/species_tree.nwk"
echo ""
echo "Next steps:"
echo "  - Visualize ANI heatmap: See Tutorial 4 for R script"
echo "  - Analyze pangenome composition: See Tutorial 4 for R script"
echo "  - Visualize species tree: See Tutorial 4 for ape/phytools R script"
echo ""
echo "✓ Comparative genomics pipeline complete"
echo "========================================="
