#!/bin/bash
set -euo pipefail

# Usage: script.sh <matrix_file> <pattern_to_remove>
# Example: script.sh ani_matrix.txt polished.fa

MATRIX_FILE="${1:-ani_matrix.txt}"
PATTERN="${2:-polished.fa}"

if [[ ! -f "$MATRIX_FILE" ]]; then
    echo "Error: Matrix file '$MATRIX_FILE' not found." >&2
    exit 1
fi

# Step 1: Extract genome names and create header
echo "Creating formatted ANI matrix..."

# Get unique genome names
cut -f1 "$MATRIX_FILE" | sed "s|.*/||; s|_${PATTERN}||" | sort -u > genome_names.txt

# Step 2: Process ANI results into a matrix format
awk -v pattern="$PATTERN" -F'\t' '
BEGIN {OFS="\t"}
{
    # Extract genome names from paths
    gsub(".*/", "", $1); gsub("_" pattern, "", $1);
    gsub(".*/", "", $2); gsub("_" pattern, "", $2);
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
}' "$MATRIX_FILE" > ani_matrix_formatted.tsv

# Step 3-5: Display and analyze
echo ""
echo "ANI Matrix:"
column -t ani_matrix_formatted.tsv

echo ""
echo "=========================================================="
echo "ANI Interpretation Guide:"
echo "=========================================================="
echo "ANI ≥ 95%  : Same species"
echo "ANI 90-95% : Closely related, possibly same genus"
echo "ANI < 90%  : Different genera"
echo "=========================================================="

echo ""
echo "Species Groups (ANI ≥ 95%):"
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
}' ani_matrix_formatted.tsv | sort -u