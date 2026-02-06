# Tutorial 4: Comparative Genomics Analysis

## Learning Objectives

By the end of this tutorial, you will be able to:
- Calculate genome similarity using Average Nucleotide Identity (ANI)
- Perform pangenome analysis to identify core and accessory genes
- Use OrthoFinder to cluster orthologous genes
- Construct phylogenomic trees from core genes
- Interpret comparative genomics results
- Understand species boundaries and genome evolution

---

## Prerequisites

### Required Software

```bash
# Activate comparative genomics environment
conda activate VTK2026_Comparative

# Verify tools are installed
fastANI --version
orthofinder --help | head -5
```

### Required Input

You must have completed Tutorial 3 (Annotation). You need:
```bash
# Polished assemblies from assembly step
ls results/assembly/*/F*_polished.fasta

# Protein sequences from annotation step
ls results/annotation/*/*.faa
```

### Directory Setup

```bash
# Navigate to project directory
cd ~/'MY-NAME-EXAMPLE'/

# Create comparative genomics directory
mkdir -p results/comparative
mkdir -p results/comparative/ani
mkdir -p results/comparative/orthofinder
mkdir -p results/comparative/pangenome
```

---

## Background: Comparative Genomics

### What is Comparative Genomics?

**Comparative genomics** analyzes similarities and differences between genomes to understand:
- **Evolutionary relationships** - How closely related are the organisms?
- **Species boundaries** - Are they the same species or different?
- **Gene content** - Which genes are shared? Which are unique?
- **Functional diversity** - What unique capabilities does each genome have?

### Key Concepts

#### 1. Average Nucleotide Identity (ANI)

**ANI** measures overall genome similarity by comparing nucleotide sequences.

**Interpretation:**
- **ANI ≥ 95%**: Same species
- **ANI 90-95%**: Closely related, possibly same genus
- **ANI < 90%**: Different genera

**Example:**
- Two *Microbacterium enclense* strains: ANI ~98-99%
- *M. enclense* vs *M. ginsengisoli*: ANI ~85-92%
- *Microbacterium* vs *Alteromonas*: ANI ~70-75%

#### 2. Pangenome Analysis

The **pangenome** is the entire gene repertoire of a group of organisms.

**Components:**
- **Core genome**: Genes present in ALL genomes (essential functions)
- **Accessory genome**: Genes in some but not all genomes (niche adaptation)
- **Unique genes**: Genes found in only one genome (strain-specific)

**Example pangenome (3 genomes):**
```
Genome A: 3,500 genes (3,000 core + 300 accessory + 200 unique)
Genome B: 3,600 genes (3,000 core + 400 accessory + 200 unique)
Genome C: 3,400 genes (3,000 core + 250 accessory + 150 unique)

Pangenome: 4,450 genes total
  Core: 3,000 genes (67%)
  Accessory: 950 genes (21%)
  Unique: 550 genes (12%)
```

#### 3. Orthologous Genes

**Orthologs** are genes in different species that evolved from a common ancestor.

**Types:**
- **1:1 orthologs**: One gene in each genome (typical for core genes)
- **Many-to-many orthologs**: Gene families with duplications
- **Paralogs**: Genes duplicated within the same genome

---

## Part 1: Average Nucleotide Identity (ANI)

### Step 1.1: Prepare Input Files

ANI analysis requires genome assemblies (FASTA format).

```bash
# Create file list for ANI analysis
cd results/comparative/ani

# Create a list of all polished assemblies
find ../../assembly -name "*_polished.fasta" > genome_list.txt

# View the list
cat genome_list.txt
```

**Expected output:**
```
../../assembly/F003_M_enclense/F003_M_enclense_polished.fasta
../../assembly/F003_M_ginsengisoli/F003_M_ginsengisoli_polished.fasta
../../assembly/F003_A_portus/F003_A_portus_polished.fasta
../../assembly/H2_M_enclense/H2_M_enclense_polished.fasta
../../assembly/H2_M_ginsengisoli/H2_M_ginsengisoli_polished.fasta
../../assembly/H2_A_portus/H2_A_portus_polished.fasta
```

### Step 1.2: Run FastANI (All-vs-All)

```bash
# Run FastANI to calculate pairwise ANI
fastANI \
    --ql genome_list.txt \
    --rl genome_list.txt \
    -o ani_matrix.txt \
    -t 8 \
    --matrix
```

**Parameter explanation:**
- `--ql`: Query genome list
- `--rl`: Reference genome list
- `-o`: Output file
- `-t`: Number of threads
- `--matrix`: Generate matrix output

**⏱️ Runtime:** ~5-10 minutes for 6 genomes

### Step 1.3: View ANI Results

```bash
# View raw ANI output
cat ani_matrix.txt
```

**Output format (tab-separated):**
```
query_genome    reference_genome    ANI     fragments_aligned    total_fragments
```

### Step 1.4: Create ANI Matrix with bash

```bash
# Step 1: Extract genome names and create header
echo "Creating formatted ANI matrix..."

# Get unique genome names
cut -f1 ani_matrix.txt | sed 's|.*/||; s|_polished.fasta||' | sort -u > genome_names.txt

# Step 2: Process ANI results into a matrix format
# This creates a tab-delimited matrix
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
}' ani_matrix.txt > ani_matrix_formatted.tsv

# Step 3: Display the matrix
echo ""
echo "ANI Matrix:"
column -t ani_matrix_formatted.tsv

# Step 4: Print interpretation guide
echo ""
echo "=========================================================="
echo "ANI Interpretation Guide:"
echo "=========================================================="
echo "ANI ≥ 95%  : Same species"
echo "ANI 90-95% : Closely related, possibly same genus"
echo "ANI < 90%  : Different genera"
echo "=========================================================="

# Step 5: Identify species groups (ANI ≥ 95%)
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
```

**Example output:**
```
                      F003_M_enclense  F003_M_ginsengisoli  F003_A_portus  ...
F003_M_enclense            100.00          87.45            72.13
F003_M_ginsengisoli         87.45         100.00            71.89
F003_A_portus               72.13          71.89           100.00
H2_M_enclense               98.76          87.32            72.08
H2_M_ginsengisoli           87.38          98.45            71.95
H2_A_portus                 72.19          71.94            99.12

Species Groups (ANI ≥ 95%):
  F003_M_enclense: H2_M_enclense
  F003_M_ginsengisoli: H2_M_ginsengisoli
  F003_A_portus: H2_A_portus
  H2_M_enclense: F003_M_enclense
  H2_M_ginsengisoli: F003_M_ginsengisoli
  H2_A_portus: F003_A_portus
```

### Step 1.5: Visualize ANI with Heatmap in R

```R
# Create ANI heatmap and dendrogram using R
# Suggestion: make an plot_ani_heatmap.R 
#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Read ANI data
data <- read.table('results/comparative/ani/ani_matrix.tsv', header=FALSE, sep='\t')
colnames(data) <- c('Query', 'Reference', 'ANI', 'Matches', 'Fragments')

# Convert to matrix format
ani_matrix <- xtabs(data[[3]] ~ data[[1]] + data[[2]])

# Convert to numeric matrix
ani_mat <- as.matrix(ani_matrix)

# Create heatmap with pheatmap
png('ani_heatmap.png', width=15, height=15, units='in', res=300)
pheatmap(ani_mat,
         display_numbers = TRUE,
         number_format = "%.2f",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("red", "yellow", "green"))(100),
         breaks = seq(70, 100, length.out=101),
         main = "Average Nucleotide Identity (ANI) Matrix",
         fontsize = 10,
         fontsize_number = 8,
         angle_col = 45)
dev.off()
cat("Heatmap saved: ani_heatmap.png\n")

# Create dendrogram
# Convert ANI to distance (100 - ANI)
dist_matrix <- 100 - ani_mat
# Ensure symmetry
dist_matrix <- (dist_matrix + t(dist_matrix)) / 2
diag(dist_matrix) <- 0

# Convert to distance object (lower triangle only)
dist_obj <- as.dist(dist_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_obj, method="average")

# Plot dendrogram
png('ani_dendrogram.png', width=10, height=6, units='in', res=300)
plot(hc, 
     main="Hierarchical Clustering Based on ANI",
     xlab="Genome",
     ylab="Distance (100 - ANI)",
     cex=0.8)
dev.off()
cat("Dendrogram saved: ani_dendrogram.png\n")

# Run the R script
#Rscript plot_ani_heatmap.R
```

---

## Part 2: Pangenome Analysis with OrthoFinder

### Step 2.1: Prepare Protein Sequences

OrthoFinder requires protein sequences from all genomes.

```bash
# Navigate to OrthoFinder directory
cd ../orthofinder

# Create directory for input proteins
mkdir -p proteins

# Copy all protein FAA files
for ANNOT_DIR in ../../annotation/*; do
    SAMPLE=$(basename $ANNOT_DIR)
    if [ -f "$ANNOT_DIR/${SAMPLE}.faa" ]; then
        cp "$ANNOT_DIR/${SAMPLE}.faa" proteins/
        echo "Copied: ${SAMPLE}.faa"
    fi
done

# Verify files
ls -lh proteins/
```

### Step 2.2: Run OrthoFinder

```bash
# Run OrthoFinder
orthofinder \
    -f proteins/ \
    -t 8 \
    -a 4 \
    -M msa \
    -S diamond \
    -o orthofinder_results
```

**Parameter explanation:**
- `-f`: Input directory with protein FAA files
- `-t`: Number of threads for sequence search
- `-a`: Number of threads for alignment
- `-M msa`: Use multiple sequence alignment for tree inference
- `-S diamond`: Use DIAMOND for faster sequence searches
- `-o`: Output directory

**⏱️ Runtime:** ~30-60 minutes for 6 genomes

**OrthoFinder steps:**
1. **DIAMOND search**: All-vs-all protein sequence comparison
2. **Clustering**: Group proteins into orthogroups
3. **Multiple alignment**: Align each orthogroup
4. **Gene trees**: Build phylogenetic tree for each orthogroup
5. **Species tree**: Infer species relationships from gene trees

### Step 2.3: Navigate Results

```bash
# Navigate to results directory
RESULTS_DIR=$(find orthofinder_results -name "Results_*" -type d)
cd $RESULTS_DIR

# List output files
ls -lh
```

**Key output files:**

| File | Description |
|------|-------------|
| `Orthogroups/Orthogroups.tsv` | Orthogroups (gene families) |
| `Orthogroups/Orthogroups.GeneCount.tsv` | Gene counts per orthogroup |
| `Phylogenetic_Hierarchical_Orthogroups/` | Hierarchical orthogroups |
| `Species_Tree/SpeciesTree_rooted.txt` | Species phylogeny (Newick) |
| `Comparative_Genomics_Statistics/` | Summary statistics |
| `Orthogroup_Sequences/` | Aligned protein sequences |
| `Single_Copy_Orthologue_Sequences/` | Core single-copy genes |

### Step 2.4: View Orthogroups

```bash
# View orthogroups table
head -20 Orthogroups/Orthogroups.tsv | column -t

# Count orthogroups
TOTAL_OG=$(wc -l < Orthogroups/Orthogroups.tsv)
echo "Total orthogroups: $((TOTAL_OG - 1))"  # Subtract header

# Count genes in each genome
cut -f2- Orthogroups/Orthogroups.tsv | head -1
```

**Example output:**
```
Orthogroup  F003_M_enclense         F003_M_ginsengisoli    F003_A_portus
OG0000000   dnaA, dnaN, dnaB        dnaA, dnaN, dnaB       dnaA, dnaN, dnaB
OG0000001   rpoA, rpoB, rpoC        rpoA, rpoB, rpoC       rpoA, rpoB, rpoC
OG0000002   hypothetical_001        hypothetical_002       -
```

### Step 2.5: View Gene Count Matrix

```bash
# View gene count summary
cat Orthogroups/Orthogroups.GeneCount.tsv | column -t

# Calculate statistics
echo ""
echo "Orthogroup Statistics:"
awk -F'\t' 'NR>1 {
    total = 0;
    for(i=2; i<=NF; i++) total += $i;
    if(total > 0) count++;
}
END {print "Total orthogroups:", count}' Orthogroups/Orthogroups.GeneCount.tsv
```

### Step 2.6: View Comparative Statistics

```bash
# View comprehensive statistics
cat Comparative_Genomics_Statistics/Statistics_Overall.tsv

# View per-species statistics
cat Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | column -t
```

**Statistics include:**
- Number of genes per genome
- Number of orthogroups
- Number of species-specific orthogroups
- Percentage of genes in orthogroups

---

## Part 3: Core and Accessory Genome Analysis

### Step 3.1: Define Core Genome with bash

**Core genome** = genes present in ALL genomes

```bash
# Navigate to pangenome directory
cd ../../../pangenome

# Find the OrthoFinder results directory
OG_FILE=$(find ../orthofinder/orthofinder_results -name "Orthogroups.GeneCount.tsv" | head -1)

echo "Using file: $OG_FILE"

# Step 1: Count number of genomes (exclude Orthogroup and Total columns)
N_GENOMES=$(head -1 "$OG_FILE" | awk -F'\t' '{print NF-2}')
echo "Number of genomes: $N_GENOMES"

# Step 2: Display genome names
echo ""
echo "Genomes:"
head -1 "$OG_FILE" | awk -F'\t' '{for(i=2; i<NF; i++) print "  - " $i}'

# Step 3: Identify CORE genes (present in ALL genomes)
echo ""
echo "Identifying core genes..."
awk -F'\t' -v n="$N_GENOMES" '
NR==1 {print; next}  # Print header
NR>1 {
    count = 0;
    for(i=2; i<NF; i++) {  # Skip Orthogroup and Total columns
        if($i > 0) count++;
    }
    if(count == n) print;
}' "$OG_FILE" > core_genes.tsv

CORE_COUNT=$(($(wc -l < core_genes.tsv) - 1))

# Step 4: Identify ACCESSORY genes (in some but not all genomes)
echo "Identifying accessory genes..."
awk -F'\t' -v n="$N_GENOMES" '
NR==1 {print $0 "\tgenome_count"; next}  # Add genome_count column
NR>1 {
    count = 0;
    for(i=2; i<NF; i++) {
        if($i > 0) count++;
    }
    if(count > 0 && count < n) print $0 "\t" count;
}' "$OG_FILE" > accessory_genes.tsv

ACCESSORY_COUNT=$(($(wc -l < accessory_genes.tsv) - 1))

# Step 5: Identify UNIQUE genes (in only ONE genome)
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

# Step 6: Calculate total orthogroups
TOTAL_OG=$(($(wc -l < "$OG_FILE") - 1))

# Step 7: Print summary
echo ""
echo "=========================================================="
echo "PANGENOME SUMMARY"
echo "=========================================================="
echo "Total orthogroups: $TOTAL_OG"
echo "Core genes (100%):      $CORE_COUNT ($(awk "BEGIN {printf \"%.1f\", 100*$CORE_COUNT/$TOTAL_OG}")%)"
echo "Accessory genes:        $ACCESSORY_COUNT ($(awk "BEGIN {printf \"%.1f\", 100*$ACCESSORY_COUNT/$TOTAL_OG}")%)"
echo "Unique genes (1 genome): $UNIQUE_COUNT ($(awk "BEGIN {printf \"%.1f\", 100*$UNIQUE_COUNT/$TOTAL_OG}")%)"
echo "=========================================================="

echo ""
echo "Saved:"
echo "  - core_genes.tsv ($CORE_COUNT orthogroups)"
echo "  - accessory_genes.tsv ($ACCESSORY_COUNT orthogroups)"
echo "  - unique_genes.tsv ($UNIQUE_COUNT orthogroups)"

# Step 8: Count unique genes per genome
echo ""
echo "=========================================================="
echo "UNIQUE GENES PER GENOME"
echo "=========================================================="

# Get genome names from header
head -1 "$OG_FILE" | awk -F'\t' '{for(i=2; i<NF; i++) print $i}' | while read GENOME; do
    # Count unique genes for this genome
    GENOME_COL=$(head -1 "$OG_FILE" | awk -F'\t' -v g="$GENOME" '{for(i=1; i<=NF; i++) if($i==g) print i}')
    
    GENOME_UNIQUE=$(awk -F'\t' -v col="$GENOME_COL" '
    NR>1 {
        count = 0;
        for(i=2; i<NF; i++) {
            if($i > 0) count++;
        }
        if(count == 1 && $col > 0) print;
    }' "$OG_FILE" | wc -l)
    
    echo "$GENOME: $GENOME_UNIQUE unique orthogroups"
    
    # Save genome-specific unique genes
    awk -F'\t' -v col="$GENOME_COL" -v g="$GENOME" '
    NR==1 {print "Orthogroup\t" g}
    NR>1 {
        count = 0;
        for(i=2; i<NF; i++) {
            if($i > 0) count++;
        }
        if(count == 1 && $col > 0) print $1 "\t" $col;
    }' "$OG_FILE" > "unique_genes_${GENOME}.tsv"
done
```

### Step 3.2: Visualize Pangenome with R


```R
# Create pangenome visualization using R
# Suggestion: make an R script plot_pangenome.R
#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

# Find and read orthogroups file dynamically
results_dir <- system("find results/comparative/orthofinder/orthofinder_results -name 'Results_*' -type d", intern=TRUE)
if (length(results_dir) > 1) {
  results_dir <- results_dir[length(results_dir)]  # Use most recent
}
og_file <- file.path(results_dir, "Orthogroups", "Orthogroups.GeneCount.tsv")
df <- read.table(og_file, header=TRUE, sep='\t', row.names=1)

# Get genome columns (exclude Total if present)
if("Total" %in% colnames(df)) {
  genome_cols <- colnames(df)[colnames(df) != "Total"]
} else {
  genome_cols <- colnames(df)
}
n_genomes <- length(genome_cols)

cat("Number of genomes:", n_genomes, "\n")

# Count genes by presence across genomes
df$genome_count <- rowSums(df[, genome_cols] > 0)

# Count orthogroups by genome presence
presence_counts <- table(df$genome_count)
presence_df <- data.frame(n_genomes_present = as.numeric(names(presence_counts)),
                          n_orthogroups = as.numeric(presence_counts))

# Define categories
core_count <- sum(presence_counts[names(presence_counts) == n_genomes])
accessory_count <- sum(presence_counts[names(presence_counts) != "1" & 
                                        as.numeric(names(presence_counts)) < n_genomes])
unique_count <- sum(presence_counts[names(presence_counts) == "1"])

# Plot 1: Gene presence distribution
p1 <- ggplot(presence_df, aes(x=n_genomes_present, y=n_orthogroups)) +
  geom_bar(stat="identity", aes(fill=factor(ifelse(n_genomes_present==1, "Unique",
                                              ifelse(n_genomes_present==n_genomes, "Core", "Accessory"))))) +
  scale_fill_manual(values=c("Core"="#2ca02c", "Accessory"="#ff7f0e", "Unique"="#d62728"),
                    name="Category") +
  labs(title="Gene Presence Distribution",
       x="Number of Genomes",
       y="Number of Orthogroups") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", hjust=0.5))

# Plot 2: Pangenome composition pie chart
comp_df <- data.frame(
  category = c("Core", "Accessory", "Unique"),
  count = c(core_count, accessory_count, unique_count)
)
comp_df$percentage <- round(100 * comp_df$count / sum(comp_df$count), 1)
comp_df$label <- paste0(comp_df$category, "\n(", comp_df$count, ")\n", comp_df$percentage, "%")

p2 <- ggplot(comp_df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("Core"="#2ca02c", "Accessory"="#ff7f0e", "Unique"="#d62728")) +
  labs(title="Pangenome Composition") +
  theme_void() +
  theme(plot.title = element_text(face="bold", hjust=0.5)) +
  geom_text(aes(label=label), position=position_stack(vjust=0.5), size=3)

# Plot 3: Total genes per genome
genes_per_genome <- colSums(df[, genome_cols] > 0)
genes_df <- data.frame(genome = names(genes_per_genome),
                       n_genes = as.numeric(genes_per_genome))

p3 <- ggplot(genes_df, aes(x=reorder(genome, n_genes), y=n_genes)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=n_genes), hjust=-0.2, size=3) +
  coord_flip() +
  labs(title="Total Genes per Genome",
       x="Genome",
       y="Number of Genes") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", hjust=0.5))

# Plot 4: Unique genes per genome
unique_per_genome <- sapply(genome_cols, function(g) {
  sum(df[, g] > 0 & df$genome_count == 1)
})
unique_df <- data.frame(genome = names(unique_per_genome),
                        n_unique = as.numeric(unique_per_genome))

p4 <- ggplot(unique_df, aes(x=reorder(genome, n_unique), y=n_unique)) +
  geom_bar(stat="identity", fill="coral") +
  geom_text(aes(label=n_unique), hjust=-0.2, size=3) +
  coord_flip() +
  labs(title="Genome-Specific Genes",
       x="Genome",
       y="Number of Unique Orthogroups") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", hjust=0.5))

# Combine all plots
png('pangenome_analysis.png', width=14, height=10, units='in', res=300)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

cat("Pangenome plots saved: pangenome_analysis.png\n")

# Run the R script
# Rscript plot_pangenome.R
```

### Step 3.3: Analyze Core Genes

```bash
# View core genes
head -20 core_genes.tsv | column -t

# Count core genes
CORE_COUNT=$(wc -l < core_genes.tsv)
echo "Core genes: $((CORE_COUNT - 1))"

# Extract core gene functions (requires linking back to annotations)
# This is simplified - in practice, you'd extract gene names/products
echo ""
echo "Core genes are typically involved in:"
echo "  - DNA replication and repair"
echo "  - Transcription and translation"
echo "  - Central metabolism"
echo "  - Cell division"
echo "  - Ribosomal proteins"
```

---

## Part 4: Phylogenomic Analysis

### Step 4.1: View Species Tree

OrthoFinder generates a species tree based on core genes.

```bash
# Navigate back to OrthoFinder results
cd ../orthofinder
RESULTS_DIR=$(find orthofinder_results -name "Results_*" -type d)
cd $RESULTS_DIR/Species_Tree

# View Newick tree
cat SpeciesTree_rooted.txt
```

**Newick format example:**
```
((F003_M_enclense:0.02,H2_M_enclense:0.02):0.15,(F003_M_ginsengisoli:0.02,H2_M_ginsengisoli:0.02):0.15,(F003_A_portus:0.25,H2_A_portus:0.25):0.40);
```

### Step 4.2: Visualize Phylogenetic Tree with R

```R
# Suggestion: make an R script plot_tree.R
# Create tree visualization using R
#!/usr/bin/env Rscript

# Load required library
library(ape)
library(phytools)

# Find and read tree file
results_dir <- system("find results/comparative/orthofinder/orthofinder_results -name 'Results_*' -type d", intern=TRUE)
if (length(results_dir) > 1) {
  results_dir <- results_dir[length(results_dir)]  # Use most recent
}
tree_file <- file.path(results_dir, "Species_Tree", "SpeciesTree_rooted.txt")

cat("Reading tree from:", tree_file, "\n")
tree <- read.tree(tree_file)

# Create tree visualization
png('species_tree.png', width=10, height=8, units='in', res=300)
par(mar=c(4,2,3,2))
plot(tree, 
     main="Species Phylogenetic Tree\n(Based on Core Gene Orthologs)",
     cex.main=1.2, font.main=2,
     edge.width=2,
     cex=0.9,
     no.margin=FALSE)
add.scale.bar(cex=0.8, lwd=2)
axisPhylo(side=1, cex.axis=0.8)
mtext("Evolutionary Distance", side=1, line=2.5, cex=0.9)
dev.off()

cat("Tree saved: species_tree.png\n")

# Print tree structure
cat("\nTree structure:\n")
cat("Number of tips:", Ntip(tree), "\n")
cat("Tip labels:", paste(tree$tip.label, collapse=", "), "\n")
cat("\nNewick format:\n")
write.tree(tree, "")
cat("\n")

# Run the R script
# Rscript plot_tree.R
```

### Step 4.3: Interpret Tree

The phylogenetic tree shows evolutionary relationships:

1. **Branch lengths**: Longer = more divergence
2. **Clades**: Groups of closely related organisms
3. **Monophyletic groups**: Species from same lineage cluster together

**Expected pattern:**
```
- M. enclense strains (F003 + H2) cluster together
- M. ginsengisoli strains cluster together  
- Alteromonas strains cluster together
- Microbacterium species separate from Alteromonas
```

---

## Part 5: Comprehensive Summary Report
<!-- NOT WORKING
### Step 5.1: Generate Combined Report

```bash
# Create comprehensive summary
cd ../comparative

cat > generate_report.sh <<'EOF'
#!/bin/bash

REPORT="comparative_genomics_report.txt"

cat > $REPORT <<REPORT_EOF
================================================================================
                    COMPARATIVE GENOMICS ANALYSIS REPORT
================================================================================

Date: $(date)
Project: Coral-Associated Bacterial Genomes

================================================================================
1. AVERAGE NUCLEOTIDE IDENTITY (ANI)
================================================================================

$(cat ani/ani_matrix_formatted.tsv)

Species Groups (ANI ≥ 95%):
REPORT_EOF

# Add ANI interpretation
python3 <<PYTHON_EOF >> $REPORT
import pandas as pd
matrix = pd.read_csv('ani/ani_matrix_formatted.tsv', sep='\t', index_col=0)
for query in matrix.index:
    related = [ref for ref in matrix.columns if query != ref and matrix.loc[query, ref] >= 95]
    if related:
        print(f"  {query}: {', '.join(related)}")
PYTHON_EOF

cat >> $REPORT <<REPORT_EOF

================================================================================
2. PANGENOME ANALYSIS
================================================================================

$(cat pangenome/core_genes.tsv | head -1)
Number of core genes: $(wc -l < pangenome/core_genes.tsv | awk '{print $1-1}')

Unique genes per genome:
REPORT_EOF

for file in pangenome/unique_genes_*.tsv; do
    GENOME=$(basename $file .tsv | sed 's/unique_genes_//')
    COUNT=$(wc -l < $file | awk '{print $1-1}')
    echo "  $GENOME: $COUNT" >> $REPORT
done

cat >> $REPORT <<REPORT_EOF

================================================================================
3. INTERPRETATION
================================================================================

Findings:
1. ANI analysis confirms species boundaries
2. Core genome represents essential housekeeping functions
3. Accessory genome reflects niche adaptation
4. Unique genes may include:
   - Horizontal gene transfer events
   - Strain-specific metabolic capabilities
   - Antibiotic resistance genes
   - Virulence factors

Recommendations for further analysis:
- Functional enrichment of accessory genes
- Analysis of unique genes for adaptive traits
- Identification of horizontally transferred genes
- Pathway analysis for metabolic differences

================================================================================
END OF REPORT
================================================================================
REPORT_EOF

cat $REPORT
EOF

chmod +x generate_report.sh
./generate_report.sh
```

### Step 5.2: View Report

```bash
# View the complete report
less comparative_genomics_report.txt

# Save a copy
cp comparative_genomics_report.txt ../../docs/
``` 

---

## Part 6: Advanced Analysis (Optional)

### Step 6.1: Functional Enrichment of Accessory Genes

```bash
# Step 1: Extract gene IDs from accessory orthogroups
echo "Extracting genes from accessory orthogroups..."

# Get the orthogroups file (not gene count, but actual gene assignments)
OG_FILE=$(find . -name "Orthogroups.tsv" | head -1)

if [ -f "$OG_FILE" ]; then
    # Create list of accessory orthogroup IDs
    cut -f1 ../pangenome/accessory_genes.tsv | tail -n +2 > accessory_og_ids.txt
    
    # Extract genes from these orthogroups
    grep -F -f accessory_og_ids.txt "$OG_FILE" > accessory_og_genes.tsv
    
    echo "Accessory orthogroup genes saved to: accessory_og_genes.tsv"
fi

# Step 2: Link to Bakta annotations
echo ""
echo "To perform functional enrichment:"
echo ""
echo "1. Extract gene IDs from accessory_og_genes.tsv"
echo "2. Match gene IDs to Bakta annotation files (*.tsv in annotation directories)"
echo "3. Extract COG/GO/EC categories for each gene"
echo "4. Count category frequencies in accessory vs core genes"
echo ""
echo "Example workflow:"
echo ""
echo "# Extract COG categories from Bakta TSV files"
echo 'for ANNOT_DIR in ../../annotation/*; do'
echo '    SAMPLE=$(basename $ANNOT_DIR)'
echo '    if [ -f "$ANNOT_DIR/${SAMPLE}.tsv" ]; then'
echo '        awk -F"\t" '"'NR>1 && \$11 != "" {print \$1"\t"\$11}'" \\'
echo '            "$ANNOT_DIR/${SAMPLE}.tsv" >> all_cog_assignments.tsv'
echo '    fi'
echo 'done'
echo ""
echo "# Count COG categories in accessory genes"
echo '# Join accessory gene IDs with COG assignments'
echo '# Use Fisher exact test or hypergeometric test for enrichment'
echo ""
echo "Recommended tools for complete analysis:"
echo "  - eggNOG-mapper: Functional annotation and enrichment"
echo "  - InterProScan: Domain and GO term annotation"
echo "  - KEGG pathway analysis: Metabolic pathway reconstruction"
echo "  - R packages: clusterProfiler, GOstats for enrichment testing"
```

### Step 6.2: Identify Horizontal Gene Transfer

```bash
echo "Horizontal Gene Transfer (HGT) Detection:"
echo ""
echo "Methods to identify HGT:"
echo "1. Phylogenetic incongruence (gene tree vs species tree)"
echo "2. Atypical GC content or codon usage"
echo "3. Presence of mobile genetic elements"
echo "4. Genes unique to one strain within a species"
echo ""
echo "Tools for HGT detection:"
echo "- IslandViewer (genomic islands)"
echo "- Alien_hunter (codon usage)"
echo "- HGT-Finder (phylogenetic)"
```-->

---

## Summary and Next Steps

### What We Learned

1. **ANI Analysis**:
   - Calculated genome-wide similarity
   - Confirmed species boundaries (ANI ≥ 95%)
   - Identified closely related strains

2. **Pangenome Analysis**:
   - Identified core genes (essential functions)
   - Found accessory genes (niche adaptation)
   - Discovered unique genes (strain-specific)

3. **Phylogenomics**:
   - Built species tree from core genes
   - Visualized evolutionary relationships
   - Confirmed taxonomic classifications

### Key Takeaways

- **Core genome** = conserved essential functions (e.g., DNA replication, translation)  
- **Accessory genome** = variable functions for adaptation (e.g., metabolism, stress response)  
- **Unique genes** = strain-specific traits (e.g., antibiotic resistance, novel pathways)  
- **ANI ≥ 95%** = same species; **<95%** = different species

### Next Steps

- **Tutorial 5**: Pathway analysis and metabolic reconstruction
- **Tutorial 6?**: Integration and biological interpretation
---

## Troubleshooting

### OrthoFinder Crashes

```bash
# If OrthoFinder runs out of memory:
# - Reduce thread count: -t 4 instead of -t 8
# - Use faster, less memory-intensive mode: -S diamond (default, recommended)
# - Skip MSA for lower memory usage: remove -M msa flag

# For more sensitive (but slower and more memory-intensive) results:
# - Use BLAST: -S blast (note: slower than diamond)
```

### Missing Dependencies

```bash
# Install missing tools
conda install -c bioconda fastani orthofinder
conda install -c conda-forge biopython matplotlib seaborn scipy pandas
```

### File Path Errors

```bash
# Verify all input files exist
find ../../assembly -name "*_polished.fasta"
find ../../annotation -name "*.faa"

# Check file permissions
ls -l proteins/
```

---

## Exercise Questions

1. **ANI Interpretation**:
   - Which genomes belong to the same species?
   - What is the ANI between *Microbacterium* and *Alteromonas*?

2. **Pangenome**:
   - What percentage of genes are core vs. accessory?
   - Which genome has the most unique genes?

3. **Evolution**:
   - Do strains from the same coral fragment cluster together?
   - What does this tell you about acquisition timing?

4. **Function**:
   - What functions might be in the accessory genome?
   - Why would unique genes be interesting for biotechnology?

---

**End of Tutorial 4**
