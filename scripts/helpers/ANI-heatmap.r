#!/usr/bin/env Rscript
# Create ANI heatmap and dendrogram using R

# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Read ANI matrix
ani_matrix <- read.table('results/comparative/ani/ani_matrix.tsv', header=TRUE, row.names=1, sep='\t')
data <- read.table('results/comparative/ani/ani_matrix.tsv', header=FALSE, sep='\t')
colnames(data) <- c('Query', 'Reference', 'ANI', 'Matches', 'Fragments')

# Convert to matrix format
ani_matrix <- xtabs(data[[3]] ~ data[[1]] + data[[2]])
# Replace NA with 0 for visualization purposes
ani_matrix[is.na(ani_matrix)] <- 0

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