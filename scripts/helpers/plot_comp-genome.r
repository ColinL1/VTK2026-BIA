#!/usr/bin/env Rscript
# Create pangenome visualization using R

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

# Find and read orthogroups file
og_file <- "results/comparative/orthofinder/OrthoFinder/Results_Feb04/Orthogroups/Orthogroups.GeneCount.tsv"
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