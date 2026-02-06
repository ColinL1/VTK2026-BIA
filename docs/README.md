# Quality Control Through De Novo Annotation

## Overview
This guide focuses on transforming raw ONT sequencing data into annotated bacterial genomes. You'll learn quality control, de novo assembly, and functional annotation without relying on reference genomes.

**Learning Objectives**:
- Assess and improve ONT read quality
- Perform de novo genome assembly
- Evaluate assembly quality metrics
- Annotate genes and predict functions
- Perform pathways investigation
- Compare genomes (ANI and orthologous)
- Explore BGCs and AMR

---

## Sturcture

1. [QC](01_QC.md)
2. [Assembly](02_ASSEMBLY.md)
3. [Annotation](03_ANNOTATION.md)
   3.1 [Online tools](03_ANNOTATION-alternatives.md)
4. [Comparisons](04_COMPARISON.md)
5. [Pathways](05_PATHWAYS.md)
6. [Gene screening](06_GENE-SCREENING.md)

---

## Detailed Workflow

### Tutorial 1: Quality Control ([docs/01_QC.md](docs/01_QC.md))

**Tools:** NanoPlot, NanoStat, Porechop_ABI, Filtlong

**Learning Goals:**
- Assess raw ONT read quality (N50, read length distribution, quality scores)
- Remove adapter sequences with Porechop_ABI
- Filter reads to target coverage (typically 100-150x)
- Remove short reads (<1 kb)
- Compare before/after filtering statistics

**Key Outputs:**
- Quality reports (NanoPlot HTML)
- Filtered reads (`*_filtered.fastq.gz`)
- Summary statistics

**Expected Time:** 1-2 hours per sample

---

### Tutorial 2: Genome Assembly ([docs/02_ASSEMBLY.md](docs/02_ASSEMBLY.md))

**Tools:** Flye, Medaka, QUAST, BUSCO

**Learning Goals:**
- Perform de novo assembly with Flye
- Polish assemblies with Medaka (error correction)
- Assess quality with QUAST (N50, L50, total length)
- Evaluate completeness with BUSCO (>90% expected)
- Visualize assembly graphs with Bandage

**Key Outputs:**
- Polished assembly (`*_polished.fasta`)
- QUAST report (assembly statistics)
- BUSCO completeness assessment
- Assembly graph visualization

**Expected Results:** 3-5 MB genomes, 1-10 contigs, >90% BUSCO completeness

**Expected Time:** 2-4 hours per sample

---

### Tutorial 3: Functional Annotation ([docs/03_ANNOTATION.md](docs/03_ANNOTATION.md))

**Tools:** Bakta

**Learning Goals:**
- Perform comprehensive gene prediction
- Assign functional annotations (COG, GO, EC numbers)
- Understand annotation file formats (GFF3, GenBank, FASTA)
- Analyze gene content and functional categories
- Extract annotation statistics

**Key Outputs:**
- Annotated GenBank file (`*.gbff`)
- Gene feature file (`*.gff3`)
- Protein sequences (`*.faa`)
- Nucleotide sequences (`*.ffn`)
- Annotation summary statistics

**Expected Results:** 3000-5000 genes, comprehensive functional assignments

**Expected Time:** 30-60 minutes per sample

---

### Tutorial 4: Comparative Genomics ([docs/04_COMPARISON.md](docs/04_COMPARISON.md))

**Tools:** FastANI, OrthoFinder

**Learning Goals:**
- Calculate Average Nucleotide Identity (ANI ≥95% = same species)
- Define species boundaries across your isolates
- Cluster orthologous genes with OrthoFinder
- Identify core genome (genes in all isolates)
- Identify accessory genome (genes in some isolates)
- Identify unique genes (genes in single isolates)
- Construct phylogenomic trees

**Key Outputs:**
- ANI matrix (`ani_matrix.tsv`)
- OrthoFinder results (orthogroups, gene trees)
- Core/accessory genome lists
- Species phylogenetic tree

**Expected Results:** Clear species groups, defined pangenome structure

**Expected Time:** 1-3 hours (depending on sample count)

---

### Tutorial 5: Pathway Analysis ([docs/05_PATHWAYS.md](docs/05_PATHWAYS.md))

**Tools:** KofamScan, KEGG-decoder

**Learning Goals:**
- Assign KEGG orthology (KO) identifiers
- Map genes to metabolic pathways
- Assess pathway completeness
- Visualize KEGG pathways
- Compare metabolic capabilities across samples

**Key Outputs:**
- KO assignments
- KEGG pathway maps
- Pathway completeness matrices
- Metabolic profile comparisons

**Expected Results:** Identification of key metabolic capabilities

**Expected Time:** 1-2 hours per sample

---

### Tutorial 6: Gene Screening ([docs/06_GENE-SCREENING.md](docs/06_GENE-SCREENING.md))

**Tools:** antiSMASH, ABRicate, AMRFinderPlus, BLAST

**Learning Goals:**
- Detect biosynthetic gene clusters (BGCs) with antiSMASH
- Screen for antimicrobial resistance (AMR) genes
- Identify nitrogen cycling genes (nif, nir, nor, nos)
- Find DMSP degradation pathways
- Detect Type VI secretion systems
- Assess vitamin biosynthesis capabilities
- Interpret results in ecological context

**Key Outputs:**
- antiSMASH BGC predictions
- ABRicate AMR gene hits
- AMRFinderPlus comprehensive AMR report
- Custom gene screening results

**Expected Results:** Identification of ecological traits and biosynthetic potential

**Expected Time:** 1-2 hours per sample (antiSMASH is time-intensive)

---
General file formats information: [INTRODUCTION](INTRODUCTION.md) | Software used in this tutorial: [SOFTWARE](SETUP.md)

---
