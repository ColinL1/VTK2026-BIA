# VTK2026 - Bacterial Isolate Assembly

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains a complete bioinformatics workflow for analyzing **Oxford Nanopore Technology (ONT)** sequencing data from bacterial isolates. Designed for the **[Voolstra lab](https://www.biologie.uni-konstanz.de/voolstra/) 2026 VTK**, this tutorial guides you through quality control, de novo genome assembly, functional annotation, comparative genomics, pathway analysis, and gene screening.

**The course focuses on going form raw ONT reads to annotated bacterial genomes without requiring reference sequences.**

---

> **Note:** For a fully automated workflow, check out the [Nextflow version](https://github.com/ColinL1/ONT-Bacterial-Assembly/) of this pipeline, which streamlines batch processing and reproducibility across multiple samples.

---

## Learning Objectives

By completing this workshop, you will be able to:

- Assess and improve ONT read quality
- Perform de novo genome assembly
- Evaluate assembly quality metrics
- Annotate genes and predict functions
- Identify biosynthetic gene clusters and AMR genes
- Compare genomes using ANI and pangenome analysis
- Reconstruct and visualize metabolic pathways
- Interpret results in an ecological context

---

## Course Structure

This workshop is organized into six main tutorials:

1. **[Quality Control](docs/01_QC.md)** - Assess, trim, and filter ONT reads
2. **[Genome Assembly](docs/02_ASSEMBLY.md)** - De novo assembly with Flye and Medaka
3. **[Functional Annotation](docs/03_ANNOTATION.md)** - Gene prediction with Bakta
4. **[Comparative Genomics](docs/04_COMPARISON.md)** - ANI analysis and pangenome characterization
5. **[Pathway Analysis](docs/05_PATHWAYS.md)** - KEGG pathway reconstruction
6. **[Gene Screening](docs/06_GENE-SCREENING.md)** - Detect BGCs, AMR, and trait-specific genes

### Additional Resources

- **[Introduction to File Formats](docs/INTRODUCTION.md)** - FASTA, FASTQ, GFF, and quality scores
- **[Setup Guide](docs/SETUP.md)** - Installation, database downloads, and environment configuration
- **[CLI Basics](https://github.com/ColinL1/CLI-exercises)** - Command-line fundamentals

---

## Quick Start

See [Setup Guide](docs/SETUP.md) for complete installation instructions.

### 2. Run the Pipeline

**Manual approach (recommended for learning):**
Follow the step-by-step tutorials in [docs/](docs/README.md)

**Automated approach (for batch processing):**
Use the provided scripts:

```bash
# Quality control
bash scripts/01_quality_control.sh <sample_name.fastq.gz> 8

# Assembly
bash scripts/02_assembly.sh <sample_name> 8

# Annotation
bash scripts/03_annotation.sh <sample_name> 8

# Gene screening
bash scripts/04_gene_screening.sh <sample_name> 8
```
<!-- 
# Comparative genomics (requires multiple samples)
bash scripts/05_comparative_genomics.sh all 8

# Pathway analysis
bash scripts/06_pathway_analysis.sh <sample_name> 8 -->

---

## Expected Outcomes

By completing this workshop, you will have generated:

1. **High-quality genome assemblies** with comprehensive quality metrics (N50, BUSCO)
2. **Functional annotations** with COG, GO, and EC classifications for all genes
3. **Gene screening results** identifying BGCs, AMR genes, and ecological traits
4. **Comparative genomics analyses** including:
   - ANI matrices showing species relationships
   - Pangenome characterization (core/accessory/unique genes)
   - Phylogenomic trees
5. **Metabolic pathway reconstructions** using KEGG annotations
6. **Publication-ready visualizations** and summary tables

### Skills Gained

- Command-line proficiency for bioinformatics analyses
- Understanding of ONT sequencing data characteristics
- Experience with industry-standard genomics tools
- Ability to interpret genomic data in ecological context
- Critical evaluation of assembly and annotation quality

---

## Key Tools and Technologies

| Tool | Purpose | Tutorial |
|------|---------|----------|
| NanoPlot | Read quality visualization | 1 |
| Porechop_ABI | Adapter trimming | 1 |
| Filtlong | Quality-based read filtering | 1 |
| Flye | De novo genome assembly | 2 |
| Medaka | Assembly polishing | 2 |
| QUAST | Assembly quality assessment | 2 |
| BUSCO | Completeness evaluation | 2 |
| Bakta | Comprehensive genome annotation | 3 |
| FastANI | Average Nucleotide Identity | 4 |
| OrthoFinder | Ortholog clustering, pangenome | 4 |
| KofamScan | KEGG orthology assignment | 5 |
| KEGG-decoder | Pathway visualization | 5 |
| antiSMASH | Biosynthetic gene cluster detection | 6 |
| ABRicate | AMR gene screening | 6 |
| AMRFinderPlus | Comprehensive AMR analysis | 6 |

---

## Prerequisites

- **Basic command-line skills** (recommended: [CLI exercises](https://github.com/ColinL1/CLI-exercises))
- **Conda/Mamba** installed ([installation guide](docs/SETUP.md))
- **~50-100 GB free disk space** for databases and results
- **ONT sequencing data** (FASTQ format)

---

## Citation and License

This workshop was designed for teaching bacterial genomics with a focus on coral-associated microbiomes at the University of Konstanz.

**License:** MIT License - see [LICENSE](LICENSE) file for details

### Tool Citations

Please cite the tools you use in your analyses.

---

**Ready to start?** Begin with the [Setup Guide](docs/SETUP.md) or dive into [Tutorial 1: Quality Control](docs/01_QC.md)!
