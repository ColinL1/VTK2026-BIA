# VTK2026-Bacterial Isolate Assembly

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains a complete bioinformatics workflow for analyzing Oxford Nanopore Technology (ONT) sequencing data from bacterial isolates. Designed for the **[Voolstra lab](https://www.biologie.uni-konstanz.de/voolstra/) 2026 VTK**, this pipeline guides through quality control, genome assembly, functional annotation, comparative genomics, and metabolic pathway analysis.

### Objectives

1. **Quality Control**: Assess and filter ONT sequencing data
2. **Assembly**: De novo genome assembly without reference genomes
3. **Annotation**: Taxon-independent functional annotation
4. **Gene Screening**: Identify secondary metabolites, AMR, and coral-specific traits
5. **Comparative Genomics**: Pangenome analysis and phylogenomics
6. **Pathway Analysis**: Reconstruct metabolic pathways and coral symbiosis functions

### Key Features

- **Complete workflow**: From raw reads to publication-ready analyses
- **Gene-focused**: Screens for relevant genes and pathways
- **Reproducible**: Conda environment with all dependencies
<!-- - **Multi-level comparisons**: Strain, species, and genus-level analyses #TODO: add back in it? --> 

---

## Detailed Documentation

- **[Setup Guide](docs/SETUP.md)**: Installation, database downloads, and troubleshooting
- **[Guide](docs/README.md)**: Quality control, assembly, and annotation, etc.

---

## Repository Structure

```
VTK2026/
├── README.md                       # This file
├── env_files                       # Conda environment specification
├── data/
│   ├── raw_reads/                  # ONT sequencing data (.fastq.gz)
│   ├── metadata/               
│   │   └── sample_info.csv         # Sample metadata (if needed)
│   └── custom_gene_databases/      # Reference sequences for screening
├── scripts/                        # Analysis scripts  (released at the end of each module) 
│   ├── 01_quality_control.sh       # QC and filtering
│   ├── 02_assembly.sh              # Flye + Medaka assembly
│   ├── 03_annotation.sh            # Bakta annotation
|   ├── 04_gene_screening.sh        # antiSMASH, ABRicate, custom genes
|   ├── 05_comparative_genomics.sh  # OrthoFinder, ANI, pangenome
|   ├── 06_pathway_analysis.sh      # KEGG pathways, visualization
|   └── helpers/
├── results/                        # All output files (generated)
│   ├── qc/                         # Quality control reports
│   ├── assembly/                   # Genome assemblies
│   ├── screening/                  # Trait and gene screening
│   ├── annotation/                 # Gene annotations
│   ├── comparative/                # Comparative genomics
│   └── pathways/                   # Pathway analysis
└── docs/
    ├── INTRODUCTION.md
    ├── README.md
    ├── SETUP.md
    ├── TUTORIAL_01_QC.md
    ├── TUTORIAL_02_ASSEMBLY.md
    ├── TUTORIAL_03_ANNOTATION.md
    └── ETC...
```

---

## Pipeline Workflow

### Day 1-2: Quality Control and Preprocessing

<!-- **Script**: `scripts/01_quality_control.sh` -->

**Tools**: NanoPlot, NanoStat, Porechop_ABI, Filtlong

**Objectives**:

- Assess raw ONT read quality and yield
- Trim adapter sequences
- Filter to 100x coverage with highest quality reads
- Remove short reads (<1 kb)
<!-- 
```bash
bash scripts/01_quality_control.sh <sample_id> <threads>
``` -->

### Day 3-4: Genome Assembly

<!-- **Script**: `scripts/02_assembly.sh` -->

**Tools**: Flye, Medaka, QUAST, BUSCO

**Objectives**:

- De novo assembly with Flye
- Consensus polishing with Medaka
- Quality assessment (N50, L50, contig count)
- Completeness evaluation with BUSCO

**Expected**: 3-5 MB genomes, 1-10 contigs, >90% BUSCO completeness
<!-- 
```bash
bash scripts/02_assembly.sh <sample_id> <threads> <genome_size>
``` -->

### Day 5: Taxon-Independent Functional Annotation

<!-- **Script**: `scripts/03_annotation.sh` -->

**Tools**: Bakta

**Objectives**:

- Predict genes and functional annotations
- Assign EC numbers, GO terms, COG categories
- Minimize "hypothetical proteins" using comprehensive database
- Generate GenBank, GFF3, and FASTA outputs
<!-- 
```bash
bash scripts/03_annotation.sh <sample_id> <threads>
``` -->

---

### Week 2: Comparative and Functional Genomics

### Day 6-7: Gene and Trait Screening

<!-- **Script**: `scripts/04_gene_screening.sh` -->

**Tools**: antiSMASH, ABRicate, AMRFinderPlus, BLAST

**Objectives**:

- Detect secondary metabolite biosynthetic gene clusters
- Screen for antimicrobial resistance genes
- Identify nitrogen cycling genes (nif, nir, nor, nos)
- Find DMSP degradation pathways (dmdA, ddd genes)
- Detect Type VI secretion systems
- Assess vitamin biosynthesis capabilities
<!-- 
```bash
bash scripts/04_gene_screening.sh <sample_id> <threads>
``` -->

### Day 8-9: Strain-Level and Species-Level Comparisons

<!-- **Script**: `scripts/05_comparative_genomics.sh` -->

**Tools**: FastANI, OrthoFinder

**Objectives**:
- Calculate Average Nucleotide Identity (ANI >95% = same species)
- Cluster orthologous genes across genomes
- Define core vs accessory genomes
- Identify source-specific genes (F003 vs H2)
- Build phylogenomic trees

**Multi-level analysis**:
- Species-level: Within *M. enclense*, *M. ginsengisoli*, *A. portus*
- Genus-level: All *Microbacterium*, all *Alteromonas*
- Cross-genus: All 12 isolates
<!-- 
```bash
# Species-level comparison
bash scripts/05_comparative_genomics.sh species <threads>

# Genus-level comparison
bash scripts/05_comparative_genomics.sh genus <threads>

# All isolates
bash scripts/05_comparative_genomics.sh all <threads>
``` -->

### Day 10: Pathway Analysis and Coral-Context Integration

<!-- **Script**: `scripts/06_pathway_analysis.sh` -->

**Tools**: KEGG Mapper, custom Python scripts

**Objectives**:
- Extract EC numbers, GO terms, COG categories
- Map to KEGG metabolic pathways
- Assess pathway completeness for coral-relevant functions:
  - Nitrogen metabolism
  - Sulfur metabolism (DMSP)
  - Vitamin biosynthesis (B12, B1, B2)
  - Oxidative stress response
- Compare metabolic profiles across samples
- Visualize pathway distributions
<!-- 
```bash
# Individual sample
bash scripts/06_pathway_analysis.sh <sample_id> <threads>

# Comparative analysis
bash scripts/06_pathway_analysis.sh comparative <threads>
``` -->

---

## Expected Outcomes

By completing this workshop, students will generate:

1. **high-quality genome assemblies** with comprehensive quality metrics
2. **Functional annotations** with COG, GO, and EC classifications
3. **Trait screening results** for coral-relevant genes and pathways
4. **Pangenome analyses** at multiple taxonomic levels
5. **Metabolic pathway reconstructions** focused on coral symbiosis
6. **Publication-ready figures** and summary tables

---

## Installation

See [docs/SETUP.md](docs/SETUP.md) for detailed installation instructions.

---

## License

MIT License - see LICENSE file for details

---

## Support

For issues or questions:
1. Check [docs/SETUP.md](docs/SETUP.md) troubleshooting section
2. Review [docs/README.md](docs/README.md) for guidance
3. Check individual tool documentation

---

## Acknowledgments

This workshop was designed for teaching bacterial genomics with a focus on coral-associated microbiomes. Special thanks to the developers of all open-source bioinformatics tools used in this pipeline.
