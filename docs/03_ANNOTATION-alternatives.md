# Alternative Annotation and Analysis Tools for Bacterial Genomics

## Overview
This guide covers annotation and taxonomic analysis tools for bacterial genomic data. Note that **BASys2**, **TYGS**, and **GGDC** are web-based analysis platforms, not archival databases. For public data deposition, submit to NCBI GenBank, ENA, or DDBJ.

## BaSys 2 (https://basys2.ca/)
**Purpose:** Bacterial genome annotation web service
- **What it does:** Performs automated annotation of bacterial genomes, identifying genes and predicting functions; results are provided to users for download but not permanently housed in a public database
- **Data requirements:**
    - Genomic sequences (FASTA format)
    - Can also accept FASTQ (raw reads) or GenBank format files

## TYGS (https://tygs.dsmz.de/)
**Purpose:** Genome analysis and taxonomic placement platform
- **What it does:** Accepts uploaded genome assemblies to perform whole-genome phylogenetic analysis and species/subspecies boundary estimates; results are analytical outputs, not public repository submissions
- **Data requirements:**
    - Complete or draft genome assemblies (user uploads for analysis)
    - Type strain designation confirmation
    - Strain metadata (isolation source, date, location)
    - 16S rRNA gene sequence (optional)

## GGDC (https://ggdc.dsmz.de/)
**Purpose:** Genome-to-Genome Distance Calculator (analysis service)
- **What it does:** Accepts uploaded sequence files to compute intergenomic distances for species delineation and phylogenetic analysis; does not permanently house data
- **Data requirements:**
    - Genome sequences in FASTA format
    - Genome assembly quality metrics
    - Corresponding metadata

## Preparation for Computational Tools (BaSys2, TYGS, GGDC)
1. Assembled genome sequences in FASTA or GenBank format
2. Optional annotation files (GFF, etc.)
3. Strain metadata if required by the tool
4. **Note:** These tools perform their own validation and taxonomic analysis; pre-running standalone taxonomic validation is not required

## Preparation for Public Database Submission (GenBank/ENA)
1. Validate genome assembly quality
2. Verify taxonomic information
3. Prepare metadata in standardized format
4. Check sequence file formats
5. Submit through respective web interfaces
