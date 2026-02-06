# Tutorial 5: Functional Annotation with KofamScan and Pathway Analysis

## Overview
This tutorial guides you through functional annotation of prokaryotic genomes using KofamScan and visualization of KEGG pathways using KEGG-decoder.

## Prerequisites
- Assembled and annotated genome(s)
- Protein FASTA files from gene prediction
- KofamScan installed
- KEGG-decoder installed

### Required Software

```bash
# Activate comparative genomics environment
conda activate VTK2026_pathways

# Verify tools are installed
exec_annotation -h # KOFAMSCAN
KEGG-decoder -h 
```

## Part 1: Running KofamScan

### 1.1 Prepare Input Files
```bash
# Ensure you have protein sequences in FASTA format
# Example: proteins.faa
```

### 1.2 Run KofamScan
```bash
# Basic KofamScan command
exec_annotation -o output.txt \
    --format detail-tsv \
    --ko-list ~/databases/kofamscan/ko_list \
    --profile ~/databases/kofamscan/profiles \
    --cpu 4 \
    proteins.faa

# For multiple genomes, create a loop
for genome in genomes/*.faa; do
        base=$(basename ${genome} .faa)
        exec_annotation -o ${base}_kofam.txt \
            --format detail-tsv \
            --cpu 4 \
            ${genome}
done
```

### 1.3 Output Files
KofamScan produces tab-separated files with KO assignments for each protein.

## Part 2: Pathway Visualization with KEGG-decoder

### 2.1 Single Genome Analysis
```bash
# Generate pathway visualization
KEGG-decoder -i kofam_output.txt -o results_visualization.pdf
```

### 2.2 Multiple Genome Comparison
```bash
# Combine multiple KofamScan outputs
KEGG-decoder -i genome1_kofam.txt,genome2_kofam.txt,genome3_kofam.txt \
    -o comparative_pathways.pdf

# Or use a list file
KEGG-decoder -i sample_list.txt -o comparative_pathways.pdf
```

### 2.3 Customizing Output
```bash
# Generate heatmap for comparison
KEGG-decoder -i *.txt -o pathways_heatmap.pdf --vizoptions static

# Export as interactive HTML
KEGG-decoder -i *.txt -o pathways.html --vizoptions interactive
```

## Part 3: Interpreting Results

### 3.1 Understanding KO Numbers
- KO (KEGG Orthology) numbers represent functional orthologs
- Each KO belongs to specific pathways and modules

### 3.2 Pathway Completeness
- KEGG-decoder shows pathway completeness as color-coded cells
- Green: complete pathway
- Yellow/Orange: partial pathway
- Red/absent: missing pathway

### 3.3 Comparative Analysis
When comparing multiple genomes:
- Identify core metabolic pathways present in all samples
- Detect unique pathways in specific genomes
- Assess metabolic potential differences

## Expected Outputs
- `*_kofam.txt`: KofamScan annotation results
- `comparative_pathways.pdf`: Visual comparison of pathways
- Heatmap showing pathway presence/absence across genomes

## Tips
- Use sufficient CPU cores (`--cpu`) to speed up analysis
- Check KofamScan database version for reproducibility
- Filter low-confidence assignments if needed

## Further Reading
- [KofamScan documentation](https://github.com/takaram/kofam_scan)
- [KEGG-decoder repository](https://github.com/bjtully/BioData/tree/master/KEGGDecoder)
- KEGG pathway database: https://www.genome.jp/kegg/pathway.html
