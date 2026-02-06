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
            --profile ~/databases/kofamscan/profiles \
            --ko-list ~/databases/kofamscan/ko_list \
            --format detail-tsv \
            --cpu 4 \
            ${genome}
done
```

### 1.3 Output Files
KofamScan produces tab-separated files with KO assignments for each protein.

**Example detail-tsv output:**
```
#	gene name	KO	thrshld	score	E-value	"KO definition"
#	---------	------	-------	------	---------	-------------
*	ECENAK_00001	K00145	118.60	125.5	5.1e-37	"N-acetyl-gamma-glutamyl-phosphate reductase [EC:1.2.1.38]"
	ECENAK_00001	K05829	412.40	105.2	7.7e-31	"[amino group carrier protein]-6-phospho-L-2-aminoadipate/5-phospho-L-glutamate reductase [EC:1.2.1.103 1.2.1.106]"
	ECENAK_00001	K12659	668.37	62.7	3.2e-18	"N-acetyl-gamma-glutamyl-phosphate reductase / acetylglutamate kinase [EC:1.2.1.38 2.7.2.8]"
	ECENAK_00001	K00133	189.07	16.7	0.00047	"aspartate-semialdehyde dehydrogenase [EC:1.2.1.11]"
*	ECENAK_00002	K00620	175.40	462.5	4e-139	"glutamate N-acetyltransferase / amino-acid N-acetyltransferase [EC:2.3.1.35 2.3.1.1]"
*	ECENAK_00003	K00930	221.03	439.7	2.6e-132	"acetylglutamate kinase [EC:2.7.2.8]"
	ECENAK_00003	K05828	221.53	170.3	9.9e-51	"[amino group carrier protein]-L-2-aminoadipate/L-glutamate 6-kinase [EC:2.7.2.17 2.7.2.19]"
	ECENAK_00003	K14682	227.37	136.1	2.8e-40	"amino-acid N-acetyltransferase [EC:2.3.1.1]"
```

**Column descriptions:**
- **gene name**: Protein/gene identifier from input FASTA
- **KO**: KEGG Orthology identifier (or `*` if no significant hit)
- **thrshld**: Score threshold for the KO
- **score**: Alignment score
- **E-value**: Expected value (statistical significance)
- **KO definition**: Functional description of the KO 

#### 1.4 mapper vs detail-tsv

KofamScan supports two main output formats:

**detail-tsv format:**
- Provides comprehensive information including scores, thresholds, and E-values
- Best for quality control and filtering results
- Contains all annotation details for each protein
- Default format shown in examples above

**mapper format:**
- Simplified format with only gene names and KO assignments
- Required input format for KEGG-decoder
- One line per KO assignment
- More compact and easier to process downstream

**Example mapper output:**
```
protein_001	K00001
protein_002	K00002
protein_003	K00003
protein_005	K00005
protein_007	K00009
protein_008	K00010
```

Note: Proteins without KO assignments are omitted in mapper format.

## Part 2: Pathway Visualization with KEGG-decoder


<details>
<summary><b>Converting detail-tsv to mapper format</b></summary>

If you have detail-tsv output but need mapper format for KEGG-decoder, you can convert it using standard Unix tools:

```bash
# Convert single file (skip header, filter out non-hits, extract columns 1-2)
grep -v '^#' genome_kofam_detail.txt | \
    grep -w '*' | \
        cut -f2,3 > genome_kofam_mapper.txt

# Batch conversion for multiple files
for file in *_detail.txt; do
    base=$(basename ${file} _detail.txt)
    grep -v '^#' ${file} | \
        grep -w '*' | \
            cut -f2,3 > ${base}_mapper.txt
done
```

**Explanation:**
- `grep -v '^#'`: Remove comment lines starting with #
- `grep -w '*'`: Filter out proteins without KO assignments
- `cut -f2,3` or `print $2"\t"$3`: Extract only gene name and KO columns

</details> 

### 2.1 Single Genome Analysis
```bash
# Generate pathway visualization
KEGG-decoder -i kofam_output.txt -o results_visualization.pdf
```
**Requires mapper formatted output**

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
KEGG-decoder -i *.txt -o pathways_heatmap.pdf -v static

# Export as interactive HTML
KEGG-decoder -i *.txt -o pathways.html -v interactive
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
