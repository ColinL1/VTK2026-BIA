# Scripts Directory

This directory contains analysis scripts that cover all the steps described in the 2026 VTK course. See [docs](docs/README.md)

**What it does**:
- Part 1. QC, Assembly, Annotation
- Part 2. Gene Screening, Comparative Genomics, Pathway Analysis
- Part 3. Visualization

---

## Week 1 Scripts

### `01_quality_control.sh`
Assess raw read quality, trim adapters, and filter reads.

**Usage**:
```bash
bash 01_quality_control.sh <sample_path> <threads>
```

**Example**:
```bash
conda activate VTK2026_QC
bash 01_quality_control.sh F003_M_enclense.fastq.gz 8
```

**Tools used**: NanoPlot, NanoStat, Porechop_ABI, Filtlong

**Outputs**:
- `results/qc/nanoplot_raw/` - Raw read QC reports
- `results/qc/nanoplot_filtered/` - Filtered read QC reports
- `results/qc/filtered/` - Filtered FASTQ files
- `results/qc/*_stats.txt` - Statistics summaries

---

### `02_assembly.sh`
De novo genome assembly and quality assessment.

**Usage**:
```bash
bash 02_assembly.sh <sample_path> <threads> <genome_size>
```

**Example**:
```bash
conda activate VTK2026_Assembly
bash 02_assembly.sh F003_M_enclense.fastq.gz 8 3.5m
```

**Genome sizes**:
- *Microbacterium*: 3-4m (3-4 Mb)
- *Alteromonas*: 4-5m (4-5 Mb)

**Tools used**: Flye, Medaka, QUAST, BUSCO

**Outputs**:
- `results/assembly/<sample>/<sample>_polished.fasta` - Final assembly
- `results/assembly/<sample>/qc/quast/` - QUAST reports
- `results/assembly/<sample>/qc/*_busco/` - BUSCO results
- `results/assembly/<sample>/flye/assembly_graph.gfa` - Assembly graph

---

### `03_annotation.sh`
Functional genome annotation with Bakta.

**Usage**:
```bash
bash 03_annotation.sh <sample_path> <threads>
```

**Example**:
```bash
conda activate VTK2026_Annotate
bash 03_annotation.sh F003_M_enclense 8
```

**Tools used**: Bakta

**Outputs**:
- `results/annotation/<sample>/<sample>.gbff` - GenBank format
- `results/annotation/<sample>/<sample>.gff3` - Gene coordinates
- `results/annotation/<sample>/<sample>.faa` - Protein sequences
- `results/annotation/<sample>/<sample>.tsv` - Annotation table
- `results/annotation/<sample>/<sample>_annotation_stats.txt` - Summary

---
<!-- 
## Part 1 Scripts

### `04_gene_screening.sh`
Screen for biosynthetic gene clusters, AMR genes, and coral-specific traits.

**Usage**:
```bash
bash 04_gene_screening.sh <sample_path> <threads>
```

**Example**:
```bash
bash 04_gene_screening.sh F003_M_enclense 8
```

**Tools used**: antiSMASH, ABRicate, AMRFinderPlus, BLAST

**Outputs**:
- `results/screening/<sample>/antismash/` - BGC predictions
- `results/screening/<sample>/amr/` - AMR gene hits
- `results/screening/<sample>/custom_genes/` - Specific genes
- `results/screening/<sample>/<sample>_screening_summary.txt` - Summary

---

### `05_comparative_genomics.sh`
Multi-level comparative genomics and pangenome analysis.

**Usage**:
```bash
bash 05_comparative_genomics.sh <analysis_type> <threads>
```

**Analysis types**:
- `species` - Compare strains within same species
- `genus` - Compare species within genus
- `all` - Compare all 12 isolates

**Examples**:
```bash
bash 05_comparative_genomics.sh species 8
bash 05_comparative_genomics.sh genus 8
bash 05_comparative_genomics.sh all 8
```

**Tools used**: FastANI, OrthoFinder

**Outputs**:
- `results/comparative/ani/` - ANI matrices
- `results/comparative/orthofinder/` - Ortholog clustering
- `results/comparative/pangenome/` - Core/accessory genome analysis
- `results/comparative/phylogenomics/` - Phylogenetic trees
- `results/comparative/comparative_summary_*.txt` - Summaries

---

### `06_pathway_analysis.sh`
Metabolic pathway reconstruction and visualization.

**Usage**:
```bash
bash 06_pathway_analysis.sh <sample_path> <threads>
bash 06_pathway_analysis.sh comparative <threads>
```

**Examples**:
```bash
# Individual sample
bash 06_pathway_analysis.sh F003_M_enclense 8

# Comparative analysis
bash 06_pathway_analysis.sh comparative 8
```

**Tools used**: Custom Python scripts, KEGG Mapper

**Outputs**:
- `results/pathways/<sample>/ec_numbers.txt` - EC numbers
- `results/pathways/<sample>/kegg_pathways.txt` - KEGG pathway info
- `results/pathways/<sample>/coral_pathways_summary.txt` - Coral-relevant pathways
- `results/pathways/comparative/` - Comparative pathway analysis -->

---

## Tips for Using Scripts

### Running Scripts in Parallel

For faster processing, you can run multiple samples in parallel using GNU Parallel or background processes:

```bash
# Using background processes (careful with resource usage!)
bash 01_quality_control.sh F003_M_enclense 4 &
bash 01_quality_control.sh F003_M_ginsengisoli 4 &
wait  # Wait for all background jobs to complete
```

### Resuming After Interruption

All scripts check for existing output files and can resume from where they left off. If a script fails:

1. Check the error message
2. Fix the issue (disk space, memory, input files)
3. Re-run the same command


### Logging

Redirect output to log files for later review:

```bash
bash 02_assembly.sh F003_M_enclense 8 3.5m 2>&1 | tee assembly_F003.log
```

---

## Troubleshooting

### Common Issues

1. **"Command not found"**
   - Ensure conda environment is activated: `conda activate coral_genomics`

2. **"Permission denied"**
   - Scripts not executable: `chmod +x scripts/*.sh scripts/*/*.sh`

3. **"No such file or directory"**
   - Check input files exist in correct locations
   - Verify paths in `data/metadata/sample_info.csv`

4. **Out of memory**
   - Reduce thread count
   - Process samples one at a time
   - Increase system swap space

5. **Disk space full**
   - Clean up intermediate files
   - Move results to larger drive
   - Remove old results before rerunning


---

## Best Practices

1. **Test on one sample first** before running batch processing
2. **Monitor resource usage** with `htop` or `top`
3. **Keep intermediate files** until confirming final results are correct
4. **Document parameter changes** in your lab notebook
5. **Back up results** regularly to prevent data loss
6. **Use version control** for any custom modifications

---

For detailed explanations and exercises, see:
- [Step-by-step Guide](../docs/README.md)
- [Setup Guide](../docs/SETUP.md)
