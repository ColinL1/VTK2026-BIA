# Scripts Directory

This directory contains analysis scripts that cover all the steps described in the 2026 VTK course. See [docs](docs/README.md)

---

## Analysis Scripts

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

### `04_gene_screening.sh`
Screen for biosynthetic gene clusters, AMR genes, and other features.

**Usage**:
```bash
conda activate VTK2026_gene-screening
bash 04_gene_screening.sh <assembly> <gff> <proteins> <threads>
```

**Example**:
```bash
bash 04_gene_screening.sh results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta \
  results/annotation/F003_M_enclense/F003_M_enclense.gff3 \
  results/annotation/F003_M_enclense/F003_M_enclense.faa 8
```

**Tools used**: antiSMASH, ABRicate, AMRFinderPlus

**Outputs**:
- `results/screening/<sample>/` - Screening results

---

> **Note**: Scripts `05_comparative_genomics.sh` and `06_pathway_analysis.sh` are not yet final.


### Logging

Redirect output to log files for later review:

```bash
bash 02_assembly.sh F003_M_enclense 8 3.5m 2>&1 | tee assembly_F003.log
```

---

## Best Practices

1. **Test on one sample first** before running batch processing
2. **Monitor resource usage** with `htop` or `top`
3. **Keep intermediate files** until confirming final results are correct
4. **Document parameter changes** in your lab notebook
5. **Back up results** regularly to prevent data loss
6. **Use version control** for any custom modifications

---

For detailed explanations and exercises, see: [Step-by-step Guide](../docs/README.md) | Software used in this tutorial: [SOFTWARE](SETUP.md) |  [Setup Guide](../docs/SETUP.md)
