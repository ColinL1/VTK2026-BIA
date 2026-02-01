# Setup Guide - Bacterial Genomics Workshop

## Prerequisites

### Required Software
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- Git (for cloning repository)

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd VTK2026
```

Or if you already have the files:
```bash
cd VTK2026
```

### 2. Create Conda Environment

This will install all required bioinformatics tools:

```bash
# Create environment from file
conda env create -f environment.yml

# Activate the environment
conda activate VTK2026
```

Environment files are provided in [env_files](env_files/)

**Note**: Initial installation may take 30-60 minutes as it downloads and installs all dependencies.

If conda is slow, consider using [mamba](https://mamba.readthedocs.io/) instead:
```bash
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

### 3. Download Required Databases

#### Bakta Database (Required for Annotation)

Bakta requires a reference database. Download the full database:

```bash
# Create a directory for databases
mkdir -p ~/databases
cd ~/databases

# Download Bakta database (approximately 30 GB)
bakta_db download --output bakta_db --type full

# Set environment variable (add this to your ~/.bashrc or ~/.zshrc)
export BAKTA_DB=~/databases/bakta_db
```

#### BUSCO Database (Auto-downloaded)

BUSCO will automatically download required lineage datasets on first use:
- `bacteria_odb10` (bacterial marker genes)

#### ABRicate Databases (Auto-updated)

ABRicate databases are updated automatically by the scripts, but you can update them manually:

```bash
# Update all ABRicate databases
abricate --setupdb
```

### 4. Prepare Your Data

#### Directory Structure

Place your ONT sequencing data in the appropriate directory:

```bash
VTK2026/
├── data/
│   ├── raw_reads/          # Place your .fastq or .fastq.gz files here
│   └── metadata/
│       └── sample_info.csv # Update with your sample information
```

#### Sample Naming Convention

Name your raw read files:

```
F003_M_enclense.fastq.gz
F003_M_ginsengisoli.fastq.gz
F003_A_portus.fastq.gz
H2_M_enclense.fastq.gz
H2_M_ginsengisoli.fastq.gz
H2_A_portus.fastq.gz
```

### 5. Test Installation

Verify that all tools are installed correctly:

```bash
# Activate environment
conda activate ENVIROMENT_NAME

# Test key tools
NanoPlot --version
flye --version
bakta --version
orthofinder -h
antismash --version

# Check Python packages
python -c "import pandas, matplotlib, seaborn, biopython; print('Python packages OK')"
```

If any command fails, the tool may not be installed correctly. Try:
```bash
conda activate VTK2026_QC
conda list | grep <tool-name>
```
<!-- ## Quick Start

### Run a Single Sample

Test the pipeline on one sample first:

```bash
# Activate environment
conda activate VTK2026

# Run QC on one sample
bash scripts/week1/01_quality_control.sh F003_M_enclense 8

# If successful, run assembly
bash scripts/week1/02_assembly.sh F003_M_enclense 8 3.5m

# Then annotation
bash scripts/week1/03_annotation.sh F003_M_enclense 8
```

### Run All Samples (Full Pipeline)

Once you've tested with a single sample, run the complete workflow:

```bash
# Make master script executable
chmod +x scripts/run_all.sh

# Run full pipeline (will process all samples sequentially)
bash scripts/run_all.sh 8  # 8 = number of threads
```

**Time estimate**: 2-3 days computational time for 12 samples -->

## Troubleshooting

### Common Issues

#### 1. Conda environment creation fails

**Solution**: Update conda and try again
```bash
conda update -n base conda
conda env create -f environment.yml
```

#### 2. Bakta database not found

**Error**: `Error: BAKTA_DB environment variable not set`

**Solution**: Set the environment variable
```bash
export BAKTA_DB=~/databases/bakta_db
# Add to ~/.bashrc to make permanent
echo 'export BAKTA_DB=~/databases/bakta_db' >> ~/.bashrc
```

#### 3. Out of memory errors

**Solution**: Reduce number of threads or process samples one at a time

#### 4. Disk space issues

**Solution**: Check available space and clean up intermediate files
```bash
df -h  # Check disk space
# Remove large intermediate files after confirming final results are OK
rm -rf results/qc/trimmed/  # Trimmed reads (keep filtered reads)
```

#### 5. Tools not found in PATH

**Solution**: Ensure conda environment is activated
```bash
conda activate VTK2026
which flye  # Should show path in conda environment
```

### Getting Help

1. Check individual tool documentation
2. Review the week-specific guides in `docs/`
3. Examine log files in `results/` directories
4. Check that input files exist and are not corrupted:
   ```bash
   ls -lh data/raw_reads/
   zcat data/raw_reads/<file>.fastq.gz | head
   ```

## Next Steps

Once setup is complete, proceed to:
- [Documentation](docs/README.md) - Quality Control, Assembly, and Annotation, etc.

### Additional Reading
- [Flye documentation](https://github.com/fenderglass/Flye)
- [Bakta documentation](https://github.com/oschwengers/bakta)
- [OrthoFinder documentation](https://github.com/davidemms/OrthoFinder)
- [antiSMASH documentation](https://docs.antismash.secondarymetabolites.org/)
- 
