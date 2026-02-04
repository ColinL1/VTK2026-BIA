# Tutorial 1: Manual Quality Control and Preprocessing

## Learning Objectives

By the end of this tutorial, you will be able to:
- Assess raw ONT sequencing data quality using NanoPlot
- Understand quality metrics (read length, N50, quality scores)
- Remove adapter sequences using Porechop_ABI
- Filter reads by length and quality using Filtlong
- Compare before/after filtering statistics

---

## Prerequisites

### Required Software

```bash
# Activate the QC conda environment
conda activate VTK2026_QC

# Verify tools are installed
NanoPlot --version
porechop_abi --version
filtlong --version
NanoStat --version
```

### Directory Setup

```bash
# Navigate to your project directory
cd ~/'MY-NAME-EXAMPLE'/ # change this to your own name. and work in here from now on. 

# Create necessary directories
mkdir -p results/qc/nanoplot_raw
mkdir -p results/qc/nanoplot_filtered
mkdir -p results/qc/trimmed
mkdir -p results/qc/filtered
```

---

## Step 1: Examine Your Raw Data

### 1.1 Locate Your FASTQ File

First, find your raw sequencing data:

```bash
# Navigate to raw data directory
cd data/raw_reads

# List available samples
ls -lh

# Example output:
# F003_M_enclense.fastq.gz
# F003_M_ginsengisoli.fastq.gz
# F003_A_portus.fastq.gz
# ... etc
```

### 1.2 Quick Inspection

Before running full QC, let's peek at the file structure:

```bash
# View first read (4 lines in FASTQ)
zcat F003_M_enclense.fastq.gz | head -n 4

# Example output:
# @read_001 runid=abc123 ch=45 start_time=2026-01-20T10:30:00Z
# ATGCGATCGATCGATCGATCGATCGATCGATCGATCG...
# +
# IIIIHHHHGGGGFFFFEEEDDDDCCCCBBBBAAAA...
```

**Understanding the FASTQ format:**
- **Line 1**: `@` + read identifier with metadata
- **Line 2**: DNA sequence (ATGC)
- **Line 3**: `+` separator
- **Line 4**: Quality scores (one character per base)

### 1.3 Count Total Reads

```bash
# Count reads (FASTQ has 4 lines per read)
zcat F003_M_enclense.fastq.gz | wc -l | awk '{print $1/4 " reads"}'

# Alternative: count header lines starting with @
# (May overcount if @ appears in quality scores)
zcat F003_M_enclense.fastq.gz | grep -c "^@"
```

### 1.4 Check File Size

```bash
# Check compressed file size
ls -lh F003_M_enclense.fastq.gz

# Check uncompressed size (without extracting)
zcat F003_M_enclense.fastq.gz | wc -c | awk '{print $1/1024/1024 " MB"}'
```

**💡 Tip:** Never decompress large FASTQ files unless absolutely necessary! Always use `zcat` to work with compressed files directly.

---

## Step 2: Raw Data Quality Assessment with NanoPlot

### 2.1 Understanding NanoPlot

**NanoPlot** is a visualization tool specifically designed for ONT data. It generates:
- Read length distribution histograms
- Quality score distributions
- Read length vs. quality scatter plots
- Summary statistics (N50, mean length, total bases)

### 2.2 Run NanoPlot on Raw Data

```bash
# Navigate your work directory
cd ~/'MY-NAME-EXAMPLE'/

# Run NanoPlot
NanoPlot -t 3 \
    --fastq data/raw_reads/examples.fastq.gz \
    --loglength \
    --plots dot \
    -o results/qc/nanoplot_raw/F003_M_enclense \
    --N50 \
    --title "F003_M_enclense - Raw Reads"
```

**Parameter explanation:**
- `-t 4`: Use 8 CPU threads (adjust based on your system)
- `--fastq`: Input file path
- `--loglength`: Use logarithmic scale for read length (better visualization)
- `--plots dot`: Generate dot plots (alternatives: kde, hex)
- `-o`: Output directory
- `--N50`: Calculate and report N50 statistic
- `--title`: Plot title

### 2.3 Understanding the Output

NanoPlot creates several files:

```bash
# List output files
ls results/qc/nanoplot_raw/F003_M_enclense/

# Key files:
# NanoPlot-report.html  - Main interactive report
# NanoStats.txt         - Text summary
# LengthvsQualityScatterPlot_dot.png
# Non_weightedHistogramReadlength.png
# Non_weightedLogTransformed_HistogramReadlength.png
```

### 2.4 View the Report

```bash
# Open HTML report in browser (only works if you have a GUI)
firefox results/qc/nanoplot_raw/F003_M_enclense/NanoPlot-report.html

# Or view text summary
cat results/qc/nanoplot_raw/F003_M_enclense/NanoStats.txt
```

### 2.5 Interpret Key Metrics

Open the text summary and look for these metrics:

```bash
cat results/qc/nanoplot_raw/F003_M_enclense/NanoStats.txt
```

**Critical statistics:**

| Metric | What it means | Good value | Your value |
|--------|---------------|------------|------------|
| **Number of reads** | Total read count | 50,000-200,000 | _________ |
| **Total bases** | Total sequencing yield | >400 Mb | _________ |
| **Mean read length** | Average read size | >5,000 bp | _________ |
| **Mean read quality** | Average Phred score | >Q10 | _________ |
| **Read length N50** | 50% of bases in reads ≥ this length | >10,000 bp | _________ |

**Calculating coverage:**
```bash
# Formula: Coverage = Total bases / Genome size
# For Microbacterium (3.5 Mb genome):
echo "scale=1; <total_bases> / 3500000" | bc

# Example: 500,000,000 bases / 3,500,000 = 142.8x coverage
```

---

## Step 3: Generate Detailed Statistics with NanoStat

### 3.1 Run NanoStat

NanoStat provides similar information to NanoPlot but in text format only (faster):

```bash
NanoStat --fastq data/raw_reads/F003_M_enclense.fastq.gz \
    --name F003_M_enclense_raw_stats.txt \
    -t 4 \
    > results/qc/F003_M_enclense_raw_stats.txt
```

### 3.2 Extract Specific Metrics

```bash
# View entire stats file
cat results/qc/F003_M_enclense_raw_stats.txt

# Extract key metrics only
grep -E "Number of reads|Total bases|Mean read length|Mean read quality|Read length N50" \
    results/qc/F003_M_enclense_raw_stats.txt
```

---

## Step 4: Adapter Trimming with Porechop_ABI

### 4.1 Understanding Adapters

**Why remove adapters?**
- ONT sequencing uses adapter sequences to capture DNA
- Adapters can appear at read ends or in the middle (chimeric reads)
- Keeping adapters can interfere with assembly

**What Porechop_ABI does:**
- Detects ONT adapter sequences
- Removes adapters from read ends
- Splits chimeric reads (adapters in the middle)
- Discards reads with internal adapters (option `--discard_middle`)

### 4.2 Run Porechop_ABI

```bash
# Create output directory
mkdir -p results/qc/trimmed

# Run adapter trimming
porechop_abi \
    -i data/raw_reads/F003_M_enclense.fastq.gz \
    -o results/qc/trimmed/F003_M_enclense_trimmed.fastq.gz \
    --threads 4\
    -abi \
    --discard_middle
```

**Parameter explanation:**
- `-i`: Input FASTQ file
- `-o`: Output FASTQ file (will be compressed with .gz)
- `--threads 4`: Number of CPU threads
- `-abi`: Use the ABI algorithm (improved adapter detection)
- `--discard_middle`: Remove chimeric reads (adapters in the middle)

### 4.3 Monitor Progress

Porechop_ABI will display progress:

```
Finding adapters...
Trimming adapters from read ends...
Splitting reads with internal adapters...
1,234 / 50,000 reads had adapters trimmed (2.5%)
123 chimeric reads discarded
```

**💡 Note:** The number of reads with adapters is typically low (1-5%) for good quality ONT data.

### 4.4 Compare Before/After

```bash
# Count reads before
echo "Raw reads:"
zcat data/raw_reads/F003_M_enclense.fastq.gz | wc -l | awk '{print $1/4}'

# Count reads after
echo "Trimmed reads:"
zcat results/qc/trimmed/F003_M_enclense_trimmed.fastq.gz | wc -l | awk '{print $1/4}'

# Note: Some reads may be discarded if they had internal adapters
```

---

## Step 5: Read Filtering with Filtlong

### 5.1 Understanding Filtlong

**Purpose:** Reduce dataset to optimal reads for assembly

**Filtlong filtering strategies:**
1. **Length-based:** Remove very short reads (<1 kb)
2. **Quality-based:** Prioritize high-quality reads
3. **Coverage-based:** Subsample to target coverage (e.g., 100x)

**Why filter?**
- Speeds up assembly (less data to process)
- Removes poor-quality reads that add noise
- Longer reads are more valuable for assembly

### 5.2 Calculate Target Bases

For optimal bacterial assembly, target **50-100x coverage**:

```bash
# Define your genome size
# Microbacterium: ~3.5 Mb = 3,500,000 bp
# Alteromonas: ~4.5 Mb = 4,500,000 bp

GENOME_SIZE=3500000  # For Microbacterium enclense
TARGET_COV=100
TARGET_BASES=$((GENOME_SIZE * TARGET_COV))

echo "Target bases: ${TARGET_BASES}"
# Output: 350,000,000 (350 Mb)
```

### 5.3 Run Filtlong

```bash
# Create output directory
mkdir -p results/qc/filtered

# Run filtering
filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    --target_bases 350000000 \
    --length_weight 10 \
    --mean_q_weight 10 \
    results/qc/trimmed/F003_M_enclense_trimmed.fastq.gz \
    | gzip > results/qc/filtered/F003_M_enclense_filtered.fastq.gz
```

**Parameter explanation:**
- `--min_length 1000`: Discard reads shorter than 1 kb
- `--keep_percent 90`: Keep best 90% of bases (removes worst 10%)
- `--target_bases 350000000`: Aim for ~350 Mb (100x coverage)
- `--length_weight 10`: Prioritize longer reads
- `--mean_q_weight 10`: Also consider quality (balanced with length)
- Input: Trimmed reads
- Output: Filtered reads (piped through gzip for compression)

### 5.4 Understanding Weight Parameters

**Length vs. Quality weighting:**
- `--length_weight 10 --mean_q_weight 10`: Balanced (default)
- `--length_weight 20 --mean_q_weight 5`: Favor length over quality
- `--length_weight 5 --mean_q_weight 20`: Favor quality over length

**For bacterial assembly:** Length is more important! Use balanced or length-heavy weighting.

### 5.5 Alternative Filtering Strategies

**Scenario 1: Low coverage (<50x) - Keep everything above threshold**
```bash
filtlong \
    --min_length 1000 \
    --min_mean_q 7 \
    results/qc/trimmed/F003_M_enclense_trimmed.fastq.gz \
    | gzip > results/qc/filtered/F003_M_enclense_filtered.fastq.gz
```

**Scenario 2: Very high coverage (>200x) - More aggressive filtering**
```bash
filtlong \
    --min_length 2000 \
    --keep_percent 80 \
    --target_bases 350000000 \
    results/qc/trimmed/F003_M_enclense_trimmed.fastq.gz \
    | gzip > results/qc/filtered/F003_M_enclense_filtered.fastq.gz
```

---

## Step 6: Post-Filtering Quality Assessment

### 6.1 Run NanoPlot on Filtered Data

```bash
NanoPlot -t 4 \
    --fastq results/qc/filtered/F003_M_enclense_filtered.fastq.gz \
    --loglength \
    --plots dot \
    -o results/qc/nanoplot_filtered/F003_M_enclense \
    --N50 \
    --title "F003_M_enclense - Filtered Reads"
```

### 6.2 Generate Filtered Statistics

```bash
NanoStat --fastq results/qc/filtered/F003_M_enclense_filtered.fastq.gz \
    --name F003_M_enclense_filtered_stats.txt \
    -t 4 \
    > results/qc/F003_M_enclense_filtered_stats.txt
```

### 6.3 Compare Raw vs. Filtered

Create a side-by-side comparison:

```bash
echo "=== RAW DATA ==="
grep -E "Number of reads|Total bases|Mean read length|Mean read quality|Read length N50" \
    results/qc/F003_M_enclense_raw_stats.txt

echo ""
echo "=== FILTERED DATA ==="
grep -E "Number of reads|Total bases|Mean read length|Mean read quality|Read length N50" \
    results/qc/F003_M_enclense_filtered_stats.txt
```

### 6.4 Expected Changes

| Metric | Expected change |
|--------|----------------|
| Number of reads | ⬇️ Decreased (subsampled) |
| Total bases | ⬇️ Decreased (to target coverage) |
| Mean read length | ⬆️ Increased (removed short reads) |
| Mean read quality | ⬆️ Increased (removed low quality) |
| Read length N50 | ⬆️ Increased (kept longest reads) |

---

## Step 7: Validate Filtered Data

### 7.1 Check File Integrity

```bash
# Verify FASTQ structure (lines should be divisible by 4)
LINES=$(zcat results/qc/filtered/F003_M_enclense_filtered.fastq.gz | wc -l)
echo "Total lines: $LINES"
echo "Divisible by 4: $(( LINES % 4 == 0 ))"  # Should output: 1 (true)
```

### 7.2 Verify Quality Filtering

```bash
# Check minimum read length
echo "Shortest read length:"
zcat results/qc/filtered/F003_M_enclense_filtered.fastq.gz | \
    awk 'NR%4==2 {print length}' | sort -n | head -n 1

# Should be ≥1000 (our --min_length threshold)
```

### 7.3 Calculate Final Coverage

```bash
# Extract total bases from stats
TOTAL_BASES=$(grep "Total bases:" results/qc/F003_M_enclense_filtered_stats.txt | awk '{print $3}')
GENOME_SIZE=3500000

echo "scale=1; $TOTAL_BASES / $GENOME_SIZE" | bc
# Should be close to 100x
```

---

## Step 8: Quality Control for All Samples

### 8.1 Process Additional Samples

Repeat Steps 2-7 for each sample:

```bash
# Example for second sample
SAMPLE="F003_M_ginsengisoli"
GENOME_SIZE=3800000
TARGET_BASES=$((GENOME_SIZE * 100))

# Step 2: Raw QC
NanoPlot -t 4 \
    --fastq data/raw_reads/${SAMPLE}.fastq.gz \
    --loglength --plots dot \
    -o results/qc/nanoplot_raw/${SAMPLE} \
    --N50 --title "${SAMPLE} - Raw Reads"

# Step 3: Raw stats
NanoStat --fastq data/raw_reads/${SAMPLE}.fastq.gz \
    -t 4 > results/qc/${SAMPLE}_raw_stats.txt

# Step 4: Adapter trimming
porechop_abi \
    -i data/raw_reads/${SAMPLE}.fastq.gz \
    -o results/qc/trimmed/${SAMPLE}_trimmed.fastq.gz \
    --threads 8 -abi --discard_middle

# Step 5: Filtering
filtlong \
    --min_length 1000 --keep_percent 90 \
    --target_bases ${TARGET_BASES} \
    --length_weight 10 --mean_q_weight 10 \
    results/qc/trimmed/${SAMPLE}_trimmed.fastq.gz \
    | gzip > results/qc/filtered/${SAMPLE}_filtered.fastq.gz

# Step 6: Filtered QC
NanoPlot -t 4 \
    --fastq results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    --loglength --plots dot \
    -o results/qc/nanoplot_filtered/${SAMPLE} \
    --N50 --title "${SAMPLE} - Filtered Reads"

NanoStat --fastq results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    -t 4 > results/qc/${SAMPLE}_filtered_stats.txt
```

### 8.2 Create Summary Table

```bash
# Create a comparison table
echo -e "Sample\tRaw_Reads\tRaw_Bases\tFiltered_Reads\tFiltered_Bases\tCoverage" > qc_summary.tsv

for sample in F003_M_enclense F003_M_ginsengisoli F003_A_portus; do
    RAW_READS=$(grep "Number of reads:" results/qc/${sample}_raw_stats.txt | awk '{print $4}')
    RAW_BASES=$(grep "Total bases:" results/qc/${sample}_raw_stats.txt | awk '{print $3}')
    FILT_READS=$(grep "Number of reads:" results/qc/${sample}_filtered_stats.txt | awk '{print $4}')
    FILT_BASES=$(grep "Total bases:" results/qc/${sample}_filtered_stats.txt | awk '{print $3}')
    
    # Calculate coverage (adjust genome size per species)
    COV=$(echo "scale=1; $FILT_BASES / 3500000" | bc)
    
    echo -e "${sample}\t${RAW_READS}\t${RAW_BASES}\t${FILT_READS}\t${FILT_BASES}\t${COV}x" >> qc_summary.tsv
done

# View the table
column -t qc_summary.tsv
```

---

## Common Issues and Troubleshooting

### Issue 1: Very Low Read Count After Filtering

**Symptoms:** <10,000 reads after filtering

**Diagnosis:**
```bash
# Check raw data yield
grep "Total bases:" results/qc/F003_M_enclense_raw_stats.txt
```

**Solutions:**
- If raw yield is low (<200 Mb): Insufficient sequencing, may need to resequence
- If raw yield is good but filtered is low: Adjust `--keep_percent` to 95 or remove `--target_bases`

### Issue 2: Very Short Mean Read Length (<2 kb)

**Symptoms:** Mean read length <2000 bp after filtering

**Possible causes:**
- Poor DNA extraction (fragmented DNA)
- Sequencing issues

**Solutions:**
- Lower `--min_length` to 500
- Check if longer reads exist in raw data (view NanoPlot histogram)

### Issue 3: Low Mean Quality (<Q8)

**Symptoms:** Mean quality below Q8

**Possible causes:**
- Old flow cell
- Sequencing chemistry issues

**Solutions:**
- Can still assemble! ONT error correction (Medaka) handles this
- If very low (<Q5), consider resequencing

### Issue 4: Filtlong Removes Too Much Data

**Symptoms:** <30x coverage after filtering

**Solution:**
```bash
# More lenient filtering
filtlong \
    --min_length 1000 \
    --min_mean_q 7 \
    --keep_percent 95 \
    results/qc/trimmed/sample_trimmed.fastq.gz \
    | gzip > results/qc/filtered/sample_filtered.fastq.gz
```

---

## Understanding Your Results

### What Makes Good ONT Data for Assembly?

✅ **Excellent:**
- Mean read length >10 kb
- Read length N50 >15 kb
- Mean quality >Q12
- Total coverage 80-150x

✅ **Good (adequate for assembly):**
- Mean read length >5 kb
- Read length N50 >8 kb
- Mean quality >Q10
- Total coverage 50-100x

⚠️ **Marginal (assembly possible but may be fragmented):**
- Mean read length >3 kb
- Read length N50 >5 kb
- Mean quality >Q8
- Total coverage 30-50x

❌ **Poor (difficult to assemble):**
- Mean read length <3 kb
- Mean quality <Q8
- Total coverage <30x

---

## Exercise: Analyze Your Data

### Task 1: Quality Assessment

Fill in this table with your sample's statistics:

| Metric | Raw Data | Filtered Data |
|--------|----------|---------------|
| Number of reads | | |
| Total bases | | |
| Mean read length | | |
| Mean quality | | |
| Read length N50 | | |
| Estimated coverage | | |

### Task 2: Interpretation Questions

1. **How much data was removed by filtering?**
   - Calculate: `(Raw bases - Filtered bases) / Raw bases × 100%`

2. **Did mean read length increase after filtering?**
   - If yes, why? If no, why not?

3. **Is your coverage adequate for assembly?**
   - Target: 50-100x for good assembly

4. **What percentage of reads were shorter than 1 kb?**
   - Hint: Compare total reads before/after `--min_length 1000`

### Task 3: Compare Samples

If you have multiple samples:

1. Which sample has the highest quality?
2. Which sample has the longest reads?
3. Do you see differences between species or strains?

---

## Next Steps

Once you have high-quality filtered reads:

✅ Proceed to **Tutorial 2: Genome Assembly**

**Files you'll need:**
- `results/qc/filtered/${SAMPLE}_filtered.fastq.gz`

**What you'll do:**
- De novo assembly with Flye
- Error correction with Medaka
- Assembly quality assessment

---

## Quick Reference Commands

### Essential QC Workflow
```bash
# 1. Raw QC
NanoPlot -t 4 --fastq raw.fastq.gz -o qc_raw/ --loglength --N50

# 2. Adapter trimming
porechop_abi -i raw.fastq.gz -o trimmed.fastq.gz --threads 8 -abi --discard_middle

# 3. Read filtering
filtlong --min_length 1000 --keep_percent 90 --target_bases 350000000 \
    trimmed.fastq.gz | gzip > filtered.fastq.gz

# 4. Filtered QC
NanoPlot -t 4 --fastq filtered.fastq.gz -o qc_filtered/ --loglength --N50
```

### Useful One-Liners
```bash
# Count reads
zcat file.fastq.gz | wc -l | awk '{print $1/4}'

# Calculate coverage
echo "scale=1; <total_bases> / <genome_size>" | bc

# Get shortest read
zcat file.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | head -1

# Get longest read
zcat file.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | tail -1

# Calculate N50 (requires seqkit)
seqkit stats -a file.fastq.gz
```

---

**Continue to:** [Tutorial 2: Genome Assembly](TUTORIAL_02_ASSEMBLY.md)
