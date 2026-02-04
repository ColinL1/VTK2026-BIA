# Tutorial 2: Manual Genome Assembly

## Learning Objectives

By the end of this tutorial, you will be able to:
- Perform de novo genome assembly using Flye
- Understand assembly graphs and repeat resolution
- Polish assemblies using Medaka
- Assess assembly quality with QUAST
- Evaluate completeness with BUSCO
- Visualize assembly graphs with Bandage

---

## Prerequisites

### Required Software

```bash
# Activate the assembly conda environment
conda activate assembly

# Verify tools are installed
flye --version
medaka --version
quast.py --version
busco --version
```

### Required Input

You must have completed Tutorial 1 (Quality Control). You need:
```bash
# Filtered reads from QC
ls results/qc/filtered/F003_M_enclense_filtered.fastq.gz
```

### Directory Setup

```bash
# Navigate to project directory
cd ~/'MY-NAME-EXAMPLE'/ # change this to your own name. 

# Create assembly directories
mkdir -p results/assembly/F003_M_enclense/flye
mkdir -p results/assembly/F003_M_enclense/medaka
mkdir -p results/assembly/F003_M_enclense/qc
```

---

## Background: De Novo Assembly

### What is De Novo Assembly?

**De novo assembly** reconstructs the genome sequence from sequencing reads **without using a reference genome**. This is essential for:
- Novel organisms (no reference available)
- Discovering structural variants
- Identifying unique genomic features

### Why No Reference for Our Samples?

Our coral-associated bacteria (*Microbacterium enclense*, *M. ginsengisoli*, *Alteromonas portus*) are:
- Lack high-quality complete genome references
- May have unique genomic features

### Long Reads vs. Short Reads

**ONT (Long Reads):**
- ✅ Span repetitive regions
- ✅ Resolve complex genomic structures
- ✅ Produce more contiguous assemblies (higher N50)
- ⚠️ Higher error rate (~5-10%)

**Illumina (Short Reads):**
- ✅ High accuracy (>99.9%)
- ⚠️ Cannot span long repeats
- ⚠️ Fragmented assemblies (many small contigs)

### Assembly Challenges

**Repeats:** Identical sequences appearing multiple times
- rRNA operons (typically 3-9 copies of 16S, 23S, 5S)
- Transposons and insertion sequences
- Tandem repeats

**Solution:** Long reads can span repeats, allowing correct placement.

---

## Step 1: De Novo Assembly with Flye

### 1.1 Understanding Flye

**Flye** is a de novo assembler optimized for long reads (ONT and PacBio).

**Algorithm:**
1. **Disjointig assembly:** Creates initial contigs from read overlaps
2. **Repeat graph construction:** Identifies and represents repeats
3. **Consensus calling:** Generates consensus sequence
4. **Repeat resolution:** Attempts to resolve repetitive regions
5. **Polishing:** Multiple rounds of error correction

**Key advantages:**
- Handles high error rates in long reads
- Excellent repeat resolution
- Meta-mode available for metagenomes

### 1.2 Estimate Genome Size

Before assembly, estimate expected genome size:

**Microbacterium species:**
- *M. enclense*: ~3.5 Mb
- *M. ginsengisoli*: ~3.8 Mb

**Alteromonas species:**
- *A. portus*: ~4.5 Mb

**Why does Flye need genome size?**
- Helps determine coverage thresholds
- Guides repeat detection
- Optimizes memory usage

### 1.3 Run Flye Assembly

```bash
# Run Flye
flye \
    --nano-hq results/qc/filtered/F003_M_enclense_filtered.fastq.gz \
    --genome-size 3.5m \
    --threads 4 \
    --out-dir results/assembly/F003_M_enclense/flye \
    --iterations 3 \
    --meta
```

**Parameter explanation:**
- `--nano-hq`: High-quality ONT reads (Q10+, basecalled with SUP model)
  - Use `--nano-raw` for older, lower-quality reads
- `--genome-size 3.5m`: Expected genome size
  - Use `m` for megabases (3.5m = 3,500,000 bp)
  - Use `g` for gigabases if needed
- `--threads 4`: Number of CPU threads
- `--out-dir`: Output directory
- `--iterations 2`: Number of polishing iterations
  - Default is 1, 2-3 recommended for better accuracy
  - More iterations = longer runtime but better quality
- `--meta`: Enables meta-mode for uneven coverage
  - Helpful for plasmids and chromosomes with different coverage
  - Prevents over-collapsing of repeats

**⏱️ Runtime:** 2-3 hours depending on coverage and genome size

### 1.4 Monitor Progress

Flye provides real-time updates. Open a new terminal to watch:

```bash
# Watch the log file
tail -f results/assembly/F003_M_enclense/flye/flye.log
```

**Assembly stages you'll see:**
1. `[INFO] Starting Flye` - Initialization
2. `[INFO] Running Disjointig` - Creating initial contigs
3. `[INFO] Assembled __ disjointigs` - Initial assembly complete
4. `[INFO] Running repeat analysis` - Identifying repeats
5. `[INFO] Generating assembly graph` - Creating repeat graph
6. `[INFO] Resolving repeats` - Attempting to resolve repeats
7. `[INFO] Polishing` - Error correction (may run multiple iterations)
8. `[INFO] Assembly statistics` - Final summary

### 1.5 Understand Assembly Output Files

```bash
# List Flye output files
ls -lh results/assembly/F003_M_enclense/flye/

# Key files:
# assembly.fasta          - Final assembled sequences (contigs)
# assembly_graph.gfa      - Assembly graph (for visualization)
# assembly_graph.gv       - Graphviz format (alternative visualization)
# assembly_info.txt       - Contig statistics
# flye.log                - Detailed log
```

### 1.6 Examine Assembly Info

```bash
# View assembly statistics
cat results/assembly/F003_M_enclense/flye/assembly_info.txt
```

**Columns explained:**
- `#seq_name`: Contig identifier
- `length`: Contig length in base pairs
- `cov.`: Coverage (read depth)
- `circ.`: Whether contig is circular (Y/N)
  - `Y` indicates complete chromosome or plasmid!
- `repeat`: Whether contig contains unresolved repeats
- `mult.`: Multiplicity (copy number)
- `alt_group`: Alternative assembly paths
- `graph_path`: Path in assembly graph

**Example interpretation:**
```
#seq_name    length    cov.  circ.  repeat  mult.
contig_1     3456789   45    Y      N       1
contig_2     125000    42    Y      N       1
contig_3     45000     48    Y      N       1
```
- `contig_1`: 3.5 Mb, circular → Complete chromosome!
- `contig_2`: 125 kb, circular → Likely plasmid
- `contig_3`: 45 kb, circular → Smaller plasmid

### 1.7 View Assembly FASTA

```bash
# View contig headers
grep "^>" results/assembly/F003_M_enclense/flye/assembly.fasta

# Count contigs
grep -c "^>" results/assembly/F003_M_enclense/flye/assembly.fasta

# Get assembly size
grep -v "^>" results/assembly/F003_M_enclense/flye/assembly.fasta | \
    tr -d '\n' | wc -c
```

### 1.8 Quick Assembly Statistics

```bash
# Use seqkit for detailed stats
seqkit stats -a results/assembly/F003_M_enclense/flye/assembly.fasta

# Output:
# file      format  type  num_seqs  sum_len  min_len  avg_len  max_len  Q1      Q2      Q3      sum_gap  N50    Q20(%)  Q30(%)
```

**Key metrics:**
- `num_seqs`: Number of contigs (lower is better)
- `sum_len`: Total assembly size
- `max_len`: Largest contig (ideally the chromosome)
- `N50`: Half of assembly is in contigs ≥ this size

---

## Step 2: Visualize Assembly Graph with Bandage

### 2.1 Understanding Assembly Graphs

**Assembly graphs** show how contigs connect:
- **Nodes:** Contigs/sequences
- **Edges:** Overlaps between contigs
- **Circular paths:** Complete chromosomes/plasmids
- **Branches:** Repeats or ambiguous regions

### 2.2 Install Bandage (if needed)

```bash
# Install via conda
conda install -c bioconda bandage

# Or download from: https://rrwick.github.io/Bandage/
```

### 2.3 View Assembly Graph

```bash
# Generate PNG image from GFA file
Bandage image \
    results/assembly/F003_M_enclense/flye/assembly_graph.gfa \
    results/assembly/F003_M_enclense/assembly_graph.png \
    --height 2000

# View the image
eog results/assembly/F003_M_enclense/assembly_graph.png
```

**Interpreting the graph:**
- **Circular loops:** Complete sequences (chromosome, plasmids)
- **Linear segments:** Incomplete or low-coverage regions
- **Connected components:** Separate replicons
- **Node thickness:** Coverage depth

### 2.4 Interactive Exploration

```bash
# Open interactive Bandage GUI
Bandage load results/assembly/F003_M_enclense/flye/assembly_graph.gfa
```

**In the GUI:**
1. Click "Draw graph" to visualize
2. Click on nodes to see details (length, coverage)
3. Use "BLAST search" to find specific sequences (e.g., 16S rRNA)

---

## Step 3: Error Correction with Medaka

### 3.1 Why Polish?

**ONT error characteristics:**
- ~5-10% error rate (mostly insertions/deletions)
- Homopolymer errors (e.g., AAAA vs. AAAAA)
- Systematic biases

**Polishing improves:**
- Base-level accuracy (>99.5% after polishing)
- Correct reading frames for genes
- Variant calling accuracy

### 3.2 Understanding Medaka

**Medaka** is a neural network-based polisher trained on ONT data.

**How it works:**
1. Aligns reads to draft assembly
2. Generates pileups (read coverage at each position)
3. Neural network predicts consensus sequence
4. Corrects errors in draft assembly

**Key concept: Medaka models**
- Each model is trained for specific basecaller/flow cell combinations
- Using the wrong model reduces polishing accuracy
- Common models:
  - `r941_min_hac_g507` - R9.4.1 flow cell, HAC basecalling
  - `r1041_e82_400bps_sup_v5.2.0` - R10.4.1 flow cell, SUP basecalling

### 3.3 Determine Your Basecaller Model

Check your sequencing run metadata:

```bash
# If you have FAST5 files, check header
# For FASTQ files, sometimes in read headers
zcat data/raw_reads/F003_M_enclense.fastq.gz | head -1

# Look for information like:
# @read_id runid=... basecall_model_version_id=...
```

**Most likely for 2026 data:**
- Flow cell: R10.4.1 or newer
- Basecaller: Dorado SUP (super-accurate)
- Model: `r1041_e82_400bps_sup_v5.2.0` or `r1041_e82_400bps_sup_v4.2.0`

**If unsure:** Use `r1041_e82_400bps_sup_v4.2.0` (compatible with most recent data)

### 3.4 Copy Draft Assembly

```bash
# Copy Flye assembly for clarity
cp results/assembly/F003_M_enclense/flye/assembly.fasta \
   results/assembly/F003_M_enclense/F003_M_enclense_flye.fasta
```

### 3.5 Run Medaka Polishing

```bash
# Run medaka_consensus
medaka_consensus \
    -i results/qc/filtered/F003_M_enclense_filtered.fastq.gz \
    -d results/assembly/F003_M_enclense/F003_M_enclense_flye.fasta \
    -o results/assembly/F003_M_enclense/medaka \
    -t 4 \
    -m r1041_e82_400bps_sup_v5.2.0
```

**Parameter explanation:**
- `-i`: Input reads (original filtered FASTQ)
- `-d`: Draft assembly to polish
- `-o`: Output directory
- `-t`: Number of threads
- `-m`: Medaka model (match your basecaller!)

**⏱️ Runtime:** 30-60 minutes

### 3.6 Monitor Medaka Progress

```bash
# Medaka runs in stages:
# 1. mini_align - Aligns reads to draft
# 2. medaka consensus - Calls consensus
# 3. medaka stitch - Assembles polished sequences
```

### 3.7 Examine Polished Assembly

```bash
# Medaka output
ls results/assembly/F003_M_enclense/medaka/

# Key file: consensus.fasta - Polished assembly

# Copy to main directory
cp results/assembly/F003_M_enclense/medaka/consensus.fasta \
   results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta
```

### 3.8 Compare Draft vs. Polished

```bash
echo "=== Draft Assembly (Flye) ==="
seqkit stats results/assembly/F003_M_enclense/F003_M_enclense_flye.fasta

echo ""
echo "=== Polished Assembly (Flye + Medaka) ==="
seqkit stats results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta
```

**Expected:** Total length and contig count should be nearly identical.

---

## Step 4: Assembly Quality Assessment with QUAST

### 4.1 Understanding QUAST

**QUAST** (QUality ASsessment Tool) evaluates genome assemblies.

**Metrics provided:**
- Total length, N50, L50
- Number of contigs
- GC content
- Largest contig
- N's per 100 kb (assembly gaps)

### 4.2 Run QUAST

```bash
# Compare both assemblies
quast.py \
    results/assembly/F003_M_enclense/F003_M_enclense_flye.fasta \
    results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta \
    -o results/assembly/F003_M_enclense/qc/quast \
    --threads 4 \
    --min-contig 500 \
    --labels "Flye,Flye+Medaka"
```

**Parameter explanation:**
- Input: Both draft and polished assemblies (for comparison)
- `-o`: Output directory
- `--threads`: Number of threads
- `--min-contig 500`: Ignore very small contigs (<500 bp)
- `--labels`: Names for comparison

### 4.3 View QUAST Report

```bash
# Open HTML report
firefox results/assembly/F003_M_enclense/qc/quast/report.html

# Or view text report
cat results/assembly/F003_M_enclense/qc/quast/report.txt
```

### 4.4 Interpret Key Metrics

**Example QUAST output:**
```
Assembly                  Flye      Flye+Medaka
# contigs (>= 0 bp)       5         5
# contigs (>= 1000 bp)    5         5
Total length (>= 0 bp)    3625480   3625325
Total length (>= 1000 bp) 3625480   3625325
# contigs                 5         5
Largest contig            3456789   3456789
Total length              3625480   3625325
GC (%)                    69.45     69.45
N50                       3456789   3456789
N90                       125000    125000
L50                       1         1
L90                       1         1
# N's per 100 kbp         0.00      0.00
```

**What to look for:**

✅ **Excellent assembly (bacterial genome):**
- Total length: 3-5 Mb (matches expected genome size)
- N50: >500 kb (ideally >1 Mb for complete chromosome)
- L50: 1-2 (one chromosome, maybe plasmids)
- Number of contigs: <10
- GC content: 40-75% (species-dependent)

✅ **Good assembly:**
- Total length: Within 10% of expected
- N50: >100 kb
- L50: <5
- Number of contigs: <20

⚠️ **Fragmented assembly:**
- N50: <50 kb
- Many contigs (>50)
- Total length much smaller or larger than expected

### 4.5 Understanding N50 and L50

**N50 definition:**
- Sort all contigs by length (largest to smallest)
- N50 = length of contig where 50% of total bases are in contigs ≥ this length

**Example:**
```
Contigs: 3,000,000 bp, 200,000 bp, 100,000 bp, 50,000 bp, 25,000 bp
Total: 3,375,000 bp
50% = 1,687,500 bp

Cumulative:
3,000,000 (>1,687,500) ← N50 = 3,000,000 bp
```

**L50:** Number of contigs needed to reach 50% of total length (lower is better)

---

## Step 5: Completeness Assessment with BUSCO

### 5.1 Understanding BUSCO

**BUSCO** (Benchmarking Universal Single-Copy Orthologs) assesses genome completeness.

**How it works:**
- Searches for conserved genes expected in all bacteria
- Database: `bacteria_odb10` (124 conserved genes)
- Reports: Complete, Duplicated, Fragmented, Missing

**Why important:**
- Assembly might look good (N50, contigs) but be incomplete
- BUSCO checks if expected genes are present

### 5.2 Download BUSCO Database

```bash
# Download bacteria lineage dataset (one-time setup)
busco --list-datasets

# Download specific dataset
busco --download bacteria_odb10

# Database installed to: ~/.busco_downloads/lineages/bacteria_odb10/
```

### 5.3 Run BUSCO

```bash
# Run BUSCO on polished assembly
busco \
    -i results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta \
    -o F003_M_enclense_busco \
    --out_path results/assembly/F003_M_enclense/qc \
    -m genome \
    -l bacteria_odb10 \
    --cpu 8 \
    --offline
```

**Parameter explanation:**
- `-i`: Input assembly
- `-o`: Output folder name
- `--out_path`: Directory to create output folder
- `-m genome`: Mode (genome, transcriptome, or proteins)
- `-l bacteria_odb10`: Lineage dataset
- `--cpu`: Number of threads
- `--offline`: Don't check for updates (faster)

**⏱️ Runtime:** 10-30 minutes

### 5.4 View BUSCO Results

```bash
# Find the summary file
BUSCO_SUMMARY=$(find results/assembly/F003_M_enclense/qc/F003_M_enclense_busco/ \
    -name "short_summary*.txt")

cat $BUSCO_SUMMARY
```

**Example output:**
```
***** Results: *****

C:98.4%[S:98.4%,D:0.0%],F:1.6%,M:0.0%,n:124

122     Complete BUSCOs (C)
122     Complete and single-copy BUSCOs (S)
0       Complete and duplicated BUSCOs (D)
2       Fragmented BUSCOs (F)
0       Missing BUSCOs (M)
124     Total BUSCO groups searched
```

### 5.5 Interpret BUSCO Results

**Categories:**
- **C (Complete):** Gene found in full length
  - **S (Single-copy):** Found once (expected for bacteria)
  - **D (Duplicated):** Found multiple times
- **F (Fragmented):** Partial gene found (possible assembly break)
- **M (Missing):** Gene not found

**Quality thresholds:**

✅ **Excellent:**
- Complete (C): >95%
- Duplicated (D): <2%
- Missing (M): <2%

✅ **Good:**
- Complete (C): >90%
- Duplicated (D): <5%
- Missing (M): <5%

⚠️ **Concerning:**
- Complete (C): <85%
- Duplicated (D): >10%
- Missing (M): >10%

**High duplication (>10%) may indicate:**
- Contamination (multiple organisms)
- Misassembly
- True gene duplication (rare)

**High missing (>10%) may indicate:**
- Incomplete assembly
- Low coverage regions
- True gene loss (unusual)

### 5.6 Examine Detailed Results

```bash
# See which genes are missing/fragmented
cat results/assembly/F003_M_enclense/qc/F003_M_enclense_busco/run_bacteria_odb10/missing_busco_list.txt

cat results/assembly/F003_M_enclense/qc/F003_M_enclense_busco/run_bacteria_odb10/fragmented_busco_list.txt
```

---

## Step 6: Generate Assembly Summary

### 6.1 Create Summary Report

```bash
# Create a summary file
SAMPLE="F003_M_enclense"
SUMMARY="results/assembly/${SAMPLE}/${SAMPLE}_assembly_summary.txt"

cat > $SUMMARY <<EOF
========================================
Assembly Summary for ${SAMPLE}
========================================

EOF

# Add QUAST results
echo "=== QUAST Statistics (Polished Assembly) ===" >> $SUMMARY
grep -A 15 "Flye+Medaka" results/assembly/${SAMPLE}/qc/quast/report.txt >> $SUMMARY

echo "" >> $SUMMARY
echo "=== BUSCO Completeness ===" >> $SUMMARY
grep "C:" results/assembly/${SAMPLE}/qc/${SAMPLE}_busco/short_summary*.txt >> $SUMMARY

echo "" >> $SUMMARY
echo "=== Flye Assembly Info ===" >> $SUMMARY
cat results/assembly/${SAMPLE}/flye/assembly_info.txt >> $SUMMARY

# View summary
cat $SUMMARY
```

### 6.2 Compare Multiple Samples

```bash
# Create comparison table
echo -e "Sample\tSize\tContigs\tN50\tBUSCO_Complete" > assembly_comparison.tsv

for SAMPLE in F003_M_enclense F003_M_ginsengisoli F003_A_portus; do
    SIZE=$(grep "Total length" results/assembly/${SAMPLE}/qc/quast/report.txt | \
        grep "Flye+Medaka" | awk '{print $NF}')
    CONTIGS=$(grep "# contigs (" results/assembly/${SAMPLE}/qc/quast/report.txt | \
        grep "Flye+Medaka" | head -1 | awk '{print $NF}')
    N50=$(grep "N50" results/assembly/${SAMPLE}/qc/quast/report.txt | \
        grep "Flye+Medaka" | awk '{print $NF}')
    BUSCO=$(grep "C:" results/assembly/${SAMPLE}/qc/${SAMPLE}_busco/short_summary*.txt | \
        grep -o "C:[0-9.]*%" | cut -d: -f2)
    
    echo -e "${SAMPLE}\t${SIZE}\t${CONTIGS}\t${N50}\t${BUSCO}" >> assembly_comparison.tsv
done

# View table
column -t assembly_comparison.tsv
```

---

## Step 7: Assembly for All Samples

### 7.1 Process Additional Samples

**Microbacterium ginsengisoli:**
```bash
SAMPLE="F003_M_ginsengisoli"
GENOME_SIZE="3.8m"

flye --nano-hq results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    --genome-size ${GENOME_SIZE} --threads 4 \
    --out-dir results/assembly/${SAMPLE}/flye \
    --iterations 2 --meta

cp results/assembly/${SAMPLE}/flye/assembly.fasta \
   results/assembly/${SAMPLE}/${SAMPLE}_flye.fasta

medaka_consensus \
    -i results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    -d results/assembly/${SAMPLE}/${SAMPLE}_flye.fasta \
    -o results/assembly/${SAMPLE}/medaka \
    -t 8 -m r1041_e82_400bps_sup_v5.2.0

cp results/assembly/${SAMPLE}/medaka/consensus.fasta \
   results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta
```

**Alteromonas portus:**
```bash
SAMPLE="F003_A_portus"
GENOME_SIZE="4.5m"  # Larger genome

flye --nano-hq results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    --genome-size ${GENOME_SIZE} --threads 4 \
    --out-dir results/assembly/${SAMPLE}/flye \
    --iterations 2 --meta

cp results/assembly/${SAMPLE}/flye/assembly.fasta \
   results/assembly/${SAMPLE}/${SAMPLE}_flye.fasta

medaka_consensus \
    -i results/qc/filtered/${SAMPLE}_filtered.fastq.gz \
    -d results/assembly/${SAMPLE}/${SAMPLE}_flye.fasta \
    -o results/assembly/${SAMPLE}/medaka \
    -t 8 -m r1041_e82_400bps_sup_v5.2.0

cp results/assembly/${SAMPLE}/medaka/consensus.fasta \
   results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta
```

### 7.2 Run QC for All Samples

```bash
# QUAST for all samples
for SAMPLE in F003_M_enclense F003_M_ginsengisoli F003_A_portus; do
    quast.py \
        results/assembly/${SAMPLE}/${SAMPLE}_flye.fasta \
        results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta \
        -o results/assembly/${SAMPLE}/qc/quast \
        --threads 4 --min-contig 500 \
        --labels "Flye,Flye+Medaka"
done

# BUSCO for all samples
for SAMPLE in F003_M_enclense F003_M_ginsengisoli F003_A_portus; do
    busco -i results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta \
        -o ${SAMPLE}_busco \
        --out_path results/assembly/${SAMPLE}/qc \
        -m genome -l bacteria_odb10 --cpu 8 --offline
done
```

---

## Common Issues and Troubleshooting

### Issue 1: Very Fragmented Assembly (>50 contigs)

**Diagnosis:**
```bash
# Check coverage
grep "cov." results/assembly/sample/flye/assembly_info.txt
```

**Possible causes:**
1. **Low coverage** (<30x)
   - Solution: Reduce filtering stringency in QC step
2. **Repetitive genome**
   - Solution: Expected for some organisms, acceptable if BUSCO >90%
3. **Contamination**
   - Solution: Check for multiple species in assembly graph

### Issue 2: Assembly Much Larger Than Expected

**Example:** Expected 3.5 Mb, got 5 Mb

**Possible causes:**
1. **Contamination** (multiple organisms)
   - Check assembly graph for disconnected components
   - Run `kraken2` for taxonomic classification
2. **Many plasmids**
   - Check `assembly_info.txt` for multiple circular contigs
3. **Duplicated regions**
   - High BUSCO duplication (>10%)

### Issue 3: Low BUSCO Completeness (<80%)

**Diagnosis:**
```bash
# Check QUAST total length
grep "Total length" results/assembly/sample/qc/quast/report.txt
```

**Solutions:**
1. If total length is low → Incomplete assembly
   - Increase coverage in filtering step
   - Check for low-coverage regions
2. If total length is correct → Wrong BUSCO lineage
   - Try `actinobacteria_odb10` for *Microbacterium*
   - Try `proteobacteria_odb10` for *Alteromonas*

### Issue 4: Medaka Crashes

**Error:** "Could not find model..."

**Solution:**
```bash
# List available models
medaka tools list_models

# Use closest match to your basecaller
# For R10.4.1 SUP: r1041_e82_400bps_sup_v4.2.0
```

### Issue 5: No Circular Contigs

**Diagnosis:**
```bash
grep "circ." results/assembly/sample/flye/assembly_info.txt
# All show "N"
```

**Possible causes:**
- Insufficient coverage at chromosome ends
- Assembly breaks at repeats
- Linear chromosomes (rare in bacteria)

**Note:** Not always a problem! Many good assemblies aren't perfectly circular.

---

## Exercise: Analyze Your Assembly

### Task 1: Assembly Metrics

Fill in this table for your sample:

| Metric | Value |
|--------|-------|
| Total assembly size | |
| Number of contigs | |
| Largest contig | |
| N50 | |
| L50 | |
| GC content | |
| BUSCO Complete (%) | |
| Circular contigs | |

### Task 2: Interpretation

1. **Does assembly size match expected genome size?**
   - *Microbacterium*: 3-4 Mb
   - *Alteromonas*: 4-5 Mb

2. **How contiguous is your assembly?**
   - Calculate: `Largest contig / Total size × 100%`
   - >90% = very good (complete chromosome)

3. **Is the assembly complete?**
   - BUSCO >90% = yes
   - BUSCO 80-90% = mostly complete
   - BUSCO <80% = incomplete

4. **How many replicons do you have?**
   - Count circular contigs in `assembly_info.txt`
   - 1 = chromosome only
   - 2-4 = chromosome + plasmids

### Task 3: Visualization

1. Open your assembly graph in Bandage
2. Identify the chromosome (largest circular component)
3. Identify plasmids (smaller circular components)
4. Take a screenshot and annotate key features

---

## Next Steps

Once you have high-quality assemblies:

✅ Proceed to **Tutorial 3: Genome Annotation**

**Files you'll need:**
- `results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta`

**What you'll do:**
- Predict genes with Bakta
- Assign functions
- Generate annotation files (GFF, GenBank)

---

## Quick Reference Commands

### Essential Assembly Workflow
```bash
# 1. Flye assembly
flye --nano-hq filtered.fastq.gz --genome-size 3.5m --threads 4 \
    --out-dir flye_out/ --iterations 2 --meta

# 2. Medaka polishing
medaka_consensus -i filtered.fastq.gz -d draft.fasta -o medaka_out/ \
    -t 8 -m r1041_e82_400bps_sup_v5.2.0

# 3. QUAST assessment
quast.py draft.fasta polished.fasta -o quast_out/ --threads 4

# 4. BUSCO completeness
busco -i polished.fasta -o busco_out -m genome -l bacteria_odb10 --cpu 8
```

### Useful One-Liners
```bash
# Assembly size
grep -v "^>" assembly.fasta | tr -d '\n' | wc -c

# Number of contigs
grep -c "^>" assembly.fasta

# Get largest contig size
bioawk -c fastx '{print $name, length($seq)}' assembly.fasta | sort -k2 -rn | head -1

# Calculate N50 (requires seqkit)
seqkit stats -a assembly.fasta | cut -f18
```

---

**Continue to:** [Tutorial 3: Genome Annotation](TUTORIAL_03_ANNOTATION.md)
