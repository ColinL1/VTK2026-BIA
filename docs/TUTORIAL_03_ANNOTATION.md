# Tutorial 3: Manual Genome Annotation

## Learning Objectives

By the end of this tutorial, you will be able to:
- Perform functional annotation using Bakta
- Understand annotation file formats (GFF, GenBank, FASTA)
- Analyze gene content and functional categories
- Extract and interpret annotation statistics
- Compare annotations across samples

---

## Prerequisites

### Required Software

```bash
# Activate annotation environment
conda activate VTK2026_Annotate

# Verify Bakta is installed
bakta --version
```

### Required Input

You must have completed Tutorial 2 (Assembly). You need:
```bash
# Polished assembly from assembly step
ls results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta
```

### Bakta Database Setup

Bakta requires a reference database. This is a one-time download (~30 GB).

```bash
# Check if database exists
if [ -z "$BAKTA_DB" ]; then
    echo "BAKTA_DB not set!"
    echo "Download database with:"
    echo "  bakta_db download --output ~/databases/bakta_db --type full"
fi

# Download Bakta database (one-time, ~1-2 hours)
mkdir -p ~/databases
bakta_db download --output ~/databases/bakta_db --type full

# Set environment variable (add to ~/.bashrc for permanence)
export BAKTA_DB=~/databases/bakta_db
```

**Database types:**
- `full`: Complete database (~30 GB) - **recommended**
- `light`: Smaller version (~3 GB) - fewer annotations

### Directory Setup

```bash
# Navigate to project directory
cd ~/'MY-NAME-EXAMPLE'/ # change this to your own name. and work in here from now on. 

# Create annotation directories
mkdir -p results/annotation/F003_M_enclense
```

---

## Background: Genome Annotation

### What is Genome Annotation?

**Annotation** is the process of identifying biological features in a genome:
1. **Structural annotation:** Where are genes located?
2. **Functional annotation:** What do those genes do?

### Why Annotation is Challenging

**The problem:**
- Genomes are just sequences: `ATGCGATCG...`
- We need to predict:
  - Gene boundaries (start/stop codons)
  - Gene function (what protein it encodes)
  - Regulatory elements (promoters, terminators)

### Prokka vs. Bakta

**Prokka** (traditional tool):
- Fast, widely used
- Relies on species-specific databases
- ⚠️ Poor performance for novel/understudied organisms
- High percentage of "hypothetical proteins"

**Bakta** (modern tool):
- Comprehensive, taxon-independent database
- Better functional annotation
- Includes COG, GO terms, EC numbers
- ✅ Better for our coral-associated bacteria!

### What Bakta Predicts

**Gene features:**
- CDS (coding sequences) - protein-coding genes
- tRNA - transfer RNA genes
- rRNA - ribosomal RNA genes (16S, 23S, 5S)
- ncRNA - other non-coding RNAs
- CRISPR arrays - adaptive immune systems
- Regulatory elements

**Functional information:**
- Gene names (e.g., *dnaA*, *rpoB*)
- Product descriptions (e.g., "DNA polymerase III")
- EC numbers (enzyme classification)
- GO terms (gene ontology)
- COG categories (functional groupings)

---

## Step 1: Run Bakta Annotation

### 1.1 Prepare Input

```bash
# Verify assembly exists and is valid
SAMPLE="F003_M_enclense"
ASSEMBLY="results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta"

# Check file
ls -lh $ASSEMBLY

# Count contigs
grep -c "^>" $ASSEMBLY

# Verify FASTA format
head -10 $ASSEMBLY
```

### 1.2 Run Bakta

```bash
# Run Bakta annotation
bakta \
    --db ~/databases/bakta/db \
    --output results/annotation/F003_M_enclense \
    --prefix F003_M_enclense \
    --threads 4 \
    --verbose \
    --keep-contig-headers \
    --compliant \
    results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta
```

**Parameter explanation:**
- `--db`: Path to Bakta database
- `--output`: Output directory
- `--prefix`: Prefix for output files
- `--threads`: Number of CPU threads
- `--verbose`: Show detailed progress
- `--keep-contig-headers`: Preserve original contig names
- `--compliant`: Generate NCBI-compliant output
- Input: Polished assembly FASTA

**⏱️ Runtime:** ~60 minutes per genome

### 1.3 Monitor Progress

Bakta runs through several stages:

```bash
# Watch progress
tail -f results/annotation/F003_M_enclense/F003_M_enclense.log
```

**Stages you'll see:**
1. `Setup` - Loading database
2. `tRNA` - Predicting tRNA genes
3. `tmRNA` - Predicting tmRNA
4. `rRNA` - Predicting rRNA (16S, 23S, 5S)
5. `ncRNA` - Other non-coding RNAs
6. `ncRNA-region` - Regulatory regions
7. `CRISPR` - CRISPR arrays
8. `CDS` - Protein-coding genes
9. `sORF` - Small ORFs
10. `gap` - Assembly gaps
11. `ori` - Replication origins
12. `Expert annotation` - Functional annotation
13. `Output` - Writing files

---

## Step 2: Explore Annotation Outputs

### 2.1 List Output Files

```bash
# Navigate to annotation directory
cd results/annotation/F003_M_enclense

# List all output files
ls -lh
```

**Output files generated:**

| File | Format | Description |
|------|--------|-------------|
| `${SAMPLE}.gbff` | GenBank | Full annotation (NCBI submission format) |
| `${SAMPLE}.gff3` | GFF3 | Genomic features (coordinates) |
| `${SAMPLE}.faa` | FASTA | Protein sequences (amino acids) |
| `${SAMPLE}.fna` | FASTA | Nucleotide sequences (CDS) |
| `${SAMPLE}.ffn` | FASTA | Gene sequences (nucleotides) |
| `${SAMPLE}.tsv` | TSV | Tab-separated summary table |
| `${SAMPLE}.txt` | Text | Human-readable summary |
| `${SAMPLE}.json` | JSON | Machine-readable metadata |
| `${SAMPLE}.log` | Log | Annotation log |
| `${SAMPLE}.png` | Image | Circular genome plot |

### 2.2 View Human-Readable Summary

```bash
# View text summary
cat F003_M_enclense.txt
```

**Example output:**
```
Genus: Microbacterium
Species: enclense

Genome Size: 3,625,480 bp
Contigs: 5
GC Content: 69.45%

Features:
  Genes: 3,456
  CDS: 3,389
  tRNA: 45
  rRNA: 9 (3x 5S, 3x 16S, 3x 23S)
  ncRNA: 8
  CRISPR: 2
```

### 2.3 View JSON Summary

```bash
# Use Python to pretty-print JSON
python3 -m json.tool F003_M_enclense.json | less
```

**JSON contains:**
- Genome statistics (size, GC, N50)
- Feature counts
- Software versions
- Database versions

---

## Step 3: Understanding GFF3 Format

### 3.1 GFF3 Structure

GFF3 (Generic Feature Format v3) stores genomic features in 9 tab-delimited columns.

```bash
# View GFF3 file
head -30 F003_M_enclense.gff3
```

**Example lines:**
```gff3
##gff-version 3
contig_1	Bakta	gene	500	1750	.	+	.	ID=GENE_0001;Name=dnaA;locus_tag=SAMPLE_00001
contig_1	Bakta	CDS	500	1750	.	+	0	ID=GENE_0001.1;Parent=GENE_0001;Name=dnaA;product=chromosomal replication initiator protein DnaA;EC_number=3.6.1.15
```

**Column breakdown:**
1. **Seqid:** `contig_1` - Contig name
2. **Source:** `Bakta` - Annotation tool
3. **Type:** `gene` or `CDS` - Feature type
4. **Start:** `500` - Start position (1-based, inclusive)
5. **End:** `1750` - End position (1-based, inclusive)
6. **Score:** `.` - Not used (quality score)
7. **Strand:** `+` or `-` - Forward or reverse strand
8. **Phase:** `0` - Reading frame (0, 1, or 2 for CDS)
9. **Attributes:** `ID=...;Name=...` - Key-value pairs

### 3.2 Extract Features from GFF3

```bash
# Count different feature types
grep -v "^#" F003_M_enclense.gff3 | cut -f3 | sort | uniq -c

# Example output:
#   3389 CDS
#   3456 gene
#     45 tRNA
#      9 rRNA
#      8 ncRNA
#      2 CRISPR
```

### 3.3 Extract Gene Information

```bash
# Get all gene names
grep "	gene	" F003_M_enclense.gff3 | grep -o "Name=[^;]*" | cut -d= -f2 | head -20

# Find specific genes
grep "dnaA" F003_M_enclense.gff3

# Extract genes on forward strand
grep "	gene	" F003_M_enclense.gff3 | awk '$7 == "+"' | wc -l

# Extract genes on reverse strand
grep "	gene	" F003_M_enclense.gff3 | awk '$7 == "-"' | wc -l
```

### 3.4 Find tRNA and rRNA

```bash
# List all tRNA genes with anticodons
grep "	tRNA	" F003_M_enclense.gff3 | grep -o "product=[^;]*"

# Count rRNA types
grep "	rRNA	" F003_M_enclense.gff3 | grep -o "product=[^;]*" | sort | uniq -c

# Expected output:
#   3 product=5S ribosomal RNA
#   3 product=16S ribosomal RNA
#   3 product=23S ribosomal RNA
```

---

## Step 4: Understanding GenBank Format

### 4.1 GenBank Structure

GenBank format combines sequence and annotation in one file.

```bash
# View GenBank file (first 100 lines)
head -100 F003_M_enclense.gbff
```

**Sections:**
1. **LOCUS** - Sequence metadata (length, molecule type, date)
2. **DEFINITION** - Brief description
3. **ACCESSION** - Database accession (empty for new genomes)
4. **FEATURES** - Annotations (genes, CDS, etc.)
5. **ORIGIN** - DNA sequence
6. **//** - End marker

### 4.2 Example GenBank Entry

```genbank
LOCUS       contig_1             3456789 bp    DNA     linear   BCT 26-JAN-2026
DEFINITION  Microbacterium enclense contig_1.
FEATURES             Location/Qualifiers
     source          1..3456789
                     /organism="Microbacterium enclense"
                     /mol_type="genomic DNA"
     gene            500..1750
                     /gene="dnaA"
                     /locus_tag="SAMPLE_00001"
     CDS             500..1750
                     /gene="dnaA"
                     /locus_tag="SAMPLE_00001"
                     /product="chromosomal replication initiator protein DnaA"
                     /protein_id="gnl|Bakta|SAMPLE_00001"
                     /translation="MSLSLWQQCLARLQDELPAQQFTTLWGKLYEVLFPFLKG..."
ORIGIN
        1 atgaaacgca ttagcatgcg ctaatgcgat cgatcgatcg atcgatcgat cgatcgatcg
       61 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg
//
```

### 4.3 Extract Information from GenBank

```bash
# Count features
grep "^     gene  " F003_M_enclense.gbff | wc -l

# Extract all gene products
grep "/product=" F003_M_enclense.gbff | cut -d'"' -f2 | head -20

# Find genes with EC numbers (enzymes)
grep "/EC_number=" F003_M_enclense.gbff | wc -l

# Extract protein sequences (translations)
grep "/translation=" F003_M_enclense.gbff | head -5
```

---

## Step 5: Analyze Protein Sequences (FAA)

### 5.1 Protein FASTA Format

```bash
# View protein sequences
head -20 F003_M_enclense.faa
```

**Example:**
```fasta
>SAMPLE_00001 chromosomal replication initiator protein DnaA
MSLSLWQQCLARLQDELPAQQFTTLWGKLYEVLFPFLKGWDNTKQYQGRPE...
>SAMPLE_00002 DNA polymerase III subunit beta
MKFTVEREHLLKPLQQVSGPLGGRPTLPILGNLLLQVADGTLSLTGTDLEA...
```

### 5.2 Protein Statistics

```bash
# Count proteins
grep -c "^>" F003_M_enclense.faa

# Calculate protein length distribution
grep -v "^>" F003_M_enclense.faa | awk '{s+=length} END {print "Total AA:", s}'

# Average protein length
grep -v "^>" F003_M_enclense.faa | \
    awk '/^>/ {if (seq) print length(seq); seq=""} !/^>/ {seq=seq $0} END {print length(seq)}' | \
    awk '{sum+=$1; n++} END {print "Average length:", sum/n, "AA"}'
```

### 5.3 Find Specific Proteins

```bash
# Search for DNA polymerase
grep -A1 "DNA polymerase" F003_M_enclense.faa

# Search for ribosomal proteins
grep -c "ribosomal protein" F003_M_enclense.faa

# Find the longest protein
grep -v "^>" F003_M_enclense.faa | \
    awk '/^>/ {if (seq) {print length(seq), header}; header=$0; seq=""} \
         !/^>/ {seq=seq $0} \
         END {print length(seq), header}' | \
    sort -rn | head -1
```

---

## Step 6: Analyze TSV Summary Table

### 6.1 TSV Structure

The TSV file is a tab-separated table with comprehensive annotation info.

```bash
# View TSV header
head -1 F003_M_enclense.tsv | tr '\t' '\n' | nl

# Example columns:
#  1	Sequence Id
#  2	Type
#  3	Start
#  4	Stop
#  5	Strand
#  6	Locus Tag
#  7	Gene
#  8	Product
#  9	DbXrefs (GO, COG, etc.)
```

### 6.2 View TSV Content

```bash
# View first 10 annotations (formatted)
head -10 F003_M_enclense.tsv | column -t

# Filter for specific feature types
grep "CDS" F003_M_enclense.tsv | head -5 | column -t
```

### 6.3 Analyze Functional Annotations

```bash
# Count hypothetical proteins
HYPOTHETICAL=$(grep -c "hypothetical protein" F003_M_enclense.tsv)
TOTAL_CDS=$(grep -c "CDS" F003_M_enclense.gff3)

echo "Hypothetical proteins: $HYPOTHETICAL / $TOTAL_CDS"
echo "Percentage: $(echo "scale=1; $HYPOTHETICAL * 100 / $TOTAL_CDS" | bc)%"

# Good annotation: <30% hypothetical
# Poor annotation: >50% hypothetical
```

### 6.4 Extract Functional Categories

```bash
# Find genes with GO terms
grep "GO:" F003_M_enclense.tsv | wc -l

# Find genes with EC numbers (enzymes)
grep "EC:" F003_M_enclense.tsv | wc -l

# Find genes with COG categories
grep "COG:" F003_M_enclense.tsv | wc -l

# Extract COG categories
grep "COG:" F003_M_enclense.tsv | grep -o "COG:[^,]*" | cut -d: -f2 | \
    cut -c1 | sort | uniq -c | sort -rn

# COG category codes:
# E: Amino acid metabolism
# G: Carbohydrate metabolism
# J: Translation, ribosomal structure
# K: Transcription
# L: Replication, recombination, repair
# etc.
```

---

## Step 7: Generate Custom Statistics

### 7.1 Create Comprehensive Summary

```bash
# Create custom statistics file
SAMPLE="F003_M_enclense"
STATS="results/annotation/${SAMPLE}/${SAMPLE}_annotation_stats.txt"

cat > $STATS <<EOF
========================================
Annotation Statistics for ${SAMPLE}
========================================

EOF

# Add genome statistics from JSON
python3 <<PYTHON_EOF >> $STATS
import json

with open('results/annotation/${SAMPLE}/${SAMPLE}.json', 'r') as f:
    data = json.load(f)

stats = data['stats']
print(f"Genome Statistics:")
print(f"  Size: {stats['size']:,} bp")
print(f"  GC content: {stats['gc']:.2f}%")
print(f"  Contigs: {stats['no_sequences']}")
print(f"  N50: {stats.get('n50', 'N/A'):,} bp")
print()
print(f"Gene Counts:")
print(f"  Total genes: {stats['no_genes']:,}")
print(f"  CDS: {stats['no_cds']:,}")
print(f"  tRNA: {stats['no_t_rna']}")
print(f"  rRNA: {stats['no_r_rna']}")
print(f"  ncRNA: {stats['no_nc_rna']}")
print(f"  CRISPR: {stats.get('no_crispr', 0)}")
PYTHON_EOF

# Add functional annotation stats
echo "" >> $STATS
echo "Functional Annotation:" >> $STATS

HYPOTHETICAL=$(grep -c "hypothetical protein" results/annotation/${SAMPLE}/${SAMPLE}.tsv || echo 0)
TOTAL_CDS=$(grep -c "CDS" results/annotation/${SAMPLE}/${SAMPLE}.gff3 || echo 1)
HYPOTHETICAL_PCT=$(echo "scale=1; $HYPOTHETICAL * 100 / $TOTAL_CDS" | bc)

echo "  Hypothetical proteins: $HYPOTHETICAL ($HYPOTHETICAL_PCT%)" >> $STATS
echo "  Functionally annotated: $((TOTAL_CDS - HYPOTHETICAL))" >> $STATS

EC_COUNT=$(grep -c "EC:" results/annotation/${SAMPLE}/${SAMPLE}.tsv || echo 0)
GO_COUNT=$(grep -c "GO:" results/annotation/${SAMPLE}/${SAMPLE}.tsv || echo 0)
COG_COUNT=$(grep -c "COG:" results/annotation/${SAMPLE}/${SAMPLE}.tsv || echo 0)

echo "  Genes with EC numbers: $EC_COUNT" >> $STATS
echo "  Genes with GO terms: $GO_COUNT" >> $STATS
echo "  Genes with COG categories: $COG_COUNT" >> $STATS

# View the summary
cat $STATS
```

### 7.2 Extract Most Common Functions

```bash
# Top 20 most common gene products
echo ""
echo "Top 20 Most Common Gene Products:"
grep "CDS" F003_M_enclense.tsv | cut -f8 | sort | uniq -c | sort -rn | head -20
```

### 7.3 Analyze COG Functional Categories

```bash
# COG category distribution
echo ""
echo "COG Functional Category Distribution:"
grep "COG:" F003_M_enclense.tsv | grep -o "COG[0-9]*" | \
    while read cog; do
        grep "$cog" ~/databases/bakta_db/cog-*.txt | cut -f3
    done | sort | uniq -c | sort -rn | head -15
```

---

## Step 8: Visualize Annotation

### 8.1 Circular Genome Plot

Bakta generates a circular genome visualization:

```bash
# View the circular plot. Only works if you have access to a GUI
eog F003_M_enclense.png
```

**Plot features:**
- **Outer ring:** Contigs
- **Second ring:** Forward strand genes (colored by COG category)
- **Third ring:** Reverse strand genes
- **Inner rings:** GC content, GC skew

### 8.2 Generate Custom Plots (Optional)

Using Python with BioPython and matplotlib:

```python
from Bio import SeqIO
import matplotlib.pyplot as plt

# Load GenBank file
records = list(SeqIO.parse("F003_M_enclense.gbff", "genbank"))

# Extract gene positions
gene_starts = []
gene_ends = []
for record in records:
    for feature in record.features:
        if feature.type == "CDS":
            gene_starts.append(int(feature.location.start))
            gene_ends.append(int(feature.location.end))

# Plot gene distribution
plt.figure(figsize=(12, 4))
plt.scatter(gene_starts, [1]*len(gene_starts), alpha=0.3, s=1)
plt.xlabel("Genome Position (bp)")
plt.title("Gene Distribution Across Genome")
plt.tight_layout()
plt.savefig("gene_distribution.png", dpi=300)
```

---

## Step 9: Compare Annotations Across Samples

### 9.1 Annotate All Samples

```bash
# Annotate remaining samples
for SAMPLE in F003_M_ginsengisoli F003_A_portus H2_M_enclense H2_M_ginsengisoli H2_A_portus; do
    echo "Annotating ${SAMPLE}..."
    
    bakta \
        --db ~/databases/bakta_db \
        --output results/annotation/${SAMPLE} \
        --prefix ${SAMPLE} \
        --threads 8 \
        --verbose \
        --keep-contig-headers \
        --compliant \
        results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta
done
```

### 9.2 Create Comparison Table

```bash
# Create comparison table
echo -e "Sample\tSize\tGenes\tCDS\ttRNA\trRNA\tHypothetical%" > annotation_comparison.tsv

for SAMPLE in F003_M_enclense F003_M_ginsengisoli F003_A_portus H2_M_enclense H2_M_ginsengisoli H2_A_portus; do
    # Extract stats from JSON
    SIZE=$(python3 -c "import json; f=open('results/annotation/${SAMPLE}/${SAMPLE}.json'); \
        print(json.load(f)['stats']['size'])")
    GENES=$(python3 -c "import json; f=open('results/annotation/${SAMPLE}/${SAMPLE}.json'); \
        print(json.load(f)['stats']['no_genes'])")
    CDS=$(python3 -c "import json; f=open('results/annotation/${SAMPLE}/${SAMPLE}.json'); \
        print(json.load(f)['stats']['no_cds'])")
    TRNA=$(python3 -c "import json; f=open('results/annotation/${SAMPLE}/${SAMPLE}.json'); \
        print(json.load(f)['stats']['no_t_rna'])")
    RRNA=$(python3 -c "import json; f=open('results/annotation/${SAMPLE}/${SAMPLE}.json'); \
        print(json.load(f)['stats']['no_r_rna'])")
    
    # Calculate hypothetical percentage
    HYPO=$(grep -c "hypothetical protein" results/annotation/${SAMPLE}/${SAMPLE}.tsv || echo 0)
    HYPO_PCT=$(echo "scale=1; $HYPO * 100 / $CDS" | bc)
    
    echo -e "${SAMPLE}\t${SIZE}\t${GENES}\t${CDS}\t${TRNA}\t${RRNA}\t${HYPO_PCT}" >> annotation_comparison.tsv
done

# View formatted table
column -t annotation_comparison.tsv
```

### 9.3 Identify Core and Accessory Genes

```bash
# Extract gene names for each sample
for SAMPLE in F003_M_enclense F003_M_ginsengisoli; do
    grep "	gene	" results/annotation/${SAMPLE}/${SAMPLE}.gff3 | \
        grep -o "Name=[^;]*" | cut -d= -f2 | sort > ${SAMPLE}_genes.txt
done

# Find common genes (core genome)
comm -12 F003_M_enclense_genes.txt F003_M_ginsengisoli_genes.txt | wc -l
echo "Core genes (shared between both Microbacterium species)"

# Find unique genes (accessory genome)
comm -23 F003_M_enclense_genes.txt F003_M_ginsengisoli_genes.txt | wc -l
echo "Unique to M. enclense"

comm -13 F003_M_enclense_genes.txt F003_M_ginsengisoli_genes.txt | wc -l
echo "Unique to M. ginsengisoli"
```

---

## Step 10: Extract Specific Gene Information

### 10.1 Find Genes of Interest

Looking for specific pathways or functions:

```bash
# Find all genes related to nitrogen metabolism
grep -i "nitrogen" F003_M_enclense.tsv | column -t

# Find antibiotic resistance genes
grep -iE "antibiotic|resistance|beta-lactam" F003_M_enclense.tsv | column -t

# Find genes related to stress response
grep -iE "stress|heat shock|cold shock" F003_M_enclense.tsv | column -t

# Find secretion systems
grep -iE "secretion|type.*secretion" F003_M_enclense.tsv | column -t
```

### 10.2 Extract Gene Sequences

```bash
# Extract specific gene by name
GENE="dnaA"

# Get nucleotide sequence
grep -A1 "$GENE" F003_M_enclense.ffn

# Get protein sequence
grep -A1 "$GENE" F003_M_enclense.faa

# Get genomic coordinates
grep "$GENE" F003_M_enclense.gff3
```

### 10.3 Create Gene Catalog

```bash
# Create a catalog of all annotated genes
cat > gene_catalog.txt <<EOF
Gene Catalog for ${SAMPLE}
============================

EOF

grep "CDS" F003_M_enclense.tsv | \
    awk -F'\t' '{print $7 "\t" $8}' | \
    sort -k1 >> gene_catalog.txt

head -20 gene_catalog.txt
```

---

## Common Issues and Troubleshooting

### Issue 1: Very High Hypothetical Percentage (>50%)

**Diagnosis:**
```bash
HYPO=$(grep -c "hypothetical protein" sample.tsv)
TOTAL=$(grep -c "CDS" sample.gff3)
echo "scale=1; $HYPO * 100 / $TOTAL" | bc
```

**Possible causes:**
1. **Novel organism** - Expected for understudied taxa
2. **Poor assembly quality** - Fragmented genes
3. **Database version** - Old Bakta database

**Solutions:**
- Update Bakta database: `bakta_db update`
- Check BUSCO completeness (should still be >90%)
- This is normal for novel organisms!

### Issue 2: Low tRNA Count (<30)

**Expected:** ~40-60 tRNA genes (one for each codon)

**Possible causes:**
- Incomplete assembly (missing regions)
- Fragmented contigs (tRNAs at contig breaks)

**Check:**
```bash
# Count tRNAs
grep "tRNA" sample.gff3 | wc -l

# List tRNA types
grep "tRNA" sample.gff3 | grep -o "product=[^;]*" | sort | uniq -c
```

### Issue 3: Unusual rRNA Count

**Expected:** 3-9 rRNA genes
- 1-3 copies of 5S
- 1-3 copies of 16S
- 1-3 copies of 23S

**Fragmented:**
```bash
# Check for fragmented rRNA
grep "rRNA" sample.gff3
# If you see many partial rRNA genes, assembly may be fragmented at rRNA operons
```

### Issue 4: Bakta Database Not Found

**Error:** `ERROR: DB not found`

**Solution:**
```bash
# Set BAKTA_DB environment variable
export BAKTA_DB=~/databases/bakta_db

# Add to ~/.bashrc for permanence
echo 'export BAKTA_DB=~/databases/bakta_db' >> ~/.bashrc
source ~/.bashrc
```

---

## Exercise: Analyze Your Annotation

### Task 1: Annotation Quality Assessment

Fill in this table for your sample:

| Metric | Value | Quality |
|--------|-------|---------|
| Total genes | | |
| CDS count | | |
| tRNA count | | (expect 40-60) |
| rRNA count | | (expect 3-9) |
| Hypothetical % | | (expect <30%) |
| Genes with EC numbers | | |
| Genes with GO terms | | |

### Task 2: Functional Analysis

1. **What are the top 5 most abundant gene functions?**
   ```bash
   grep "CDS" sample.tsv | cut -f8 | sort | uniq -c | sort -rn | head -5
   ```

2. **How many genes are involved in:**
   - Energy metabolism: `grep -i "ATP\|energy" sample.tsv | wc -l`
   - DNA replication: `grep -i "DNA.*replication\|polymerase" sample.tsv | wc -l`
   - Translation: `grep -i "ribosom\|translation" sample.tsv | wc -l`

3. **Does your genome encode CRISPR systems?**
   ```bash
   grep "CRISPR" sample.gff3
   ```

### Task 3: Comparative Analysis

If you annotated multiple samples:

1. Compare gene counts between species
2. Identify genes unique to each species
3. Calculate core genome size (shared genes)

---

## Next Steps

With completed annotations, you can now:

✅ **Week 2 analyses:**
- Screen for coral-relevant genes
- Perform comparative genomics
- Build pangenomes
- Reconstruct metabolic pathways

**Files you'll use:**
- `${SAMPLE}.gff3` - Gene coordinates
- `${SAMPLE}.faa` - Protein sequences
- `${SAMPLE}.tsv` - Functional annotations
- `${SAMPLE}.gbff` - Complete annotation

---

## Quick Reference Commands

### Essential Annotation Workflow
```bash
# 1. Run Bakta
bakta --db $BAKTA_DB --output annot_out/ --prefix sample \
    --threads 8 --verbose --compliant assembly.fasta

# 2. View summary
cat sample.txt

# 3. Extract statistics
grep -c "CDS" sample.gff3
grep -c "hypothetical" sample.tsv

# 4. Find specific genes
grep "gene_name" sample.gff3
```

### Useful One-Liners
```bash
# Count feature types
grep -v "^#" sample.gff3 | cut -f3 | sort | uniq -c

# Hypothetical percentage
echo "scale=1; $(grep -c hypothetical sample.tsv) * 100 / $(grep -c CDS sample.gff3)" | bc

# Extract gene products
grep "CDS" sample.tsv | cut -f8 | sort | uniq -c | sort -rn | head -20

# Find longest protein
bioawk -c fastx '{print $name, length($seq)}' sample.faa | sort -k2 -rn | head -1
```

---

**Continue to:** [TBD](./TUTORIAL_04_COMPARING.md)
