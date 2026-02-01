# Introduction to Sequencing File Formats

## Overview

This guide will help you understand the file formats you'll encounter when working with sequencing data. Understanding these formats is crucial for quality control, assembly, and downstream analysis.

## Table of Contents
1. [FASTA Format](#fasta-format)
2. [FASTQ Format](#fastq-format)
3. [GFF and GenBank Formats](#gff-and-genbank-formats)
4. [Quality Scores](#quality-scores)
5. [ONT-Specific Considerations](#ont-specific-considerations)
6. [Common File Operations](#common-file-operations)
7. [Practical Tips](#practical-tips)

---

## FASTA Format

### What is FASTA?

FASTA is the simplest and most universal sequence format in bioinformatics. It stores nucleotide or protein sequences in plain text format without quality information.

### Structure

```
>sequence_identifier description (optional)
ATGCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGA
```

**Components:**
- **Header line**: Starts with `>` followed by a unique identifier
- **Sequence lines**: The actual nucleotide or amino acid sequence (can span multiple lines)

### Example

```fasta
>contig_1 length=1500 coverage=45.2
ATGAAACGCATTAGCATGCGCTAATGCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>contig_2 length=2300 coverage=38.7
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

### When You'll Use FASTA

- **Assembled genomes**: Final output from assemblers
- **Reference sequences**: When comparing your assembly to known genomes
- **Gene sequences**: Annotated coding sequences (CDS)
- **Protein sequences**: Translated amino acid sequences

---

## FASTQ Format

### What is FASTQ?

FASTQ is the standard format for storing raw sequencing reads along with their quality scores. This is typically the first file format you'll encounter from the sequencer.

### Structure

Each sequence entry consists of **four lines**:

```
@sequence_identifier description
ATGCGATCGATCGATCGATCGATCG
+
IIIIIIHHHGGGFFFDDDBBBAAAA
```

1. **Line 1**: Header line starting with `@` (identifier and optional description)
2. **Line 2**: The raw nucleotide sequence
3. **Line 3**: Separator line starting with `+` (sometimes repeats the identifier, usually empty)
4. **Line 4**: Quality scores (one character per base, same length as line 2)

### Example

```fastq
@read_001 runid=a1b2c3d4 sampleid=bacteria_01
ATGAAACGCATTAGCATGCGCTAATGCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIHHHHHGGGGFFFFEEEDDDCCCCBBBBAAA@@@???>>>===<<<;;;
@read_002 runid=a1b2c3d4 sampleid=bacteria_01
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIHHHHHGGGGGFFFFEEEEDDDCCCCBBBBAAAA@
```

### Important Notes

- **One read = 4 lines**: Always! This is critical for parsing
- **Case sensitive**: Both sequence and quality strings
- **No line breaks**: Unlike FASTA, FASTQ sequences cannot span multiple lines
- **Matching lengths**: Sequence (line 2) and quality (line 4) must be identical length

---

## GFF and GenBank Formats

### What are Annotation Formats?

After assembling your bacterial genome (FASTA), the next step is **annotation** - identifying genes, their locations, and their functions. This information is stored in specialized formats like GFF and GenBank.

### GFF Format (General Feature Format)

**GFF3** is a tab-delimited text format that stores genomic features (genes, CDS, tRNA, etc.) and their coordinates.

#### Structure

GFF files have **9 columns** separated by tabs:

```
seqid  source  type  start  end  score  strand  phase  attributes
```

**Column descriptions:**
1. **seqid**: Sequence/contig identifier (matches FASTA header)
2. **source**: Program that generated the annotation (e.g., Prokka, PGAP)
3. **type**: Feature type (gene, CDS, tRNA, rRNA, etc.)
4. **start**: Start position (1-based, inclusive)
5. **end**: End position (1-based, inclusive)
6. **score**: Confidence score (often "." for not applicable)
7. **strand**: `+` (forward) or `-` (reverse)
8. **phase**: For CDS, indicates reading frame (0, 1, or 2)
9. **attributes**: Semicolon-separated key-value pairs (ID, Name, product, etc.)

#### Example

```gff3
##gff-version 3
contig_1	Prokka	gene	500	1750	.	+	.	ID=gene001;Name=dnaA;product=chromosomal replication initiator protein DnaA
contig_1	Prokka	CDS	500	1750	.	+	0	ID=cds001;Parent=gene001;product=chromosomal replication initiator protein DnaA;EC_number=3.6.1.15
contig_1	Prokka	gene	2100	2450	.	-	.	ID=gene002;Name=rpsT;product=30S ribosomal protein S20
contig_1	Prokka	CDS	2100	2450	.	-	0	ID=cds002;Parent=gene002;product=30S ribosomal protein S20
contig_1	Prokka	tRNA	3000	3075	.	+	.	ID=trna001;product=tRNA-Ala;anticodon=GCC
```

#### Key Features

- **Hierarchical**: Genes have child features (CDS, exons)
- **Attributes**: Flexible key-value system for metadata
- **Coordinates**: 1-based (position 1 is the first base)
- **Human-readable**: Easy to parse and view in text editors

#### Common Feature Types in Bacterial Genomes

| Type | Description | Example |
|------|-------------|----------|
| `gene` | Complete gene | Entire gene region |
| `CDS` | Coding sequence | Protein-coding region |
| `tRNA` | Transfer RNA | tRNA genes |
| `rRNA` | Ribosomal RNA | 16S, 23S, 5S rRNA |
| `ncRNA` | Non-coding RNA | Regulatory RNAs |
| `repeat_region` | Repetitive elements | CRISPR arrays |

#### Viewing GFF Files

```bash
# View header and first few features
head -n 20 annotation.gff

# Count genes
grep -c "\tgene\t" annotation.gff

# Count CDS (protein-coding genes)
grep -c "\tCDS\t" annotation.gff

# Extract all tRNA features
grep "\ttRNA\t" annotation.gff

# Get gene names
grep "\tgene\t" annotation.gff | grep -o "Name=[^;]*" | cut -d= -f2
```

### GenBank Format

**GenBank** is a rich text format that combines sequence data with annotations. It's the standard format for NCBI submissions and downloads.

#### Structure

GenBank files have two main sections:
1. **Header**: Metadata about the sequence
2. **Features**: Annotations (similar to GFF)
3. **Sequence**: The actual DNA sequence (at the end)

#### Example

```genbank
LOCUS       contig_1                5000 bp    DNA     linear   BCT 26-JAN-2026
DEFINITION  Bacterial sp. contig_1, complete sequence.
ACCESSION   
VERSION     
KEYWORDS    .
SOURCE      Bacterial sp.
  ORGANISM  Bacterial sp.
            Bacteria.
FEATURES             Location/Qualifiers
     source          1..5000
                     /organism="Bacterial sp."
                     /mol_type="genomic DNA"
     gene            500..1750
                     /gene="dnaA"
                     /locus_tag="GENE001"
     CDS             500..1750
                     /gene="dnaA"
                     /locus_tag="GENE001"
                     /product="chromosomal replication initiator protein DnaA"
                     /protein_id="gnl|Prokka|GENE001"
                     /translation="MSLSLWQQCLARLQDELP..."
     gene            complement(2100..2450)
                     /gene="rpsT"
                     /locus_tag="GENE002"
     CDS             complement(2100..2450)
                     /gene="rpsT"
                     /locus_tag="GENE002"
                     /product="30S ribosomal protein S20"
                     /translation="MAKKTKRLRNKGSRRRLSQKK..."
     tRNA            3000..3075
                     /product="tRNA-Ala"
                     /anticodon=(pos:3034..3036,aa:Ala,seq:gcc)
ORIGIN
        1 atgaaacgca ttagcatgcg ctaatgcgat cgatcgatcg atcgatcgat cgatcgatcg
       61 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg
      121 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg
//
```

#### Key Features

- **Rich metadata**: Organism info, references, comments
- **Translations**: Includes amino acid sequences for proteins
- **Standardized**: Required format for NCBI GenBank submissions
- **Integrated**: Sequence and annotation in one file
- **Location syntax**: Special notation for features

#### GenBank Location Syntax

```
500..1750              Simple range (forward strand)
complement(2100..2450) Reverse strand
join(100..200,300..400) Feature spans multiple regions (for eukaryotes)
<1..500                Feature extends before sequence start
1000..>5000            Feature extends beyond sequence end
```

#### File Extensions

- `.gbk`, `.gb`, `.genbank` - GenBank format
- `.gbff` - GenBank flat file
- `.gff`, `.gff3` - GFF format

### Converting Between Formats

You'll often need to convert between formats:

```bash
# GFF to GenBank (using seqret from EMBOSS)
seqret -sequence genome.fasta -feature -fformat gff -osformat genbank -outseq output.gbk

# GenBank to FASTA (extract sequences)
seqret -sequence annotation.gbk -osformat fasta -outseq sequences.fasta

# GenBank to protein FASTA
seqret -sequence annotation.gbk -osformat fasta -feature -fformat cds -outseq proteins.faa

# Using BioPython (Python)
python -c "from Bio import SeqIO; SeqIO.convert('file.gbk', 'genbank', 'file.fasta', 'fasta')"
```

### Quick Statistics from Annotation Files

```bash
# Count features in GFF
grep -v "^#" annotation.gff | cut -f3 | sort | uniq -c

# Example output:
#   1523 CDS
#   1550 gene
#     45 tRNA
#      3 rRNA
#     12 ncRNA

# Extract gene products (what do the genes do?)
grep "\tCDS\t" annotation.gff | grep -o "product=[^;]*" | cut -d= -f2 | sort | uniq -c | sort -rn | head

# Get genome statistics from GenBank
grep "^LOCUS" annotation.gbk
```


---

## Quality Scores

### Phred Quality Scores

Quality scores represent the probability that a base call is incorrect. They use ASCII characters to encode numeric values.

**Formula:**
```
Q = -10 × log₁₀(P)
```

Where:
- `Q` = Phred quality score
- `P` = Probability of incorrect base call

### Common Quality Scores

| Phred Score | ASCII Character | Error Probability | Accuracy |
|-------------|----------------|-------------------|----------|
| 10 | `+` | 1 in 10 | 90% |
| 20 | `5` | 1 in 100 | 99% |
| 30 | `?` | 1 in 1,000 | 99.9% |
| 40 | `I` | 1 in 10,000 | 99.99% |

### Encoding

Most modern sequencing platforms use **Phred+33** encoding:
- ASCII character value - 33 = Quality score
- For example: `I` (ASCII 73) → 73 - 33 = Quality score 40

### Reading Quality Scores

```
Sequence: ATGCGATCG
Quality:  IIIHHHGGG

I = 40 (excellent)
H = 39 (excellent)
G = 38 (very good)
```

**Rule of thumb:**
- **Q ≥ 30**: High quality, trustworthy base
- **Q 20-29**: Good quality, generally acceptable
- **Q 10-19**: Low quality, consider filtering
- **Q < 10**: Very low quality, should be filtered

---

## Common File Operations

### Viewing Files

**Never open large FASTQ files directly!** They can be gigabytes in size.

```bash
# View first 4 lines (1 read) of FASTQ
head -n 4 file.fastq

# View first 8 lines (2 reads) of FASTQ
head -n 8 file.fastq

# View first 2 sequences of FASTA
head -n 4 file.fasta

# Count number of reads in FASTQ
grep -c "^@" file.fastq  # Caution: may overcount if @ appears in quality

# Better way to count reads
wc -l file.fastq | awk '{print $1/4}'

# View compressed files
zcat file.fastq.gz | head -n 4
```

### Checking File Integrity

```bash
# Count sequences in FASTA
grep -c "^>" file.fasta

# Verify FASTQ structure (should be divisible by 4)
wc -l file.fastq

# Check if file is gzipped
file file.fastq.gz
```

### Basic Statistics

```bash
# Get sequence lengths from FASTA
grep -v "^>" file.fasta | awk '{print length}'

# Calculate N50 (use specialized tools)
# seqkit stats file.fasta
# or
# assembly-stats file.fasta
```

### Useful Tools

- **seqkit**: Swiss army knife for FASTA/FASTQ manipulation
- **seqtk**: Toolkit for processing sequences
- **NanoPlot**: ONT-specific quality control
- **FastQC**: General sequencing quality control
- **filtlong**: Quality filtering for long reads

---

## Practical Tips

### 1. Always Compress Your Files

```bash
# Compress
gzip file.fastq  # Creates file.fastq.gz

# Decompress
gunzip file.fastq.gz  # Creates file.fastq

# Work with compressed files directly (preferred)
zcat file.fastq.gz | head
```

**Compression ratio**: FASTQ files compress ~70-80% (10 GB → 2-3 GB)

### 2. Quality Control First

Before attempting assembly:
- Check read length distribution
- Examine quality score distribution
- Look for adapter contamination
- Filter low-quality reads

```bash
# Example with NanoPlot
NanoPlot -t 4 --fastq reads.fastq.gz -o qc_report/
```

### 3. Filtering Strategies

For ONT bacterial assembly:
- **Minimum length**: 1,000 bp (filters very short reads)
- **Minimum quality**: Q7-Q10 (filters poorest reads)
- **Target coverage**: 50-100× for good assembly

```bash
# Example with filtlong
filtlong --min_length 1000 --min_mean_q 7 reads.fastq.gz > filtered.fastq
```

### 4. Understanding Coverage

**Coverage** = (Total bases sequenced) / (Genome size)

For a 5 Mb bacterial genome:
- 250 Mb of data = 50× coverage
- 500 Mb of data = 100× coverage

**Rule of thumb**: 50-100× coverage for bacterial genomes with ONT

### 5. Common Pitfalls

**Don't:**

- Open huge FASTQ files in text editors (will crash)
- Forget that FASTQ is 4 lines per read
- Ignore quality scores during assembly
- Mix up FASTQ (`@`) and FASTA (`>`) headers
- Delete raw data before confirming assembly quality

**Do:**

- Always work with compressed files when possible
- Use command-line tools for inspection
- Keep raw data backed up
- Document your quality filtering parameters
- Check file formats before running tools

### 6. Sanity Checks

Before running assembly:

```bash
# 1. File exists and is not empty
ls -lh reads.fastq.gz

# 2. File is properly gzipped
zcat reads.fastq.gz | head -n 4

# 3. FASTQ structure is valid (lines divisible by 4)
zcat reads.fastq.gz | wc -l  # Should be divisible by 4

# 4. Sequences and quality scores match in length
zcat reads.fastq.gz | head -n 4 | awk 'NR==2 {s=length} NR==4 {q=length} END {print s==q}'
# Should print "1" for true
```

---

## Quick Reference

### File Format Comparison

| Feature | FASTA | FASTQ |
|---------|-------|-------|
| Purpose | Sequences only | Sequences + quality |
| Header | `>` | `@` |
| Lines per entry | 2+ (header + sequence) | 4 (fixed) |
| Quality info | No | Yes |
| Typical use | Assemblies, references | Raw reads |
| File size | Smaller | Larger |

### Common File Extensions

- `.fasta`, `.fa`, `.fna` - FASTA nucleotide
- `.faa` - FASTA amino acid
- `.fastq`, `.fq` - FASTQ
- `.fastq.gz`, `.fq.gz` - Compressed FASTQ
- `.fast5` - ONT raw signal (binary)
- `.bam` - Aligned reads (binary)

### Workflow Overview

```
ONT Sequencer
    ↓
Pod5 (raw signals)
    ↓ [basecalling]
FASTQ (raw reads)
    ↓ [quality control]
FASTQ (filtered reads)
    ↓ [assembly]
FASTA (contigs/genome)
    ↓ [annotation]
GFF/GenBank (annotated genome)
```

---

## Additional Resources

- **ONT Community**: https://community.nanoporetech.com/
- **seqkit manual**: https://bioinf.shenwei.me/seqkit/
- **Understanding FASTQ**: https://en.wikipedia.org/wiki/FASTQ_format
- **Phred quality scores**: https://en.wikipedia.org/wiki/Phred_quality_score

---

## ONT-Specific Considerations
### ONT FASTQ Headers

ONT FASTQ files contain rich metadata in headers:

```fastq
@read_id runid=abc123 read=12345 ch=45 start_time=2026-01-26T10:30:00Z
```

**Common fields:**
- `runid`: Unique flow cell identifier
- `read`: Read number
- `ch`: Channel number (MinION has 512 channels)
- `start_time`: When the read was sequenced
- `barcode`: Sample barcode (if multiplexed)

### Basecalling

ONT sequencers produce **raw signal data** (FAST5 files) that must be "basecalled" to convert electrical signals to DNA sequences (FASTQ).

**Tools:**
- **Guppy**: ONT's official basecaller
- **Dorado**: Newer, faster basecaller
- **Bonito**: Research basecaller

You'll typically receive FASTQ files that have already been basecalled.

---
