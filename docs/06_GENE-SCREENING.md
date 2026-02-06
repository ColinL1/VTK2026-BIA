# Tutorial 6: Gene Screening and Trait Detection

## Learning Objectives

By the end of this tutorial, you will be able to:

- Understand the difference between secondary metabolite biosynthesis and antimicrobial resistance
- Detect biosynthetic gene clusters (BGCs) using antiSMASH
- Screen for antimicrobial resistance genes with ABRicate
- Perform comprehensive AMR analysis with AMRFinderPlus
- Interpret gene screening results for ecological significance
- Generate screening summaries and reports

---

## Prerequisites

### Required Software

```bash
# Activate the gene screening conda environment
conda activate VTK2026_gene-screening

# Verify tools are installed
antismash --version
abricate --version
amrfinder --version
```

### Required Input

You must have completed:

- **Tutorial 2 (Assembly):** Polished genome assembly
- **Tutorial 3 (Annotation):** Protein sequences and GenBank files

Required files:

```bash
# Polished assembly
ls results/assembly/F003_M_enclense/F003_M_enclense_polished.fasta

# Protein FASTA (from annotation)
ls results/annotation/F003_M_enclense/F003_M_enclense/F003_M_enclense.faa

# GenBank file (optional, for antiSMASH)
ls results/annotation/F003_M_enclense/F003_M_enclense/F003_M_enclense.gbff
```

### Directory Setup

```bash
# Navigate to project directory
cd ~/VTK2026-BIA/  # change to your project name

# Create screening directories
mkdir -p results/screening/F003_M_enclense/antismash
mkdir -p results/screening/F003_M_enclense/amr
```

---

## Background: Secondary Metabolites vs. Antimicrobial Resistance

### Two Distinct Biological Functions

**Secondary Metabolites (Natural Products):**

- **Purpose:** Chemical defense, communication, competition
- **Examples:** Antibiotics, antifungals, pigments, toxins
- **Biosynthesis:** Complex multi-gene clusters (BGCs)
- **In marine bacteria:** Often novel compounds with biotech potential

**Antimicrobial Resistance (AMR):**

- **Purpose:** Protect bacteria from antimicrobial drugs
- **Examples:** Beta-lactamases, efflux pumps, target modifications
- **Genetics:** Single genes or operons

### The Tools: What Each Does

| Tool | Purpose | Targets | Output |
|------|---------|---------|--------|
| **antiSMASH** | Secondary metabolite BGC detection | NRPS, PKS, RiPPs, terpenes, etc. | BGC regions, predicted products |
| **ABRicate** | Quick AMR gene screening | AMR genes (multiple databases) | Gene hits, % identity/coverage |
| **AMRFinderPlus** | Comprehensive AMR analysis | AMR, stress, virulence genes + mutations | Detailed gene annotations |

**Key distinction:**

- ✅ **antiSMASH:** "What can this bacterium *make*?" (biosynthetic potential)
- ✅ **AMRFinderPlus:** "What can this bacterium *resist*?" (survival mechanisms)

---

## Step 1: Biosynthetic Gene Cluster Detection with antiSMASH

### 1.1 Understanding antiSMASH

**antiSMASH** (antibiotics & Secondary Metabolite Analysis SHell) is the gold standard for BGC prediction.

**What are Biosynthetic Gene Clusters?**

- Co-localized genes encoding enzymes for metabolite production
- Typical size: 10-100 kb (10-50 genes)
- Organized in operons with shared regulation

**Common BGC types:**

- **NRPS** (Non-Ribosomal Peptide Synthetases): Produce peptides without ribosomes
  - Examples: Penicillin, vancomycin, daptomycin
- **PKS** (Polyketide Synthases): Produce polyketides
  - Examples: Erythromycin, tetracycline, rapamycin
- **RiPPs** (Ribosomally synthesized and Post-translationally modified Peptides)
  - Examples: Nisin, lantibiotics, thiopeptides
- **Terpenes:** Produce terpenoid compounds
  - Examples: Carotenoids, hopanoids
- **Bacteriocins:** Antimicrobial peptides
- **Siderophores:** Iron-chelating compounds

### 1.2 How antiSMASH Works

**Algorithm:**

1. **Gene prediction:** Identifies open reading frames (ORFs)
2. **Core biosynthetic gene detection:** Searches for signature enzymes (NRPS, PKS, etc.)
3. **Cluster boundary definition:** Determines BGC start/end
4. **Annotation:** Predicts substrate specificity, tailoring enzymes
5. **Comparison:** Matches to known BGC database (MIBiG)
6. **Prediction:** Estimates chemical structure of products

**Database comparisons:**

- **MIBiG:** Minimum Information about a Biosynthetic Gene cluster
- **ClusterBlast:** Compare to known BGCs
- **KnownClusterBlast:** Compare to characterized BGCs
- **SubClusterBlast:** Compare to BGC components

### 1.3 Prepare Input Files

antiSMASH can use either:

- **FASTA file** (genome sequence) → antiSMASH does gene calling
- **GenBank file** (pre-annotated) → Uses existing annotations

```bash
# We'll use the polished assembly (FASTA)
SAMPLE="F003_M_enclense"
ASSEMBLY="results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta"

# Verify file exists
ls -lh $ASSEMBLY
```

**💡 Tip:** Using GenBank from Prokka/Bakta provides better annotations, but FASTA is sufficient.

### 1.4 Run antiSMASH

```bash
# Run antiSMASH
antismash \
    --genefinding-tool none \
    --genefinding-gff3 ${GFF} \
    --output-dir results/screening/${SAMPLE}/antismash \
    --output-basename ${SAMPLE} \
    --cpus 8 \
    --taxon bacteria \
    --cb-general \
    --cb-subclusters \
    --cb-knownclusters \
    --asf \
    --pfam2go \
    ${ASSEMBLY}
```

**Parameter explanation:**

- `--genefinding-tool prodigal`: Use Prodigal for gene prediction
  - Use this when input is FASTA
  - Omit if using GenBank file (uses existing annotations)
- `--genefinding-gff3`: path to bakta gff3 file
- `--output-dir`: Output directory for results
- `--output-basename`: Prefix for output files
- `--cpus 8`: Number of threads
- `--taxon bacteria`: Organism type (bacteria/fungi)
- `--cb-general`: Enable general ClusterBlast comparison
- `--cb-subclusters`: Enable SubClusterBlast (compare subregions)
- `--cb-knownclusters`: Enable KnownClusterBlast (known BGCs)
- `--asf`: Active Site Finder (predict substrate specificity)
- `--pfam2go`: Add Gene Ontology annotations

**Additional useful options:**

```bash
# For deeper analysis (slower):
--clusterhmmer        # Compare to BGC HMM profiles
--tigrfam             # Add TIGRFam annotations
--smcog-trees         # Generate phylogenetic trees

# For faster screening:
--minimal             # Skip most extra analyses
```

**⏱️ Runtime:** 30 minutes - 2 hours depending on genome size and analysis depth

### 1.5 Monitor antiSMASH Progress

antiSMASH provides progress updates:

```bash
# Watch the terminal output
# You'll see stages like:
# - Preparing input sequences
# - Running gene prediction
# - Detecting biosynthetic gene clusters
# - Running cluster analyses
# - Generating HTML output
```

**Example output:**

```
antiSMASH 7.0
[INFO] Preparing input sequences...
[INFO] Running gene prediction with prodigal...
[INFO] Detecting biosynthetic gene clusters...
[INFO] Found 8 candidate clusters
[INFO] Running ClusterBlast...
[INFO] Running KnownClusterBlast...
[INFO] Generating HTML visualization...
[INFO] antiSMASH complete!
```

### 1.6 Examine antiSMASH Output

```bash
# List output files
ls results/screening/${SAMPLE}/antismash/

# Key files:
# F003_M_enclense.json          - Machine-readable results
# F003_M_enclense.gbk           - GenBank with BGC annotations
# F003_M_enclense.zip           - All results packaged
# index.html                     - Interactive HTML report
# regions.js                     - BGC region data
```

### 1.7 View the Interactive HTML Report

```bash
# Open the HTML report in browser
firefox results/screening/${SAMPLE}/antismash/index.html
```

### 1.8 Navigate the antiSMASH Report

**Main sections:**

1. **Overview:**
   - Shows genome with BGC locations
   - Color-coded by BGC type
   - Click on regions to zoom in

2. **Region Details:**
   - Gene organization diagram
   - Core biosynthetic genes highlighted
   - Substrate predictions
   - Similarity to known BGCs

3. **ClusterBlast Results:**
   - Comparison to similar BGCs
   - % similarity to known clusters
   - Links to references

4. **Downloadable Data:**
   - Sequences (nucleotide, protein)
   - Annotations
   - Full reports

### 1.9 Interpret BGC Results

**Example BGC report:**

```
Region 1.1 (1-45 kb): NRPS
  Core biosynthetic genes: 3
  Similarity: 65% to biosynthetic gene cluster for bacitracin
  Predicted core structure: Unknown
```

**What to look for:**

✅ **High-confidence BGC:**
>
- >70% similarity to known cluster → Likely produces similar compound
- Complete core biosynthetic genes
- Presence of resistance genes (self-protection)

⚠️ **Novel/Orphan BGC:**

- <40% similarity to known clusters → Potentially novel compound!
- Unusual gene organization
- **Exciting for natural product discovery!**

**BGC abundance expectations:**

| Genome Size | Typical # BGCs | Note |
|-------------|----------------|------|
| 3-4 Mb (like *Microbacterium*) | 3-8 | Smaller genomes, fewer clusters |
| 4-6 Mb (like *Alteromonas*, *Marinobacter*) | 6-15 | Medium genomes |
| >8 Mb (like *Streptomyces*) | 20-40 | Specialized producers |

### 1.10 Extract BGC Summary from JSON

```python
import json
import os

sample = "F003_M_enclense"
json_file = f"results/screening/{sample}/antismash/{sample}.json"

if os.path.exists(json_file):
    with open(json_file) as f:
        data = json.load(f)
    
    records = data.get('records', [])
    total_clusters = sum(len(r.get('areas', [])) for r in records)
    
    print(f"Total BGCs detected: {total_clusters}\n")
    
    # Count by type
    cluster_types = {}
    for record in records:
        for area in record.get('areas', []):
            for product in area.get('products', []):
                cluster_types[product] = cluster_types.get(product, 0) + 1
    
    if cluster_types:
        print("BGC types:")
        for ctype, count in sorted(cluster_types.items()):
            print(f"  {ctype}: {count}")
else:
    print(f"antiSMASH JSON not found: {json_file}")
```

**Example output:**

```
Total BGCs detected: 7

BGC types:
  NRPS: 2
  T1PKS: 1
  arylpolyene: 1
  siderophore: 2
  terpene: 1
```

### 1.11 Export BGC Sequences

```bash
# antiSMASH includes BGC sequences in GenBank format
# Extract specific BGC region (example: region 1)

# List all regions
ls results/screening/${SAMPLE}/antismash/*.region*.gbk

# You can convert to FASTA if needed:
from Bio import SeqIO

for record in SeqIO.parse("results/screening/${SAMPLE}/antismash/${SAMPLE}.gbk", "genbank"):
    # Extract regions marked as BGC
    print(f">{record.id}")
    print(record.seq)
```

---

## Step 2: Quick AMR Screening with ABRicate

### 2.1 Understanding ABRicate

**ABRicate** performs rapid screening for antimicrobial resistance genes using curated databases.

**How it works:**

1. BLAST search against AMR gene databases
2. Filter by identity (%) and coverage (%)
3. Report matching genes with metadata

**Advantages:**

- ⚡ Fast (minutes, not hours)
- 🗄️ Multiple curated databases
- 📊 Simple tab-delimited output
- 🔍 Good for initial screening

**Limitations:**

- Nucleotide-level only (misses divergent homologs)
- No point mutation detection
- Basic annotations

### 2.2 ABRicate Databases

ABRicate includes several databases:

| Database | Source | Focus | # Genes |
|----------|--------|-------|---------|
| **CARD** | Comprehensive Antibiotic Resistance Database | Broad AMR coverage | ~6,000 |
| **NCBI** | National Center for Biotechnology Information | Curated bacterial AMR | ~5,000 |
| **ResFinder** | DTU, Denmark | Clinical resistance | ~3,000 |
| **ARG-ANNOT** | Antibiotic Resistance Gene-ANNOTation | French curated | ~2,000 |
| **MEGARes** | Veterinary medicine focus | Livestock AMR | ~7,000 |

**Best practice:** Screen against multiple databases for comprehensive coverage.

### 2.3 Update ABRicate Databases

```bash
# Update all databases (one-time setup, or periodic refresh)
abricate --setupdb

# List available databases
abricate --list
```

**Expected output:**

```text
DATABASE    SEQUENCES  DBTYPE  DATE
card        6253       nucl    2024-12-15
ncbi        5386       nucl    2024-11-20
resfinder   3077       nucl    2024-10-10
...
```

### 2.4 Run ABRicate Screening

We'll screen against three major databases:

```bash
SAMPLE="F003_M_enclense"
ASSEMBLY="results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta"
AMR_DIR="results/screening/${SAMPLE}/amr"

# Screen against CARD
echo "Screening against CARD database..."
abricate \
    --db card \
    --minid 75 \
    --mincov 50 \
    ${ASSEMBLY} \
    > ${AMR_DIR}/${SAMPLE}_card.tab

# Screen against NCBI
echo "Screening against NCBI database..."
abricate \
    --db ncbi \
    --minid 75 \
    --mincov 50 \
    ${ASSEMBLY} \
    > ${AMR_DIR}/${SAMPLE}_ncbi.tab

# Screen against ResFinder
echo "Screening against ResFinder database..."
abricate \
    --db resfinder \
    --minid 75 \
    --mincov 50 \
    ${ASSEMBLY} \
    > ${AMR_DIR}/${SAMPLE}_resfinder.tab
```

**Parameter explanation:**

- `--db`: Database to search
- `--minid 75`: Minimum % identity (75% is permissive)
  - Use 90% for high-confidence hits
  - Use 75-80% for divergent homologs
- `--mincov 50`: Minimum % coverage of reference gene
  - Use 80% for complete genes
  - Use 50% for partial/fragmented genes

**Thresholds:**

- **High confidence:** ≥90% identity, ≥80% coverage
- **Moderate confidence:** 80-90% identity, 50-80% coverage
- **Low confidence:** <80% identity or <50% coverage (may be false positive)

**⏱️ Runtime:** 1-5 minutes per database

### 2.5 View ABRicate Results

```bash
# View CARD results
column -t ${AMR_DIR}/${SAMPLE}_card.tab

# Count hits
echo "CARD hits: $(tail -n +2 ${AMR_DIR}/${SAMPLE}_card.tab | wc -l)"
echo "NCBI hits: $(tail -n +2 ${AMR_DIR}/${SAMPLE}_ncbi.tab | wc -l)"
echo "ResFinder hits: $(tail -n +2 ${AMR_DIR}/${SAMPLE}_resfinder.tab | wc -l)"
```

**Example output:**

```
FILE                                  SEQUENCE  START    END      GENE        COVERAGE    %COVERAGE  %IDENTITY  DATABASE  ACCESSION
F003_M_enclense_polished.fasta        contig_1  234567   235890   tetA        1-1323/1323  100.00    95.23      card      ARO:3000162
F003_M_enclense_polished.fasta        contig_1  456789   457234   aadA        1-445/445    100.00    87.64      card      ARO:3000001
```

**Columns explained:**

- `SEQUENCE`: Contig where gene was found
- `START/END`: Genomic coordinates
- `GENE`: Resistance gene name
- `COVERAGE`: Alignment coverage (bp aligned / total gene length)
- `%COVERAGE`: Percentage of reference gene covered
- `%IDENTITY`: Nucleotide identity to reference
- `ACCESSION`: Database accession number

### 2.6 Combine ABRicate Results

```bash
# Merge all database results
echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" \
    > ${AMR_DIR}/${SAMPLE}_amr_all.tab

cat ${AMR_DIR}/${SAMPLE}_*.tab | grep -v "^#" | grep -v "^FILE" \
    >> ${AMR_DIR}/${SAMPLE}_amr_all.tab

# View combined results
column -t ${AMR_DIR}/${SAMPLE}_amr_all.tab | less -S
```

### 2.7 Interpret ABRicate Results

**Common resistance mechanisms found:**

| Gene/Class | Resistance To | Mechanism | Concern Level |
|------------|---------------|-----------|---------------|
| **tetA/B/C** | Tetracycline | Efflux pump | Moderate |
| **aadA** | Aminoglycosides | Adenylation | Moderate |
| **bla*** | Beta-lactams | Hydrolysis | High |
| **erm*** | Macrolides | rRNA methylation | Moderate |
| **vanA/B** | Vancomycin | Target modification | Critical |
| **qnr*** | Fluoroquinolones | Target protection | High |

**Expected in marine bacteria:**

- **Low numbers (0-5 genes):** Typical environmental bacteria
- **Moderate (5-10):** May have intrinsic resistance or acquired genes
- **High (>10):** Investigate contamination or multi-resistant strain

**Intrinsic vs. Acquired resistance:**

- **Intrinsic:** Species naturally resistant (e.g., *Lactobacillus* to vancomycin)
- **Acquired:** Resistance genes on plasmids/mobile elements (concerning)

**💡 Tip:** Cross-reference plasmid detection (from assembly) with AMR genes to identify mobile resistance.

---

## Step 3: Comprehensive AMR Analysis with AMRFinderPlus

### 3.1 Understanding AMRFinderPlus

**AMRFinderPlus** is NCBI's reference tool for AMR gene detection.

**Advantages over ABRicate:**

- Protein-level search (more sensitive than nucleotide BLAST)
- Detects point mutations conferring resistance
- Curated annotations with mechanism information
- Organism-aware (considers intrinsic resistance)
- Linked to extensive NCBI metadata

**What AMRFinderPlus detects:**

1. **AMR genes:** Beta-lactamases, efflux pumps, etc.
2. **Resistance mutations:** Point mutations in gyrase, ribosome, etc.
3. **Virulence factors:** Toxins, adhesins
4. **Stress response genes:** Metal resistance, biocide tolerance

### 3.2 How AMRFinderPlus Works

**Algorithm:**

1. **HMM search:** Searches protein sequences against curated HMMs
2. **BLAST search:** Nucleotide and protein BLAST against reference database
3. **Point mutation detection:** Identifies known resistance SNPs
4. **Hierarchy resolution:** Prioritizes best-matching genes
5. **Annotation:** Adds resistance class, mechanism, subclass

**Database:**

- Curated by NCBI
- Updated regularly (quarterly releases)
- Organism-specific intrinsic resistance markers

### 3.3 Update AMRFinderPlus Database

```bash
# Update AMRFinderPlus database (do this periodically)
amrfinder --update

# Check database version
amrfinder --database_version
```

**Expected output:**

```
Database version: 2024-11-18.1
```

### 3.4 Run AMRFinderPlus

AMRFinderPlus uses **protein sequences** for optimal sensitivity:

```bash
SAMPLE="F003_M_enclense"
PROTEIN_FASTA="results/annotation/${SAMPLE}/${SAMPLE}/${SAMPLE}.faa"
AMR_DIR="results/screening/${SAMPLE}/amr"

# Run AMRFinderPlus
amrfinder \
    --plus \
    --protein ${PROTEIN_FASTA} \
    --threads 8 \
    --output ${AMR_DIR}/${SAMPLE}_amrfinder.tsv
```

**Parameter explanation:**

- `--plus`: Include point mutations and additional virulence/stress genes
  - Omit for AMR genes only (faster)
- `--protein`: Input protein FASTA (from Prokka/Bakta)
- `--threads`: Number of CPU threads
- `--output`: Output TSV file

**Optional parameters:**

```bash
# Specify organism for intrinsic resistance filtering
--organism Escherichia        # Skip intrinsic E. coli resistance

# Use nucleotide input (less sensitive)
--nucleotide assembly.fasta

# Set thresholds
--ident_min 0.9               # Minimum identity (default: varies by gene)
--coverage_min 0.5            # Minimum coverage (default: varies by gene)

# Output formats
--plus --annotation_format prokka  # If using Prokka annotations
```

**⏱️ Runtime:** 5-20 minutes

### 3.5 View AMRFinderPlus Results

```bash
# View results (tab-delimited)
column -t ${AMR_DIR}/${SAMPLE}_amrfinder.tsv | less -S

# Count total hits
tail -n +2 ${AMR_DIR}/${SAMPLE}_amrfinder.tsv | wc -l
```

**Example output:**

```
Protein identifier  Gene symbol  Sequence name  Scope      Element type  Element subtype  Class                   Subclass          Method      Target length  Reference length  % Coverage  % Identity  Alignment length  Accession  Name                                           HMM id  HMM description
PROKKA_00123        tetA         contig_1       core       AMR           AMR              TETRACYCLINE            TETRACYCLINE      BLASTP      1200           1200              100.00      95.3        1200              WP_000123456  tetracycline resistance protein TetA          ...     ...
PROKKA_00234        aadA1        contig_1       core       AMR           AMR              AMINOGLYCOSIDE          STREPTOMYCIN      HMM         264            264               100.00      88.6        264               WP_000234567  aminoglycoside nucleotidyltransferase AadA1   ...     ...
```

**Key columns:**

- `Protein identifier`: Locus tag from annotation
- `Gene symbol`: Standard gene name (e.g., tetA, blaKPC)
- `Element type`: AMR / VIRULENCE / STRESS
- `Class`: Antibiotic/compound class (e.g., TETRACYCLINE, BETA-LACTAM)
- `Subclass`: Specific drug (e.g., CEFTAZIDIME, METHICILLIN)
- `Method`: BLASTP / BLASTX / HMM / PARTIAL / POINT_MUTATION
- `% Coverage / % Identity`: Alignment quality metrics

### 3.6 Filter High-Confidence Hits

```bash
# Extract high-confidence AMR genes (≥90% identity, ≥80% coverage)
awk -F'\t' 'NR==1 || ($12>=80 && $13>=90)' \
    ${AMR_DIR}/${SAMPLE}_amrfinder.tsv \
    > ${AMR_DIR}/${SAMPLE}_amrfinder_high_confidence.tsv

# Count high-confidence hits
echo "High-confidence AMR genes:"
tail -n +2 ${AMR_DIR}/${SAMPLE}_amrfinder_high_confidence.tsv | wc -l
```

### 3.7 Summarize by Resistance Class

```bash
# Count genes by antibiotic class
echo "AMR genes by class:"
tail -n +2 ${AMR_DIR}/${SAMPLE}_amrfinder.tsv | \
    cut -f7 | sort | uniq -c | sort -rn
```

**Example output:**

```
AMR genes by class:
  3 TETRACYCLINE
  2 AMINOGLYCOSIDE
  1 BETA-LACTAM
  1 MACROLIDE
```

### 3.8 Identify Point Mutations

```bash
# Extract point mutations (if using --plus)
grep "POINT_MUTATION" ${AMR_DIR}/${SAMPLE}_amrfinder.tsv | \
    column -t
```

**Example:**

```
Protein identifier  Gene symbol  Mutation       Class              Method          
PROKKA_01234        gyrA         S83L           FLUOROQUINOLONE    POINT_MUTATION
PROKKA_02345        rpoB         S531L          RIFAMPIN           POINT_MUTATION
```

**Point mutations are critical:**

- Cannot be removed (chromosomal)
- Often confer high-level resistance
- Important for tracking resistance evolution

### 3.9 Compare ABRicate vs. AMRFinderPlus

```bash
# Extract gene names from both tools
echo "ABRicate genes:"
tail -n +2 ${AMR_DIR}/${SAMPLE}_amr_all.tab | cut -f5 | sort -u

echo ""
echo "AMRFinderPlus genes:"
tail -n +2 ${AMR_DIR}/${SAMPLE}_amrfinder.tsv | cut -f2 | sort -u
```

**Why might results differ?**

- **AMRFinderPlus detects more:** Protein-level search is more sensitive
- **ABRicate may find extra:** Uses multiple databases with different curation
- **Best practice:** Use AMRFinderPlus as primary, ABRicate as supplementary

---

## Step 4: Generate Screening Summary

### 4.1 Create Comprehensive Summary

```bash
SAMPLE="F003_M_enclense"
SUMMARY_FILE="results/screening/${SAMPLE}/${SAMPLE}_screening_summary.txt"

cat > ${SUMMARY_FILE} <<EOF
========================================
Gene Screening Summary for ${SAMPLE}
========================================
Date: $(date +"%Y-%m-%d")

EOF

# Add antiSMASH summary
echo "========================================" >> ${SUMMARY_FILE}
echo "SECONDARY METABOLITE BIOSYNTHESIS (antiSMASH)" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}

if [ -f "results/screening/${SAMPLE}/antismash/${SAMPLE}.json" ]; then
    python3 - << PYTHON_EOF >> ${SUMMARY_FILE}
import json
try:
    with open("results/screening/${SAMPLE}/antismash/${SAMPLE}.json", 'r') as f:
        data = json.load(f)
    records = data.get('records', [])
    total_clusters = sum(len(r.get('areas', [])) for r in records)
    print(f"Total BGCs detected: {total_clusters}\n")
    
    # Count by type
    cluster_types = {}
    for record in records:
        for area in record.get('areas', []):
            for product in area.get('products', []):
                cluster_types[product] = cluster_types.get(product, 0) + 1
    
    if cluster_types:
        print("BGC types:")
        for ctype, count in sorted(cluster_types.items()):
            print(f"  - {ctype}: {count}")
    else:
        print("No BGC types detected")
except Exception as e:
    print(f"Error parsing antiSMASH results: {e}")
PYTHON_EOF
else
    echo "antiSMASH results not found" >> ${SUMMARY_FILE}
fi

echo "" >> ${SUMMARY_FILE}

# Add AMR summary
echo "========================================" >> ${SUMMARY_FILE}
echo "ANTIMICROBIAL RESISTANCE (AMRFinderPlus)" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}

if [ -f "${AMR_DIR}/${SAMPLE}_amrfinder.tsv" ]; then
    AMR_COUNT=$(tail -n +2 "${AMR_DIR}/${SAMPLE}_amrfinder.tsv" | wc -l)
    echo "Total AMR/stress/virulence genes: ${AMR_COUNT}" >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
    
    # Summarize by class
    echo "By resistance class:" >> ${SUMMARY_FILE}
    tail -n +2 "${AMR_DIR}/${SAMPLE}_amrfinder.tsv" | \
        cut -f7 | sort | uniq -c | sort -rn | \
        awk '{print "  - " $2 ": " $1}' >> ${SUMMARY_FILE}
    
    echo "" >> ${SUMMARY_FILE}
    
    # Point mutations
    MUTATION_COUNT=$(grep -c "POINT_MUTATION" "${AMR_DIR}/${SAMPLE}_amrfinder.tsv" || echo "0")
    echo "Point mutations detected: ${MUTATION_COUNT}" >> ${SUMMARY_FILE}
    
    if [ "$MUTATION_COUNT" -gt 0 ]; then
        echo "Mutations:" >> ${SUMMARY_FILE}
        grep "POINT_MUTATION" "${AMR_DIR}/${SAMPLE}_amrfinder.tsv" | \
            awk -F'\t' '{print "  - " $2 " (" $3 "): " $7}' >> ${SUMMARY_FILE}
    fi
else
    echo "AMRFinderPlus results not found" >> ${SUMMARY_FILE}
fi

echo "" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}
echo "OUTPUT FILES" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}
echo "antiSMASH report: results/screening/${SAMPLE}/antismash/index.html" >> ${SUMMARY_FILE}
echo "AMRFinderPlus: ${AMR_DIR}/${SAMPLE}_amrfinder.tsv" >> ${SUMMARY_FILE}
echo "ABRicate (combined): ${AMR_DIR}/${SAMPLE}_amr_all.tab" >> ${SUMMARY_FILE}

echo "✓ Summary generated: ${SUMMARY_FILE}"
```

### 4.2 Display Summary

```bash
# Display the summary
cat ${SUMMARY_FILE}
```

**Example output:**

```
========================================
Gene Screening Summary for F003_M_enclense
========================================
Date: 2026-02-06

========================================
SECONDARY METABOLITE BIOSYNTHESIS (antiSMASH)
========================================
Total BGCs detected: 7

BGC types:
  - NRPS: 2
  - siderophore: 2
  - T1PKS: 1
  - arylpolyene: 1
  - terpene: 1

========================================
ANTIMICROBIAL RESISTANCE (AMRFinderPlus)
========================================
Total AMR/stress/virulence genes: 8

By resistance class:
  - TETRACYCLINE: 3
  - AMINOGLYCOSIDE: 2
  - BETA-LACTAM: 1
  - MACROLIDE: 1
  - STRESS: 1

Point mutations detected: 0

========================================
OUTPUT FILES
========================================
antiSMASH report: results/screening/F003_M_enclense/antismash/index.html
AMRFinderPlus: results/screening/F003_M_enclense/amr/F003_M_enclense_amrfinder.tsv
ABRicate (combined): results/screening/F003_M_enclense/amr/F003_M_enclense_amr_all.tab
```

---

## Step 5: Biological Interpretation

### 5.1 Cross-Referencing BGC and AMR

**Look for connections:**

1. **Self-resistance genes near BGCs:**

   ```bash
   # Check if AMR genes are near BGCs (same contig, close coordinates)
   # Example: Aminoglycoside resistance near aminoglycoside-producing BGC
   ```

2. **Mobile genetic elements:**

   ```bash
   # If AMR genes on plasmids → likely acquired
   # If BGCs on plasmids → horizontal transfer potential
   ```

3. **Ecological interpretation:**
   - BGCs + low AMR → Natural product producer, pristine environment
   - BGCs + high AMR → Possible pollution-impacted site
   - Low BGCs + high AMR → Concerning for AMR spread

### 5.2 Practical Applications

**For natural product discovery:**

1. Prioritize novel BGCs (<40% similarity)
2. Focus on NRPS/PKS for antibiotic potential
3. Consider marine-specific modifications (halogenation, etc.)

**For AMR surveillance:**

1. Track mobile resistance genes (on plasmids)
2. Monitor critical resistance (carbapenems, colistin)
3. Compare across spatial/temporal samples

**For coral health:**

1. Beneficial BGCs: Siderophores, antimicrobials (protect coral)
2. Concerning: Virulence factors in opportunistic pathogens

---

## Step 6: Visualize Screening Results (Optional)

### 6.1 Generate BGC Distribution Plot

```bash
# Create simple bar plot of BGC types
python3 - << 'EOF'
import json
import matplotlib.pyplot as plt

sample = "F003_M_enclense"
json_file = f"results/screening/{sample}/antismash/{sample}.json"

with open(json_file) as f:
    data = json.load(f)

cluster_types = {}
for record in data.get('records', []):
    for area in record.get('areas', []):
        for product in area.get('products', []):
            cluster_types[product] = cluster_types.get(product, 0) + 1

if cluster_types:
    plt.figure(figsize=(10, 6))
    plt.bar(cluster_types.keys(), cluster_types.values(), color='steelblue')
    plt.xlabel('BGC Type')
    plt.ylabel('Count')
    plt.title(f'Biosynthetic Gene Clusters in {sample}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'results/screening/{sample}/bgc_distribution.png', dpi=300)
    print(f"Plot saved: results/screening/{sample}/bgc_distribution.png")
EOF
```

---

## Advanced Topics

### Using the Pipeline Script

For automated screening of multiple samples:

```bash
# Run the full pipeline script
bash scripts/04_gene_screening.sh F003_M_enclense 8
```

The script performs all steps:

1. antiSMASH BGC detection
2. ABRicate screening (CARD, NCBI, ResFinder)
3. AMRFinderPlus comprehensive analysis
4. Summary generation

### Screening Multiple Samples

```bash
# Loop through all samples
for SAMPLE in F003_M_enclense H2_Alteromonas_portus F003_Marinobacter_adhaerens; do
    echo "Processing ${SAMPLE}..."
    bash scripts/04_gene_screening.sh ${SAMPLE} 8
done

# Combine results across samples
python3 scripts/helpers/combine_screening_results.py results/screening/
```

---

## Troubleshooting

### antiSMASH Issues

**Problem:** antiSMASH fails with "No genes found"

```bash
# Solution: Check input format
head results/assembly/${SAMPLE}/${SAMPLE}_polished.fasta

# Ensure FASTA headers are simple (no special characters)
# Re-run with --genefinding-tool prodigal
```

**Problem:** antiSMASH runs out of memory

```bash
# Solution: Run with --minimal flag for faster, lighter analysis
antismash --minimal --output-dir results/screening/${SAMPLE}/antismash ${ASSEMBLY}
```

**Problem:** ClusterBlast fails

```bash
# Solution: Disable database comparisons if needed
antismash --skip-clusterblast --skip-knownclusters ${ASSEMBLY}
```

### AMRFinderPlus Issues

**Problem:** "Database not found"

```bash
# Solution: Update database
amrfinder --update --database /path/to/amrfinder/db
```

**Problem:** Few or no hits detected

```bash
# Solution 1: Check input file format
head ${PROTEIN_FASTA}

# Ensure proper FASTA format (headers start with >)

# Solution 2: Lower thresholds (use with caution)
amrfinder --protein ${PROTEIN_FASTA} --ident_min 0.7 --coverage_min 0.4
```

**Problem:** High number of hits (>20)

```bash
# Investigate potential contamination
# Run CheckM2 to verify single organism
checkm2 predict --input ${ASSEMBLY} --output-directory results/qc/checkm2
```

### ABRicate Issues

**Problem:** Database update fails

```bash
# Solution: Manually download databases
cd $(which abricate | xargs dirname)/../db
git clone https://github.com/tseemann/abricate-db.git
```

**Problem:** No hits in any database

```bash
# This is actually normal for some environmental bacteria!
# Verify with AMRFinderPlus (more sensitive)
```

---

## Summary

You have completed comprehensive gene screening for:

1. **Secondary metabolite biosynthesis** (antiSMASH)
   - Detected BGCs for natural product discovery
   - Identified biosynthetic potential
   - Compared to known compound databases

2. **Antimicrobial resistance** (ABRicate, AMRFinderPlus)
   - Screened for AMR genes
   - Detected point mutations
   - Assessed resistance profiles

3. **Ecological interpretation**
   - Distinguished BGCs from AMR
   - Evaluated environmental context
   - Identified potential applications

### Key Takeaways

| Feature | antiSMASH (BGCs) | AMRFinderPlus (AMR) |
|---------|------------------|---------------------|
| **Purpose** | Natural product discovery | Resistance surveillance |
| **Target** | Multi-gene biosynthetic clusters | Single resistance genes + mutations |
| **Output** | Predicted chemical structures | Resistance phenotypes |
| **Ecological significance** | Competitive advantage, signaling | Survival under antibiotic pressure |
| **Application** | Drug discovery, biotechnology | Public health, ecology |

---

## Additional Resources

### antiSMASH

- Documentation: <https://docs.antismash.secondarymetabolites.org/>
- Publication: Blin et al. (2023) Nucleic Acids Research
- MIBiG database: <https://mibig.secondarymetabolites.org/>

### AMRFinderPlus

- Documentation: <https://github.com/ncbi/amr/wiki>
- Publication: Feldgarden et al. (2021) Scientific Reports
- NCBI Pathogen Detection: <https://www.ncbi.nlm.nih.gov/pathogens/>

### ABRicate

- GitHub: <https://github.com/tseemann/abricate>
- CARD: <https://card.mcmaster.ca/>
- ResFinder: <https://cge.food.dtu.dk/services/ResFinder/>

### Marine Natural Products

- MarinLit database: <http://pubs.rsc.org/marinlit/>
- Review: Carroll et al. (2021) "Marine natural products"
- Coral microbiome: <https://coraltraits.org/>

### BiG-SCAPE: BGC Clustering and Family Analysis

- GitHub: <https://github.com/medema/BiG-SCAPE>
- Publication: Navarro-Muñoz et al. (2020) Nature Chemistry Biology
- Tutorial: <https://bigscape-corason.secondarymetabolites.org/>

---
