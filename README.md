# Nottingham MRes: LSDV Genomic Analysis Pipeline

Complete bioinformatic workflows for sequencing quality assessment, de novo assembly, reference construction, phylogenetic analysis, and recombination detection of Lumpy Skin Disease Virus (LSDV) genomes.

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
  - [Environment Setup](#environment-setup)
  - [Creating Conda/Micromamba Environments](#creating-condamicromamba-environments)
- [Scripts and Workflows](#scripts-and-workflows)
  - [QC Scripts](#qc-scripts)
  - [Assembly Scripts](#assembly-scripts)
  - [Analysis Scripts](#analysis-scripts)
  - [Figure Scripts](#figure-scripts)
- [Usage Guide](#usage-guide)
- [Input/Output Specifications](#inputoutput-specifications)
- [Troubleshooting](#troubleshooting)
- [Data Availability](#data-availability)

---

## Project Overview

This repository contains custom bioinformatic scripts developed for comprehensive analysis of LSDV whole-genome sequences, including:

- **Quality Control & Filtering**: Assessment and trimming of raw sequencing reads
- **De Novo Assembly**: Genome assembly from filtered reads across multiple sequencing workflows
- **Reference-Guided Consensus Generation**: Creation of high-quality consensus sequences using hybrid references
- **Phylogenetic Analysis**: Maximum likelihood tree reconstruction and phylogenetic inference
- **Recombination Detection**: Identification of recombination signals in viral genomes
- **Data Visualization**: Publication-quality figures for sequencing efficiency and genome coverage

---

## Repository Structure

```
Nottingham-MRes/
├── qc_scripts/
│   ├── remove_rawdedupPlus.slurm
│   ├── raw_FiltPlus.slurm
│   ├── PyCount_Rds2.slurm
│   └── forBBMAP_KX.slurm
├── assembly_scripts/
│   ├── Con_Gen_05_Pur4.slurm
│   └── Spades_Iso.slurm
├── analysis_scripts/
│   ├── GenB_search.py
│   ├── run_raxml_all_and_support.sh
│   ├── simple_snp_heatmap_i.py
│   └── 3seq_cor.sh
├── figure_scripts/
│   ├── Enr_cov_contig_plot_working.slurm
│   ├── nested_bar_plot.py
│   └── viral_enrichment_panel2_4.py
├── environments/
│   ├── environment_core_analysis.yml
│   ├── environment_python_analysis.yml
│   ├── environment_read_mapping.yml
│   ├── environment_recombination.yml
│   └── environment_visualization.yml
├── supplementary_data/
│   └── [Reference genomes, supporting data]
└── README.md
```

---

## Installation

### Environment Setup

This project uses **5 separate conda/micromamba environments** to manage dependencies across different analysis stages. Multiple environments are used to avoid package conflicts and maintain compatibility with specific tool versions.

### Creating Conda/Micromamba Environments

#### Prerequisites

Ensure you have conda or micromamba installed:

```bash
# For conda
conda --version

# For micromamba (recommended for faster installation)
micromamba --version
```

#### Installing Environments

Create each environment from its corresponding `.yml` file:

```bash
# Environment 1: Core QC & Assembly
conda env create -f environments/environment_core_analysis.yml

# Environment 2: Python Visualization & Analysis
conda env create -f environments/environment_python_analysis.yml

# Environment 3: Read Mapping & Assembly Validation
conda env create -f environments/environment_read_mapping.yml

# Environment 4: Recombination & Advanced Phylogenetics
micromamba create -f environments/environment_recombination.yml  # Created with micromamba

# Environment 5: Advanced Figure Generation
conda env create -f environments/environment_visualization.yml
```

#### Activating Environments

```bash
# Using conda
conda activate environment_core_analysis
conda activate environment_python_analysis
conda activate environment_read_mapping
conda activate environment_visualization

# Using micromamba
micromamba activate environment_recombination
```

#### Environment Descriptions

| Environment | Python | Key Tools | Purpose |
|---|---|---|---|
| **environment_core_analysis** | 3.10 | FastQC, Fastp, SPAdes, QUAST, RAxML-NG, ModelTest-NG, BLAST, BWA, Bedtools, MultiQC | QC, assembly, basic phylogenetics |
| **environment_python_analysis** | 3.10 | BioPython, Pandas, NumPy, Matplotlib, PyQt6 | Python analysis, visualization |
| **environment_read_mapping** | 3.7 | BBMap, BWA, Samtools, Minimap2, RagTag, Picard, Mummer | Read mapping, assembly validation |
| **environment_recombination** | 3.11 | 3Seq, HyPhy, OpenMPI | Recombination detection, advanced phylogenetics |
| **environment_visualization** | 3.13 | BioPython, Pandas, NumPy, Matplotlib, PySide6, Pillow | Advanced figure generation |

---

## Scripts and Workflows

### QC Scripts

#### 1. remove_rawdedupPlus.slurm
**Purpose**: Quality control filtering of paired-end reads WITH deduplication

**Environment**: `environment_core_analysis`

**Dependencies**: Fastp v0.24.0

**Input Files**:
- Paired-end raw FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)

**Output Files**:
- Filtered FASTQ files: `Filtered_fastq2/*_R1_filtered.fastq.gz`, `*_R2_filtered.fastq.gz`
- JSON stats: `Filtered_fastq2/*_fastp_stats.json`
- HTML reports: `Filtered_fastq2/*_fastp_report.html`

**Usage**:
```bash
conda activate environment_core_analysis

# Modify INPUT_DIR path in script
INPUT_DIR="/path/to/raw/fastq/files"
mkdir -p Filtered_fastq2

sbatch remove_rawdedupPlus.slurm
```

**Parameters Explained**:
- `--dedup`: Remove PCR duplicates
- `--correction`: Perform base correction
- `--trim_poly_g`: Trim polyG tails
- `--length_required 110`: Keep reads ≥110 bp
- `--thread 10`: Use 10 threads

---

#### 2. raw_FiltPlus.slurm
**Purpose**: Quality control filtering of paired-end reads WITHOUT deduplication

**Environment**: `environment_core_analysis`

**Dependencies**: Fastp v0.24.0

**Input Files**:
- Paired-end raw FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)

**Output Files**:
- Filtered FASTQ files: `Filtered_fastq/*_R1_filtered.fastq.gz`, `*_R2_filtered.fastq.gz`
- JSON stats: `Filtered_fastq/*_fastp_stats.json`
- HTML reports: `Filtered_fastq/*_fastp_report.html`

**Usage**:
```bash
conda activate environment_core_analysis

# Modify INPUT_DIR path in script
INPUT_DIR="/path/to/raw/fastq/files"
mkdir -p Filtered_fastq

sbatch raw_FiltPlus.slurm
```

**Differences from remove_rawdedupPlus.slurm**:
- No `--dedup` flag (deduplication not applied)
- Useful for comparing deduplication effects on downstream analysis

---

#### 3. PyCount_Rds2.slurm
**Purpose**: Count total reads in filtered FASTQ files

**Environment**: `environment_python_analysis`

**Dependencies**: Python 3.10, gzip module

**Input Files**:
- Filtered paired-end FASTQ files (`*_R1_filtered.fastq.gz`, `*_R2_filtered.fastq.gz`)

**Output Files**:
- `read_counts.txt`: Tab-delimited file with sample names and total reads

**Usage**:
```bash
conda activate environment_python_analysis

# Run from directory containing filtered FASTQ files
sbatch PyCount_Rds2.slurm

# View results
cat read_counts.txt
```

**Output Format**:
```
Base Name          Total Reads
sample_05          15,234,562
sample_07          12,456,789
```

---

#### 4. forBBMAP_KX.slurm
**Purpose**: Read binning to reference LSDV genome (KX894508.1)

**Environment**: `environment_read_mapping`

**Dependencies**: BBMap v39.06

**Input Files**:
- Filtered paired-end FASTQ files (`*_R1_filtered.fastq.gz`, `*_R2_filtered.fastq.gz`)
- Reference FASTA: LSDV complete genome (KX894508.1)

**Output Files**:
- Mapped reads: `mappedKX/*_mapped.fastq`
- Unmapped reads: `unmappedKX/*_unmapped.fastq`

**Usage**:
```bash
conda activate environment_read_mapping

# Modify reference path in script
REFERENCE="/path/to/reference/KX894508.fasta"

# Run from directory containing filtered FASTQ files
sbatch forBBMAP_KX.slurm
```

**Key Parameters**:
- `nodisk`: Keeps reference index in memory (requires ~8GB RAM)
- Default sensitivity and specificity settings

---

### Assembly Scripts

#### 5. Con_Gen_05_Pur4.slurm
**Purpose**: Generate consensus sequences from BAM alignment files with IUPAC ambiguity codes

**Environment**: `environment_core_analysis`

**Dependencies**: Bedtools v2.31.1, bcftools v1.10.2

**Input Files**:
- BAM alignment file (sorted): `alignment.sorted.bam`
- Reference FASTA: `reference.fasta`
- Reference index: `reference.fasta.fai`

**Output Files**:
- Regions BED file: `regions.bed`
- VCF files (variants, filtered, normalized): `variants.vcf.gz`, `filtered_variants.vcf.gz`, `normalized_variants.vcf.gz`
- Final consensus: `consensus.fa`

**Usage**:
```bash
conda activate environment_core_analysis

# Modify paths in script
REFERENCE="/path/to/reference.fasta"
INPUT_BAM="/path/to/alignment.sorted.bam"
OUTPUT_DIR="/path/to/output/directory"

sbatch Con_Gen_05_Pur4.slurm
```

**Key Parameters**:
- `MIN_COVERAGE=100`: Minimum coverage for consensus region
- `MIN_ALLELE_FREQ=0.5`: Minimum allele frequency (50%)
- `MIN_QUAL=20`: Minimum base quality threshold
- `-I` flag in bcftools consensus: Enables IUPAC ambiguity encoding (heterozygous sites preserved)

**Output Interpretation**:
- Positions with coverage ≥100 are included
- Heterozygous sites (allele freq 50-100%) are encoded with IUPAC codes (e.g., R for A/G, Y for C/T)
- Homozygous alternate sites use standard nucleotide codes

---

#### 6. Spades_Iso.slurm
**Purpose**: De novo genome assembly using SPAdes isolate mode

**Environment**: `environment_core_analysis`

**Dependencies**: SPAdes v4.0.0

**Input Files**:
- Interleaved paired-end FASTQ files: `sample.fastq.gz`

**Output Files**:
- Assembly directory: `Spades_Iso/sample/`
  - Main contig file: `contigs.fasta`
  - Scaffolds: `scaffolds.fasta`
  - Assembly graph: `assembly_graph.fastg`
  - SPAdes log: `spades.log`

**Usage**:
```bash
conda activate environment_core_analysis

# Modify input directory in script
INPUT_DIR="/path/to/interleaved/fastq/files"
OUTPUT_DIR="/path/to/output/directory"

sbatch Spades_Iso.slurm
```

**Key Parameters**:
- `--isolate`: Use isolate mode (assumes high coverage, low complexity)
- `-k 21,33,55,77,99,127`: k-mer sizes for progressive assembly
- `--threads 10`: Use 10 threads
- `--memory 250`: Allocate 250GB RAM

**Notes**:
- Requires interleaved FASTQ format (use `interleave.py` from BBTools if needed)
- K-mer range suitable for viral genomes (~150kb)

---

### Analysis Scripts

#### 7. GenB_search.py
**Purpose**: Automated retrieval of LSDV genome sequences from NCBI GenBank

**Environment**: `environment_python_analysis`

**Dependencies**: BioPython, internet connection

**Input Files**:
- Text file with GenBank accession numbers: `Genbank_forTree.txt` (one accession per line)

**Output Files**:
- Individual FASTA files: `Genbank_Genomes/{accession}.fasta`

**Usage**:
```bash
conda activate environment_python_analysis

# Create input file
echo "NC_003027.1" > Genbank_forTree.txt
echo "KX894508.1" >> Genbank_forTree.txt

# Run script
python GenB_search.py
```

**Input File Format**:
```
NC_003027.1
KX894508.1
MG191649.1
```

**Output Structure**:
```
Genbank_Genomes/
├── NC_003027.1.fasta
├── KX894508.1.fasta
└── MG191649.1.fasta
```

**Notes**:
- Requires valid email address in script (used by NCBI)
- Respects NCBI rate limits automatically
- Creates output directory if it doesn't exist

---

#### 8. run_raxml_all_and_support.sh
**Purpose**: Build maximum likelihood phylogenetic trees with bootstrap support

**Environment**: `environment_core_analysis` (note: references `raxml-ng-env`, ensure RAxML-NG is installed in environment_core_analysis)

**Dependencies**: RAxML-NG v1.2.2, ModelTest-NG v0.1.7

**Input Files**:
- Aligned FASTA sequences: `alignment.fasta`

**Output Files**:
- Best tree: `LSDVFULL_GTRG4.raxml.bestTree`
- Supported tree: `LSDVFULL_GTRG4_support.raxml.support`
- Bootstrap trees: `LSDVFULL_GTRG4.raxml.bootstraps`
- Optimization logs: `LSDVFULL_GTRG4.raxml.log`

**Usage**:
```bash
conda activate environment_core_analysis

# Modify input file in script
MSA_FILE="/path/to/alignment.fasta"
OUTPUT_PREFIX="analysis_prefix"

sbatch run_raxml_all_and_support.sh
```

**Key Parameters**:
- `--model GTR+G+F`: GTR substitution model with gamma rate heterogeneity
- `--threads 10`: Use 10 threads
- Automatic bootstrap replicates (typically 100)

**Output Interpretation**:
- `.raxml.support`: Tree with bootstrap support values at nodes (0-100)
- `.raxml.bootstraps`: Individual bootstrap trees (for further analysis)
- `.raxml.bestTree`: Topology with highest likelihood

---

#### 9. simple_snp_heatmap_i.py
**Purpose**: Calculate pairwise SNP distances and generate publication-quality heatmap

**Environment**: `environment_python_analysis`

**Dependencies**: BioPython, NumPy, Pandas, Matplotlib, Seaborn

**Input Files**:
- Aligned FASTA file: `alignment.fasta`

**Output Files**:
- Heatmap PNG: `{output_prefix}.png`
- Distance matrix CSV: `{output_prefix}_distances.csv`

**Usage**:
```bash
conda activate environment_python_analysis

python simple_snp_heatmap_i.py alignment.fasta output_prefix
```

**Example**:
```bash
python simple_snp_heatmap_i.py combined_genomes.fasta figure4_snp_distances
# Outputs: figure4_snp_distances.png and figure4_snp_distances_distances.csv
```

**Output Files**:
- PNG heatmap with annotated SNP counts
- CSV distance matrix for further analysis

**Features**:
- Counts gaps and differences (important for indel analysis)
- Automatic sample name cleaning for readability
- Summary statistics printed to console
- Red-Yellow-Blue color scheme for publication quality

---

#### 10. 3seq_cor.sh
**Purpose**: Detect recombination signals using 3Seq algorithm

**Environment**: `environment_recombination` (micromamba)

**Dependencies**: 3Seq v1.8.0, OpenMPI

**Input Files**:
- Aligned FASTA sequences: `alignment.fasta`

**Output Files**:
- Main results: `{prefix}.3s.rec` (recombinant triplets)
- Log file: `{prefix}.3s.log`

**Usage**:
```bash
# Using micromamba
micromamba activate environment_recombination

# Modify input file and output prefix in script
sbatch 3seq_cor.sh
```

**Parameters in Script**:
- `-d`: Use only distinct sequences (removes duplicates)
- `-t0.05`: P-value threshold (0.05 = 5%)
- Automatic 1000×1000×1000 p-value lookup table

**Output Format** (3s.rec file):
```
Sequence1  Sequence2  Sequence3  P-value  Signal Type
seq_05     seq_07     ref_KX     0.001    Recombination signal
```

**Interpretation**:
- Each line represents a triplet of sequences
- P-values < 0.05 indicate significant recombination signal
- Dunn-Sidak correction applied for multiple testing

---

### Figure Scripts

#### 11. Enr_cov_contig_plot_working.slurm
**Purpose**: Generate coverage and contig alignment visualization plots

**Environment**: `environment_visualization`

**Dependencies**: Pandas, NumPy, Matplotlib, LineCollection

**Input Files**:
- Coverage data TSV: `coverage_data.tsv` (columns: sample, start, coverage)
- Contig alignment TSV: `contigs.tsv` (columns: sample, start, end, length)

**Output Files**:
- PDF plot: `contig_coverage_plot.pdf`

**Usage**:
```bash
conda activate biofigures

# Modify file paths in script
COVERAGE_DIR="/path/to/coverage/data"
CONTIG_DIR="/path/to/contig/data"
OUTPUT_DIR="/path/to/output"

sbatch Enr_cov_contig_plot_working.slurm
```

**Input TSV Format** (coverage):
```
sample	start	coverage
05	1	150
05	100	145
07	1	200
```

**Input TSV Format** (contigs):
```
sample	start	end	length
05	1000	5000	4000
05	6000	8500	2500
07	1500	4000	2500
```

**Output**:
- High-resolution PDF with coverage tracks and contig alignments
- Publication-ready figure with ITR regions highlighted

---

#### 12. nested_bar_plot.py
**Purpose**: Generate sequencing efficiency comparison across library preparation strategies

**Environment**: `environment_python_analysis`

**Dependencies**: Matplotlib, Argparse

**Input Files**:
- Tab-delimited data file with sequencing metrics

**Output Files**:
- PNG figure: `sequencing_efficiency.png` (or specified output)

**Usage**:
```bash
conda activate environment_python_analysis

python nested_bar_plot.py input_data.tsv --output sequencing_efficiency.png
```

**Input Format**:
```
pipeline	raw	filtered	mapped	unique
CL-SG	1000000	850000	1000	900
CL-CAP	1000000	900000	900000	850000
CC-LYS	1000000	920000	30100	28000
CC-PUR	1000000	940000	522400	490000
```

**Output**:
- Nested horizontal bars showing read retention at each processing stage
- Log-scale x-axis for visualization of highly disparate values
- Percentage labels for unique (viral-enriched) reads

---

#### 13. viral_enrichment_panel2_4.py
**Purpose**: Visualize viral enrichment efficiency across sequencing workflows

**Environment**: `environment_python_analysis`

**Dependencies**: Matplotlib, NumPy, Argparse

**Input Files**:
- None (hardcoded values for LSDV workflows)

**Output Files**:
- PNG figure: `viral_enrichment_panel2_4.png` (or specified output)

**Usage**:
```bash
conda activate environment_python_analysis

python viral_enrichment_panel2_4.py --output viral_enrichment_panel2_4.png
```

**Output**:
- Vertical bar plot comparing viral enrichment across four workflows
- Log-scale y-axis to visualize wide range of enrichment values
- Color-coded by workflow type (Clinical vs. Culture)

---

## Usage Guide

### Typical Workflow

#### Step 1: Quality Control & Read Binning
```bash
conda activate environment_core_analysis

# Run fastp filtering (choose with or without deduplication)
sbatch remove_rawdedupPlus.slurm

# Map reads to LSDV reference
conda activate environment_read_mapping
sbatch forBBMAP_KX.slurm

# Count reads at each stage
conda activate environment_python_analysis
sbatch PyCount_Rds2.slurm
```

#### Step 2: De Novo Assembly
```bash
conda activate environment_core_analysis

# Prepare interleaved FASTQ files from mapped reads
sbatch Spades_Iso.slurm
```

#### Step 3: Reference Construction & Consensus
```bash
conda activate environment_core_analysis

# Generate consensus sequence from BAM alignment to hybrid reference
sbatch Con_Gen_05_Pur4.slurm
```

#### Step 4: Phylogenetic Analysis
```bash
conda activate environment_core_analysis

# Retrieve reference sequences
conda activate environment_python_analysis
python GenB_search.py

# Build phylogenetic trees with bootstrap support
conda activate environment_core_analysis
sbatch run_raxml_all_and_support.sh

# Calculate pairwise distances
conda activate environment_python_analysis
python simple_snp_heatmap_i.py alignment.fasta output_prefix
```

#### Step 5: Recombination Detection
```bash
micromamba activate environment_recombination
sbatch 3seq_cor.sh
```

#### Step 6: Figure Generation
```bash
# Coverage and contig plots
conda activate environment_visualization
sbatch Enr_cov_contig_plot_working.slurm

# Efficiency and enrichment plots
conda activate environment_python_analysis
python nested_bar_plot.py input_data.tsv --output sequencing_efficiency.png
python viral_enrichment_panel2_4.py --output viral_enrichment_panel2_4.png
```

---

## Input/Output Specifications

### File Format Requirements

#### FASTQ Files
- **Format**: Standard Illumina FASTQ format (gzip compressed recommended)
- **Naming**: `samplename_R1.fastq.gz` and `samplename_R2.fastq.gz` for paired-end
- **Quality Encoding**: Phred+33 (modern Illumina default)

#### FASTA Files
- **Format**: Standard FASTA format
- **Sequence IDs**: Use informative names without spaces or special characters
- **Line Length**: Typically 80 characters per line (not enforced)

#### Alignment Files
- **BAM**: Binary alignment format (indexed with `.bai` files when required)
- **Must be sorted**: `samtools sort -o output.sorted.bam input.bam`

#### Reference Genomes
- **Format**: FASTA
- **Index Required**: `.fai` index file (generate with `samtools faidx reference.fasta`)

### Directory Structure Expectations

```
project_directory/
├── raw_fastq/              # Input raw sequencing reads
├── filtered_fastq/         # Output from fastp QC filtering
├── mapped_reads/           # Output from BBMap binning
├── assemblies/             # Output from SPAdes
├── alignments/             # BAM files
├── consensus_sequences/    # Output from consensus generation
├── phylogenetic_analysis/  # RAxML outputs
└── figures/                # Generated plots and visualizations
```

---

## Troubleshooting

### Common Issues

#### 1. Conda/Micromamba Environment Not Found
```bash
# Problem: "conda: command not found"

# Solution: Initialize conda
/path/to/miniconda/bin/conda init bash
source ~/.bashrc

# For micromamba
eval "$(micromamba shell hook --shell bash)"
```

#### 2. Out of Memory (OOM) Errors
- **SPAdes assembly**: Increase `--memory` parameter (default 250GB)
- **BBMap reference indexing**: Use `nodisk` option to keep index in memory
- **Solution**: Monitor available memory with `free -h` and adjust SLURM `--mem` parameter

#### 3. NCBI Connection Timeout
```bash
# Problem: GenB_search.py times out downloading sequences

# Solution: Check internet connection and NCBI status
# Retry with rate limiting or download sequences manually

# Manual alternative:
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=KX894508.1&rettype=fasta" \
  -O KX894508.1.fasta
```

#### 4. Missing Reference Index
```bash
# Problem: "reference.fasta.fai not found"

# Solution: Generate index
samtools faidx /path/to/reference.fasta
```

#### 5. Mismatched Read Counts (R1 vs R2)
- **Problem**: PyCount_Rds2.py reports different counts for paired reads
- **Cause**: Reads discarded during filtering (adapter trimming, length filtering)
- **Solution**: This is normal; verify both reads remain paired with BBTools `repair.sh`

#### 6. SPAdes Assembly Fails
```bash
# Problem: "spades.py: command not found"

# Solution: Verify SPAdes installation
conda activate environment_core_analysis
spades.py --version

# If not installed, update environment
conda env update -f environments/environment_base.yml
```

#### 7. RAxML-NG Tree Building Slow
- **Cause**: Large number of sequences (>100) or very long alignments
- **Solution**: Increase `--threads` parameter in run_raxml_all_and_support.sh
- **Alternative**: Reduce bootstrap replicates or use faster model

#### 8. 3Seq No Results
```bash
# Problem: 3s.rec file empty or contains no recombination signals

# Possible causes:
# - Very similar sequences (low variation)
# - Alignment issues (poorly aligned regions)
# - High p-value threshold (try -t0.1 for less stringent)

# Solution: Check alignment quality with alignment viewer
# Verify sequences aren't too similar for recombination detection
```

---

## Data Availability

**Raw Sequencing Data**

Raw sequencing reads for the 2022 Mongolian LSDV outbreak samples will be deposited in the NCBI Sequence Read Archive (SRA). 

**Genome Sequences**

All five LSDV complete genome sequences analyzed in this study are now available in this repository's `/genomes` folder:

- Three new 2022 Mongolian outbreak assemblies (FASTA and GenBank annotation formats)
- Two reassembled 2021 Mongolian samples (FASTA format)

**Sequence Alignments**

Sequence alignments used for all phylogenetic and comparative genomics analyses are available in the `/alignments` folder:

- **Phylogenetic analysis**: `Fig_3.6A_FULL.fasta` (full-length genome alignment), `Fig_3.6B_GPCR.fasta` (GPCR region alignment), `Fig_3.7_Regional.fasta` (regional alignment)
- **SNP distance analysis**: `Fig_3.8_HeatMap.fasta` (alignment for pairwise distance calculations)
- **Recombination analysis**: `3seq_CORE.fasta` (core genome region), `3seq_5P_Accessory.fasta` (5' accessory region), `3seq_5Pext_Accessory.fasta` (extended 5' region), `3seq_3P_Accessory.fasta` (3' accessory region)
- **Phylogenetic reconstruction of the core genome**: `Fig_3.9A_CORE_BP32.fasta` (central recombinant segment) and `Fig_3.9B_CORE_BP32CAT.fasta` (concatenated flanking breakpoint regions)
- **Phylogenetic reconstruction of the 3' accessory genome**: `Fig_3.10_3P_BP8.fasta` (3' accessory recombinant segment)

These alignments enable independent verification of all analyses presented in Results sections 3.3.1–3.3.3.

**Submission Status**

These sequences are in the process of being formally submitted to NCBI GenBank and will be assigned accession numbers upon acceptance.

**Accession Numbers**

GenBank accession numbers and SRA project identifiers will be updated in this README upon public release of data.


---

## Citation

If you use these scripts in your research, please cite:

```
Wright, C., et al. (2026). LSDV Genomic Analysis Pipeline. 
GitHub repository: https://github.com/CarolineWright/Nottingham-MRes
```

---

## License

[Insert]

---

## Contact

For questions or issues regarding these scripts, please contact: caroline.wright@pirbright.ac.uk

---

## Acknowledgments

Scripts developed as part of MRes research at the University of Nottingham.

### Generative AI Disclosure

This project utilized Claude AI (Anthropic's Claude.ai platform) in the following capacities:

- **README Generation**: The comprehensive README.md documentation was collaboratively generated through iterative refinement with Claude, ensuring accuracy, clarity, and completeness in describing all 13 bioinformatic scripts and 5 conda environments.

- **Script Formatting & Standardization**: All 13 bioinformatic scripts were formatted and standardized using Claude, with specific HPC paths replaced with generalized placeholders (`/path/to/...`) to enhance portability and usability across different systems.

- **Documentation Structure**: Claude assisted in organizing complex methodological workflows into clear, hierarchical documentation with consistent formatting.

**Important Note**: All scripts were reviewed and validated by the user for scientific accuracy and functionality prior to publication. Claude's role was primarily in documentation generation, code formatting, and standardization, not in generating novel scientific methodology or analysis approaches.

Tools and software used:
- FastQC, Fastp (QC assessment)
- SPAdes, QUAST (Assembly)
- BBMap, BWA, Samtools (Alignment)
- RAxML-NG, ModelTest-NG (Phylogenetics)
- 3Seq (Recombination detection)
- Matplotlib, Seaborn (Visualization)
