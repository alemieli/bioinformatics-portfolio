# RNA-seq Data Analysis Pipeline

A complete bioinformatics workflow for analyzing RNA sequencing data, from raw reads to differential gene expression analysis.

## Overview

This project demonstrates a standard RNA-seq analysis pipeline following best practices for transcriptomic data processing. The workflow processes paired-end Illumina sequencing data through quality control, read alignment, gene quantification, and statistical analysis to identify differentially expressed genes.

## Dataset

**Source:** Public RNA-seq dataset from [Bioinfogp CNB-CSIC](https://bioinfogp.cnb.csic.es/files/samples/rnaseq/)

**Organism:** *Saccharomyces cerevisiae S288c*

**Reference** Hakkaart X, Liu Y, Hulst M, El Masoudi A, Peuscher E, Pronk J, van Gulik W, Daran-Lapujade P. Physiological responses of Saccharomyces cerevisiae to industrially relevant conditions: Slow growth, low pH, and high CO2 levels. Biotechnol Bioeng. 2020 Mar;117(3):721-735. doi: 10.1002/bit.27210. Epub 2020 Jan 22.

**Design:** Paired-end sequencing, multiple treatment groups with biological replicates

## Workflow

### 1. Quality Control (FastQC)
- Assessment of raw read quality
- Identification of adapter contamination
- Evaluation of base quality scores and GC content

### 2. Read Trimming (Trimmomatic)
- Removal of adapter sequences
- Trimming of low-quality bases
- Filtering reads below minimum length threshold

### 3. Read Alignment (HISAT2)
- Mapping trimmed reads to the yeast reference genome (R64-1-1.115)
- Generation of alignment files in SAM/BAM format

### 4. Read Quantification (featureCounts)
- Counting reads mapped to genomic features (genes)
- Generation of gene count matrix for downstream analysis

### 5. Differential Expression Analysis (DESeq2 in R)
- Normalization of count data
- Statistical testing for differential expression
- Generation of visualizations (PCA plots, heatmaps, volcano plots)

## Key Results

- Identified differentially expressed genes (DEGs) between experimental conditions
- Generated quality control reports confirming high data quality
- Performed pairwise comparisons with statistical significance testing (FDR < 0.05)
- Created visual representations of gene expression patterns

## Technologies Used

**Command Line Tools:**
- FastQC
- Trimmomatic
- HISAT2
- Samtools
- featureCounts (Subread package)

**R/Bioconductor:**
- DESeq2
- pheatmap
- ggplot2

## Repository Structure

```
RNA-seq-tutorial/
├── README.md
└── scripts/
    ├── 01_quality_control.sh
    ├── 02_trimming.sh
    ├── 03_alignment.sh
    ├── 04_counting.sh
    └── 05_DESeq2_analysis.R
```

## References

This pipeline is based on standard RNA-seq analysis protocols:
- Shouib, R. et al. (2025). A Guide to Basic RNA Sequencing Data Processing and Transcriptomic Analysis. Bio-protocol 15(9): e5295.
- Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.

## Notes

This project was completed as part of my bioinformatics training to develop skills in transcriptomic data analysis and prepare for work in computational biology and genomics research.
