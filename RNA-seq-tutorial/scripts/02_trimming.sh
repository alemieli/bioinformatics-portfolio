#!/bin/bash

# ============================================================================
# RNA-seq Analysis Pipeline - Step 2: Read Trimming
# ============================================================================
# Author: Alessandro Mieli
# Purpose: Remove adapters and low-quality bases using Trimmomatic
# ============================================================================

# For paired-end data
# Adjust sample names and adapter file as needed

for sample in 1M_SRR93364{68..76}; do
    trimmomatic PE -threads 4 -phred33 \
    ${sample}_1.fastq ${sample}_2.fastq. \
    ${sample}_1.trimmed.fastq ${sample}_1.unpaired.fastq \
    ${sample}_2.trimmed.fastq ${sample}_2.unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
done

# ============================================================================
# Parameters explained:
# PE: paired-end mode
# -threads 4: use 4 CPU threads for faster processing
# -phred33: quality score encoding
# ILLUMINACLIP: remove Illumina adapters (TruSeq3-PE.fa)
# LEADING:3: remove low-quality bases from start (quality < 3)
# TRAILING:3: remove low-quality bases from end (quality < 3)
# SLIDINGWINDOW:4:15: scan reads with 4-base window, cut when avg quality < 15
# MINLEN:36: discard reads shorter than 36 bases after trimming
# ============================================================================

# Run FastQC again on trimmed files to verify quality improvement
mkdir fastqc_trimmed_reports_Sample1

fastqc -o fastqc_trimmed_Sample1_reports *.trimmed.fastq

echo "Analysis complete for *.trimmed.fastq"
