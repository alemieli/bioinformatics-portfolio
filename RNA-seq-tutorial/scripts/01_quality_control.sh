#!/bin/bash

# ============================================================================
# RNA-seq Analysis Pipeline - Step 1: Quality Control
# ============================================================================
# Author: Alessandro Mieli
# Purpose: Assess quality of raw FASTQ files using FastQC
# ============================================================================

# Create output directory for QC reports
mkdir -p fastqc_reports

# Run FastQC on all FASTQ files
# -o: output directory
# *.fastq: process all fastq files in current directory

fastqc -o fastqc_reports *.fastq

# Alternative: Run FastQC on specific paired-end files
# fastqc -o fastqc_reports Sample1_R1.fastq Sample1_R2.fastq

echo "Quality control complete. Check fastqc_reports/ for HTML reports."

# ============================================================================
# Expected outputs:
# - HTML reports with quality metrics for each sample
# - Per-base sequence quality plots
# - GC content distribution
# - Adapter contamination assessment
# ============================================================================
