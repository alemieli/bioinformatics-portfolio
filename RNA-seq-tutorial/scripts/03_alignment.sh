#!/bin/bash

# ============================================================================
# RNA-seq Analysis Pipeline - Step 3: Read Alignment
# ============================================================================
# Author: Alessandro Mieli
# Purpose: Align trimmed reads to reference genome using HISAT2
# ============================================================================

# Download and extract yeast genome index (RS64-1-1)
# wget https://fungi.ensembl.org/Saccharomyces_cerevisiae/Info/Index
# tar -xzf Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
# Building the index

hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa sacCer_index

# Align paired-end reads to genome
# -x: path to genome index
# -1: forward reads (R1)
# -2: reverse reads (R2)
# -S: output SAM file

for sample in 1M_SRR93364{68..76}; do
    hisat2 -x sacCer_index \
    -1 ${sample}_1.trimmed.fastq \
    -2 ${sample}_2.trimmed.fastq \
    -S ${sample}_aligned.sam
done

# ============================================================================
# Convert SAM to BAM, sort, and index for downstream analysis
# ============================================================================

# Convert SAM to BAM (binary, compressed format)
samtools view -S -b ${sample}_aligned.sam > ${sample}_aligned.bam

# Sort BAM file by genomic coordinates
samtools sort ${sample}_aligned.bam -o ${sample}_sorted_aligned.bam

# Index sorted BAM file for fast access
samtools index ${sample}_sorted_aligned.bam

echo "Alignment complete. Output: Sample_sorted_aligned.bam"

# Optional: Remove intermediate files to save space
# rm Sample_aligned.sam Sample_aligned.bam
