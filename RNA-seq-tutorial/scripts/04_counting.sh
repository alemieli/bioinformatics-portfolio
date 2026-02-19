#!/bin/bash

# ============================================================================
# RNA-seq Analysis Pipeline - Step 4: Gene Quantification
# ============================================================================
# Author: Alessandro Mieli
# Purpose: Count reads mapped to genes using featureCounts
# ============================================================================

# Download genome annotation file (GTF) if not already available
# https://fungi.ensembl.org/Saccharomyces_cerevisiae/Info/Index
# gunzip Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz

# Count reads for a single sample
# -p: paired-end data
# -a: annotation file (GTF format)
# -o: output file name

for sample in 1M_SRR93364{68..76}; do
    featureCounts -p \
    -a Saccharomyces_cerevisiae.R64-1-1.115.gtf \
    -o ${sample}_counts.txt \
    ${sample}_sorted_aligned.bam
done

# ============================================================================
# For multiple samples, process all at once:
# ============================================================================

# featureCounts -p \
#   -a Saccharomyces_cerevisiae.R64-1-1.115.gtf \
#   -o all_samples_counts.txt \
#   Sample1_sorted_aligned.bam \
#   Sample2_sorted_aligned.bam \
#   Sample3_sorted_aligned.bam

echo "Gene counting complete. Output: Sample1_counts.txt"

# ============================================================================
# Post-processing:
# The output file contains gene counts along with metadata (chr, start, end, etc.)
# For DESeq2 analysis, keep only GeneID and count columns
# This can be done in Excel or with command-line tools:
# cut -f1,7 Sample1_counts.txt > Sample1_counts_clean.txt
# ============================================================================
