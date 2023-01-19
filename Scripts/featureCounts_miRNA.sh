#!/bin/bash
# Script to to generate a counts matrix from the BAM files and the SAF miRNA
# annotation file.

# Path to the featureCounts program
FEATURECOUNTS="/home/reedkm/shared/RIS_Projects/Software/subread-2.0.3-Linux-x86_64/bin/featureCounts"
# Path to the SAF file
SAF="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Turkey_miRNAs.saf"
# Path to the GTF file
GTF="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Reference_Genome/anno.gtf"
# Path to the output file
MICRO_OUTFILE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Counts/miRNA_Expression.bt.txt"
MRNA_OUTFILE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Counts/miRNA_Expression.bt.mRNA_Counts.txt"
# Path to the directory of BAM Links
BAM_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_BT_NiceNames"

# Set some parameters for the counting
MIN_MAPQ="10"
STRAND="1"
THREADS="12"

# Count the microRNAs
cd "${BAM_DIR}"
"${FEATURECOUNTS}" \
    -Q "${MIN_MAPQ}" \
    -s "${STRAND}" \
    -T "${THREADS}" \
    -F SAF \
    -a "${SAF}" \
    -O \
    -o "${MICRO_OUTFILE}" \
    *

# Count the mRNAs
"${FEATURECOUNTS}" \
    -Q "${MIN_MAPQ}" \
    -s "${STRAND}" \
    -T "${THREADS}" \
    -a "${GTF}" \
    -F GTF \
    -O \
    -o "${MRNA_OUTFILE}" \
    *
