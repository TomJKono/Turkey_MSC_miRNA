#!/bin/bash
# Script to to generate a counts matrix from the BAM files and the SAF miRNA
# annotation file.

# Path to the featureCounts program
FEATURECOUNTS="/home/reedkm/shared/RIS_Projects/Software/subread-2.0.3-Linux-x86_64/bin/featureCounts"
# Path to the SAF file
SAF="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Turkey_miRNAs.saf"
# Path to the output file
OUTFILE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Counts/miRNA_Expression.txt"
# Path to the directory of BAM Links
BAM_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_NiceNames"

# Set some parameters for the counting
MIN_MAPQ="10"
STRAND="1"
THREADS="12"

# Count it up!
cd "${BAM_DIR}"
"${FEATURECOUNTS}" \
    -Q "${MIN_MAPQ}" \
    -s "${STRAND}" \
    -t "${THREADS}" \
    -F SAF \
    -a "${SAF}" \
    -o "${OUTFILE}" \
    *
