#!/bin/bash
# Script to to generate a counts matrix from the BAM files and the SAF miRNA
# annotation file.

# Path to the featureCounts program
FEATURECOUNTS="/home/reedkm/shared/RIS_Projects/Software/subread-2.0.3-Linux-x86_64/bin/featureCounts"
# Path to the SAF file
SAF="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Turkey_miRNAs_NewGenome.saf"
# Path to the output file
MICRO_OUTFILE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/Counts/miRNA_Expression.bt.NewGenome.txt"
# Path to the directory of BAM Links
BAM_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_BT_NewGenome_NiceNames"

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
