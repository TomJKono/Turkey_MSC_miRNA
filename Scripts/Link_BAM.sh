#!/bin/bash
# Silly little script to make symlinks to BAM files. This is so that when we
# run featureCounts, the output matrix has nice column headers, rather than
# long ugly file paths.

# Directory to where the BAM files are
#   These are for STAR
# SRC_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM"
#   These are for Bowtie against the new genome
SRC_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_NewGenome"
# Directory to where the links will be made
#   For STAR
# DEST_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_NiceNames"
#   For Bowtie against new genome
DEST_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_BT_NewGenome_NiceNames"

mkdir -p "${DEST_DIR}"
cd "${DEST_DIR}"

# Make the links!
for bfile in $(find "${SRC_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*_sorted.bam')
do
    # for STAR filenames
    # sname=$(basename "${bfile}" | sed -e 's/_Aligned.sortedByCoord.out.bam//g')
    # for Bowtie filenames
    sname=$(basename "${bfile}" | sed -e 's/_sorted.bam//g')
    ln -s "${bfile}" "./${sname}"
done
