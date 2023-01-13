#!/bin/bash
# Silly little script to make symlinks to BAM files. This is so that when we
# run featureCounts, the output matrix has nice column headers, rather than
# long ugly file paths.

# Directory to where the BAM files are
SRC_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM"
# Directory to where the links will be made
DEST_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_NiceNames"

mkdir -p "${DEST_DIR}"
cd "${DEST_DIR}"

# Make the links!
for bfile in $(find "${SRC_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.bam')
do
    sname=$(basename "${bfile}" | sed -e 's/_Aligned.sortedByCoord.out.bam//g')
    ln -s "${bfile}" "./${sname}"
done
