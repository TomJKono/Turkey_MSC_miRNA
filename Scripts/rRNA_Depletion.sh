#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem-per-cpu=4gb
#SBATCH -p agsmall
#SBATCH -t 12:00:00
#SBATCH -A reedkm
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e rRNA_Depletion_%j.err
#SBATCH -o rRNA_Depletion_%j.out

# Use BBDuk to remove reads that have a high k-mer match to a known rRNA in
# our SILVA rRNA reference.
module load java/openjdk-17.0.2
BBDUK="/home/reedkm/shared/RIS_Projects/Software/bbmap/bbduk.sh"

# Define path to input FASTQ
FASTQ_IN="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_Reads"
# Define path to output FASTQ
FASTQ_OUT="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_rRNA_Dep_Reads"
# Define path to output reports
REP_OUT="/home/reedkm/shared/RIS_Projects/Reed_Project_016/rRNA_Dep_Reports"

# Define the parameters for rRNA trimming
SILVA_REF="/home/msistaff/public/CHURP_Deps/v0/db/SILVA_132_LSU_SSU_Ref_Dedup_Kmers_min100.fasta.gz"
MIN_LENGTH="15"
K_LENGTH="15"
EDIT_DIST="0"
THREADS="12"

mkdir -p "${FASTQ_OUT}"
for r1 in $(find "${FASTQ_IN}" -mindepth 1 -maxdepth 1 -type f -name '*.fastq.gz')
do
    sname=$(basename "${r1}" | sed -e 's/\.fastq\.gz//g')
    _JAVA_OPTIONS="-Xmx45g -Xms24g" "${BBDUK}" \
        in="${r1}" \
        out="${FASTQ_OUT}/${sname}_rRNA-dep.fastq.gz" \
        outm="${FASTQ_OUT}/${sname}_rRNA-match.fastq.gz" \
        k="${K_LENGTH}" \
        minlength="${MIN_LENGTH}" \
        editdistance="${EDIT_DIST}" \
        statscolumns="5" \
        stats="${REP_OUT}/${sname}_rRNA-stats.txt" \
        ref="${SILVA_REF}" \
        threads="${THREADS}" \
        prealloc="t"
done
