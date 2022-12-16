#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH -t 4:00:00
#SBATCH -A reedkm
#SBATCH -p agsmall
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIl
#SBATCH -e Combine_miRNA_Reads_%j.err
#SBATCH -o Combine_miRNA_Reads_%j.out

# Combine the trimmed miRNA reads from the multiple samples into one single
# joint file. This is to be use as input for miRDeep2, which will use the
# mapping of all reads to th ereference genome to predict novel miRNAs.

# Paths to FASTQ input and output directories
FASTQ_IN="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_rRNA_Dep_Reads"
FASTQ_OUT="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRDeep2_Run"

# We want to exclude these two samples from the analysis. This is defined as a
# regular expression because we will use 'grep' to exclude these from a 'find'
# command.
EXCLUDE="(NCT_48_38_38_R_smRNA|NCT_72_38_R_smRNA)"

# Combine the files
cat $(find "${FASTQ_IN}" -mindepth 1 -maxdepth 1 -type f -name '*_trimmed_rRNA-dep.fastq.gz' | grep -vE "${EXCLUDE}" | sort -V) > "${FASTQ_OUT}/combined_reads.fq.gz"
