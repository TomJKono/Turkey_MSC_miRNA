#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 2:00:00
#SBATCH -p agsmall
#SBATCH -A reedkm
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=konox006@umn.edu
#SBATCH -e Build_Kallisto_Index_%j.err
#SBATCH -o Build_Kallisto_Index_%j.out

# Build a Kallisto index of the novel and known mature miRNA sequences
# identified in the turkey miRNA sequencing experiment. The input file is a
# FASTA with combined novel and known mature miRNA sequences discovered by
# miRDeep2.

module load gcc/7.2.0

# Define paths to program and input file
KALLISTO="/home/reedkm/shared/RIS_Projects/Software/kallisto/kallisto"
MIRNA_IN="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Kallisto_Idx/Turkey_miRNAs.fa"

# Set the K-mer length
K_LEN="15"

cd $(dirname "${MIRNA_IN}")
"${KALLISTO}" index -i "Turkey_miRNAs_k15" -k "${K_LEN}" "${MIRNA_IN}"
