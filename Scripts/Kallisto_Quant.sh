#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 12:00:00
#SBATCH -p agsmall
#SBATCH -A reedkm
#SBATCH --array=0-21
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Kallisto_Quant_%j.err
#SBATCH -o Kallisto_Quant_%j.out

# Run the quantification algorithm from Kallisto on the turkey miRNA samples.
# Be sure to exclude the two samples that are not for miRNA analysis! Note that
# kallisto can only process one sample at a time, so we are running this with
# an array.

module load gcc/7.2.0

# Define paths to program, miRNA reference index, and output
KALLISTO="/home/reedkm/shared/RIS_Projects/Software/kallisto/kallisto"
IDX="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Kallisto_Idx/Turkey_miRNAs_k15"
OUT_BASE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Kallisto_Quant"

# Define path to input reads and regex for exclusion
READS_IN="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_rRNA_Dep_Reads"
EXCLUDE="(NCT_48_38_38_R_smRNA|NCT_72_38_R_smRNA)"

# Define parameters for the quantification
#   NOTE, these would be better estimated from the bioanalyzer!
FRAG_LEN="50"
FRAG_LEN_SD="10"

# Build an array of input filenames
FASTQ_ARR=($(find "${READS_IN}" -mindepth 1 -maxdepth 1 -type f -name '*_trimmed_rRNA-dep.fastq.gz' | grep -vE "${EXCLUDE}" | sort -V))
CURR_SAMP="${FASTQ_ARR[${SLURM_ARRAY_TASK_ID}]}"
SNAME=$(basename "${CURR_SAMP}" | sed -e 's/_trimmed_rRNA-dep.fastq.gz//g')

# Build a full output directory from the sample name and output base
OUT_DIR="${OUT_BASE}/${SNAME}"

# Run the quant!
mkdir -p "${OUT_DIR}"
"${KALLISTO}" quant \
    -i "${IDX}" \
    -o "${OUT_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    --single \
    -l "${FRAG_LEN}" \
    -s "${FRAG_LEN_SD}" \
    "${CURR_SAMP}"
