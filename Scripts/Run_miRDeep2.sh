#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 48:00:00
#SBATCH -A reedkm
#SBATCH -p agsmall
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Run_miRDeep2_%j.err
#SBATCH -o Run_miRDeep2_%j.out

# Run the miRDeep2 pipeline to discover novel miRNAs in this dataset

# Load up the Conda environment for miRDeep2
module load python3/3.8.3_anaconda2020.07_mamba
source "${HOME}/.bashrc"
conda activate /home/riss/konox006/conda_envs/miRDeep2_env

# Set paths to resources needed for analysis
RUN_BASE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRDeep2_Run"
REF_IDX="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Reference_Genome/Bowtie_Idx/bt_idx"
REF_FA="${RUN_BASE}/Turkey5.1_miRDeep2_Names.fa"
CHICKEN_MATURE="${RUN_BASE}/mature_gg_miRDeep2_Names.fa"
CHICKEN_HAIRPIN="${RUN_BASE}/hairpin_gg_miRDeep2_Names.fa"
INPUT_READS="${RUN_BASE}/combined_reads.fq"
COLLAPSED_READS="${RUN_BASE}/combined_reads_Collapsed.fa"
ALN_FILE="${RUN_BASE}/reads_vs_ref.arf"

# Set the minimum length to process
MIN_LEN="18"

# Run the analysis! Start with mapping miR reads to the genome and producing
# the proprietary ARF file
cd "${RUN_BASE}"
mapper.pl \
    "${INPUT_READS}" \
    -e -h -j \
    -l "${MIN_LEN}" \
    -m \
    -p "${REF_IDX}" \
    -s "${COLLAPSED_READS}" \
    -t "${ALN_FILE}" \
    -v \
    -o "${SLURM_CPUS_PER_TASK}"

# Then run the miRDeep2 algorithm on the mapped reads file
miRDeep2.pl \
    "${COLLAPSED_READS}" \
    "${REF_FA}" \
    "${ALN_FILE}" \
    "${CHICKEN_MATURE}" \
    none \
    "${CHICKEN_HAIRPIN}" \
    -P \
    2> "${RUN_BASE}/miRDeep2_$(date +'%F_%T')_report.log"
