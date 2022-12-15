#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 12:00:00
#SBATCH -A reedkm
#SBATCH -p agsmall
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Cutadapt_%j.err
#SBATCH -o Cutadapt_%j.out

# Run cutadapt on the R1 reads from the miRNA dataset. I am not sure why R2
# were generated, but the kit only has usable information in the R1 files.

# Load the python module and source our Conda environment
module load python3/3.8.3_anaconda2020.07_mamba
source "${HOME}/.bashrc"
conda activate /home/riss/konox006/conda_envs/cutadapt_env

# Define the path to the input reads
FASTQ_IN="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Raw_Data"
# Define the path to the output
FASTQ_OUT="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_Reads"

# Define the run parameters of the cutadapt run. These are lifted directly from
# the Takara SMARTer smRNAseq data analysis document
MIN_LEN="15"
CROP="3"
ADAPTER_SEQ="AAAAAAAAAA"

mkdir -p "${FASTQ_OUT}"

# Run the cutadapt prog on all R1s
#   Note that we need the trailing slash on the 'find' directory because the
#   FASTQ_IN path is actually a symlink.
for r1 in $(find "${FASTQ_IN}/" -mindepth 1 -maxdepth 1 -name '*R1_001.fastq.gz')
do
    # Use the file name to identify the sample name
    sname=$(basename "${r1}" | sed -e 's/_R1_001\.fastq\.gz//g')
    cutadapt \
        -j "${SLURM_CPUS_PER_TASK}" \
        -m "${MIN_LEN}" \
        -u "${CROP}" \
        -a "${ADAPTER_SEQ}" \
        -o "${FASTQ_OUT}/${sname}_trimmed.fastq.gz" \
        "${r1}"
done
