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
#SBATCH -e Bowtie_Mapping_%j.err
#SBATCH -o Bowtie_Mapping_%j.out

# Map the trimmed and rRNA-depleted FASTQ files to the turkey genome with
# bowtie1.

module load samtools/1.14
module load parallel/20210822

# Path to bowtie program
BOWTIE="/home/reedkm/shared/RIS_Projects/Software/bowtie-1.3.1-linux-x86_64/bowtie"
# Path to reference index
REF_IDX="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Reference_Genome/Bowtie_Idx/bt_idx"

# Path to alignment directory
OUT_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/miRNA_Expression/BAM_BT"
# Path to reads directory
READS_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_rRNA_Dep_Reads"

# Define parameters for the mapping
SEED_LEN="15"
SEED_MISMATCH="0"
MAX_MISMATCH_QUAL_SUM="80"
MAX_NUM_ALN="5"

for r1 in $(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*_trimmed_rRNA-dep.fastq.gz')
do
    sname=$(basename "${r1}" | sed -e 's/_trimmed_rRNA-dep.fastq.gz//g')
    "${BOWTIE}" \
        -p "${SLURM_CPUS_PER_TASK}" \
        -q \
        -n "${SEED_MISMATCH}" \
        -e "${MAX_MISMATCH_QUAL_SUM}" \
        -l "${SEED_LEN}" \
        -m "${MAX_NUM_ALN}" \
        --best \
        --strata \
        --un "${OUT_DIR}/${sname}_un.fastq" \
        -S \
        -x "${REF_IDX}" \
        "${r1}" \
        "${OUT_DIR}/${sname}.sam"
done

# Then, convert all the same files to bam files
cd "${OUT_DIR}"
for s in $(find . -mindepth 1 -maxdepth 1 -type f -name '*.sam')
do
    echo "samtools view -hb ${s} > ${s/.sam/.bam}"
done | parallel
rm *.sam

# And compress all of the unmapped reads
parallel gzip ::: $(find . -mindepth 1 -maxdepth 1 -type f -name '*_un.fastq')
