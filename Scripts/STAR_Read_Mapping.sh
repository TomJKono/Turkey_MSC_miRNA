#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=4gb
#SBATCH --tmp=48gb
#SBATCH -t 24:00:00
#SBATCH -p agsmall
#SBATCH -A reedkm
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e STAR_Read_Mapping_%j.err
#SBATCH -o STAR_Read_Mapping_%j.out

# Map the trimmed and rRNA-depleted reads with STAR.
module load star/2.7.1a
module load samtools/1.14

# Define the paths for the STAR genome index and output locations
REF_BASE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Reference_Genome"
STAR_DIR="${REF_BASE}/STAR_Idx"
STAR_TEMP="/scratch.local/STAR_tmp"
OUT_BASE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/BAM"
READS_DIR="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Trimmed_rRNA_Dep_Reads"

# Define the parameters for STAR mapping here:
#   Some of these were borrowed from the ENCODE small RNA mapping pipeline
SEED_LMAX="10"
MULTIMAP_MAX="20"
INTRON_MAX="1"
MAX_MISMATCH_PROP="0.03"
MIN_NORM_SCORE="0"
MIN_MATCH_PROP="0"
MIN_TOTAL_MATCH="15"
OUT_SAM_UNMAPPED="Within"
MIN_SPLICE_OVERHANG="1000"
BAM_SORT_RAM="90000000000"

# Make the output directory and map the reads
mkdir -p "${OUT_BASE}"
for r1 in $(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*_trimmed_rRNA-dep.fastq.gz')
do
    sname=$(basename "${r1}" | sed -e 's/_trimmed_rRNA-dep\.fastq\.gz//g')
    STAR \
        --outTmpDir "${STAR_TEMP}" \
        --outFileNamePrefix "${OUT_BASE}/${sname}_" \
        --runMode alignReads \
        --genomeDir "${STAR_DIR}" \
        --readFilesIn "${r1}" \
        --readFilesCommand zcat \
        --runThreadN "${SLURM_CPUS_PER_TASK}" \
        --seedSearchStartLmax "${SEED_LMAX}" \
        --outFilterMultimapNmax "${MULTIMAP_MAX}" \
        --alignIntronMax "${INTRON_MAX}" \
        --outFilterMismatchNoverLmax "${MAX_MISMATCH_PROP}" \
        --outFilterScoreMinOverLread "${MIN_NORM_SCORE}" \
        --outFilterMatchNminOverLread "${MIN_MATCH_PROP}" \
        --outFilterMatchNmin "${MIN_TOTAL_MATCH}" \
        --genomeLoad NoSharedMemory \
        --outSAMunmapped "${OUT_SAM_UNMAPPED}" \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM "${BAM_SORT_RAM}" \
        --alignSJDBoverhangMin "${MIN_SPLICE_OVERHANG}"
    # And, annoyingly, we have to remove the temp dir because STAR dies if
    # the temp dir already exists.
    rm -rf "${STAR_TEMP}"
done
