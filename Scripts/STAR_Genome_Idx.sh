#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 12:00:00
#SBATCH -A reedkm
#SBATCH -p agsmall
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e STAR_Genome_Idx_%j.err
#SBATCH -o STAR_Genome_Idx_%j.out

# Index the genome for mapping with STAR. We will use some special options to
# make sure we can get the index is suitable for our small RNA reads.
module load star/2.7.1a

# Define paths to the reference genome index and annotation file to use
#   This is "Turkey_5.1" from Ensembl release 108.
REF_BASE="/home/reedkm/shared/RIS_Projects/Reed_Project_016/Reference_Genome"
REF_FA="${REF_BASE}/Meleagris_gallopavo.Turkey_5.1.dna_rm.toplevel.fa.gz"
REF_GTF="${REF_BASE}/Meleagris_gallopavo.Turkey_5.1.108.chr.gtf.gz"
# Set the output directory for the STAR index files
STAR_DIR="${REF_BASE}/STAR_Idx"

# Set the parameters of the STAR indexing
SPLICE_OVERHANG="1"

# Unzip the genome because STAR can't read zipped references
gzip -cd "${REF_FA}" > "${REF_BASE}/genome.fa"

# Make the output directory and run the command
mkdir -p "${STAR_DIR}"
STAR \
    --runMode genomeGenerate \
    --genomeDir "${STAR_DIR}" \
    --genomeFastaFiles "${REF_BASE}/genome.fa" \
    --sjdbGTFfile "${REF_GTF}" \
    --sjdbOverhang "${SPLICE_OVERHANG}" \
    --runThreadN "${SLURM_CPUS_PER_TASK}"
