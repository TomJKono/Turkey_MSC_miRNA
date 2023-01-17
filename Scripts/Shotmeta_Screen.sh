#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=15gb
#SBATCH -t 96:00:00
#SBATCH --mail-user=konox006@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A reedkm
#SBATCH -o Shotmeta_Screen_%j.out
#SBATCH -e Shotmeta_Screen_%j.err
#SBATCH -p ag2tb

# Run Todd's shotmeta workflow for metagenomic profiling.

echo "["`date`"] Script start."

#######################################################################
# EDIT BELOW HERE
#######################################################################

# This is the folder you want to create, which will contain
# all of the results from running this pipeline. It should NOT
# already exist. If this folder already exists, this script will
# fail.
PROJ_DIR="/scratch.global/konox006/Reed_Project_016_Shotmeta"

# This is the folder that contains all the fastq file you want
# to analyze. The script will look for any files that end in
# .fastq or .fastq.gz.
FASTQ_DIR="/home/reedkm/data_release/umgc/novaseq/210825_A00223_0632_AHJHJ5DRXY/Reed_Project_016"

#######################################################################
# EDIT ABOVE HERE
#######################################################################






# Test whether the PROJ_DIR already exists
if [ -d ${PROJ_DIR} ]; then
    # The PROJ_DIR already exists
    echo "The shotmeta project directory ${PROJ_DIR} already exists. Delete it with this command 'rm -r ${PROJ_DIR}' or specify a new PROJ_DIR."
    #exit 100
else
    echo "This shotmeta project directory is: ${PROJ_DIR}."
fi





# ---------------------------------------------------------------------
# Make symlinks with correct filenames
# ---------------------------------------------------------------------

mkdir -p $PROJ_DIR/input_fastqs
cd $PROJ_DIR/input_fastqs


ln -s ${FASTQ_DIR}/*.fastq* .

# Rename any files that have spaces -- change to underscores
find -L . -type f -name "* *" | while read file; do mv "$file" "${file// /_}"; done



# ---------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------


mkdir -p $PROJ_DIR/shotmeta
cd $PROJ_DIR/shotmeta

#export MODULEPATH=/home/lmnp/knut0297/software/modulesfiles:$MODULEPATH
#module load shotmeta/1.2
# This is our local branch of the Slurm version of shotmeta. Add it to the end
# of the PATH variable so we can execute the shotmeta scripts.
export PATH=${PATH}:/home/reedkm/shared/RIS_Projects/Software/shotmeta
# We also need to load a Perl module
module load perl/5.30.1
SHOTMETA_DB=/scratch.global/konox006/shotmeta_db_20191202

# Test it
# shotmeta -t -u -d ${SHOTMETA_DB} $(find $PROJ_DIR/input_fastqs/*.fastq*)


# Run it
shotmeta -u -d ${SHOTMETA_DB} $(find $PROJ_DIR/input_fastqs/*.fastq*)



# ---------------------------------------------------------------------
# Job summary info
# ---------------------------------------------------------------------

echo "["$(date)"] Script end."

#if [ ! -z ${PBS_JOBID+x} ]; then
#    qstat -f $PBS_JOBID
#fi



