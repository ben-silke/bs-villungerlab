#!/bin/bash
#SBATCH --output /nobackup/lab_villunger/bsilke/logs/star_index_gencode_%j.log # log file location (stdout), %j stands for unique job ID
#SBATCH --job-name=star_index_gencode_human
#SBATCH --partition=tinyq # job queue where the job is submitted to
#SBATCH --qos=tinyq # qos must match the partition
#SBATCH --nodes=1 # number of physical nodes
#SBATCH --ntasks=1 # 1 task
#SBATCH --cpus-per-task=8 # 1 task on 1 CPU
#SBATCH --time=02:00:00 # Job time is max 2 hours (it's cancelled if it's still executing after 2 hours)
# Optional parameters
#SBATCH --mem=400000 # using 400gb of memory
#SBATCH --error /nobackup/lab_villunger/bsilke/logs/star_index_gencode_%j.err # error log file location (stderr), %j stands for unique job ID
#SBATCH --mail-type=end # send an email when this job ends
#SBATCH --mail-user=bsilke@cemm.at # email your CeMM account

# (Optional) Examples of the Slurm environmental variables. Provide additional information
echo "======================"
echo $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NAME
echo $SLURM_JOB_PARTITION
echo $SLURM_NTASKS
echo $SLURM_NPROCS
echo $SLURM_JOB_ID
echo $SLURM_JOB_NUM_NODES
echo $SLURM_NODELIST
echo $SLURM_CPUS_ON_NODE
echo "======================"

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load STAR

# *** run the job ***
date

STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir "/nobackup/lab_villunger/bsilke/gtf_star_index_gencode_tscrpt" \dd
# this doesnt work u need to cat "gencode.v44.transcripts.fa"
--genomeFastaFiles "/nobackup/lab_villunger/bsilke/references/star/GRCh38.primary_assembly.genome.fa"  \
--sjdbGTFfile "/nobackup/lab_villunger/bsilke/references/star/gencode.v44.annotation.gtf" \
--limitGenomeGenerateRAM 350000000000
date

seff $SLURM_JOB_ID

# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/ .
# https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/