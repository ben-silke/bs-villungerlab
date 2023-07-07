#!/bin/bash
#SBATCH --output /nobackup/villunger/drizzotto/star_index_%j.log # log file location (stdout), %j stands for unique job ID
#SBATCH --job-name=STAR_index_human
#SBATCH --partition=tinyq # job queue where the job is submitted to
#SBATCH --qos=longq # qos must match the partition
#SBATCH --nodes=1 # number of physical nodes
#SBATCH --ntasks=1 # 1 task
#SBATCH --cpus-per-task=1 # 1 task on 1 CPU
#SBATCH --time=02:00:00 # Job time is max 2 hours (it's cancelled if it's still executing after 2 hours)
# Optional parameters
#SBATCH --mem=32000 # using 32Gb of memory
#SBATCH --error /nobackup/villunger/drizzotto/star_index_%j.err # error log file location (stderr), %j stands for unique job ID
#SBATCH --mail-type=end # send an email when this job ends
#SBATCH --mail-user=drizzotto@cemm.at # email your CeMM account

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
# May need to load star
module load STAR

# set the temporary folder to your scratch location (avoid using node local /tmp)
export TMPDIR=/nobackup/villunger/drizzotto/Rtmp

# *** run the job ***
date

# Download the ensembl genome
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/homo_sapiens "$external_drive/STAR/ensembl_genome" 

# Also need ncRNA  -- https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/ncrna/
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/homo_sapiens/ncrna "$external_drive/STAR/ensembl_genome/nc" 

# Download the annotation
rsync -av rsync://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens /nobackup/villunger/drizzotto/STAR/references/annotation



STAR \
--runThreadN 2 \
--runMode genomeGenerate \
--genomeDir /nobackup/villunger/drizzotto/STAR/references \
--genomeFastaFiles /Volumes/bs_external/villunger/references/STAR/ensembl_genome/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz /Volumes/bs_external/villunger/references/STAR/ensembl_genome/nc/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz\
# /Volumes/bs_external/villunger/references/STAR/ensemble_annotation/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gff3.gz
--sjdbGTFfile /Volumes/bs_external/villunger/references/ensemble_annotation/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gff3.gz \ # path to file which contains reference
# --outFileNamePrefix /nobackup/villunger/drizzotto/STAR/output \#
date
seff $SLURM_JOB_ID