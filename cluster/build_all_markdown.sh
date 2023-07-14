#!/bin/bash
#SBATCH --output /nobackup/lab_villunger/bsilke/logs/markdown_creation_%j.log # log file location (stdout), %j stands for unique job ID
#SBATCH --job-name=markdown_creation
#SBATCH --partition=tinyq # job queue where the job is submitted to
#SBATCH --qos=tinyq # qos must match the partition
#SBATCH --nodes=1 # number of physical nodes
#SBATCH --ntasks=1 # 1 task
#SBATCH --cpus-per-task=1 # 1 task on 1 CPU
#SBATCH --time=12:00:00 #
# Optional parameters
#SBATCH --mem=10000 # using 10gb of memory
#SBATCH --error /nobackup/lab_villunger/bsilke/logs/markdown_creation_%j.err # error log file location (stderr), %j stands for unique job ID
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
module load R
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
# *** run the job ***
date

Rscript r/src/salmon_DESeq2/markdown_files/DHCB_r1to6_create_data.R
Rscript r/src/salmon_DESeq2/markdown_files/ZM_r1to6_create_data.R
Rscript r/src/salmon_DESeq2/markdown_files/Etop_r1to6_create_data.R
Rscript r/src/salmon_DESeq2/markdown_files/Nutl_r1to6_create_data.R
Rscript r/src/salmon_DESeq2/markdown_files/Noc_r1to6_create_data.R


sleep 100

Rscript r/src/salmon_DESeq2/markdown_files/knit_all_files.R

date

seff $SLURM_JOB_ID