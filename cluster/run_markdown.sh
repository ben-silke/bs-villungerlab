#!/bin/bash
#SBATCH --output /nobackup/lab_villunger/bsilke/logs/run_markdown_%j.log # log file location (stdout), %j stands for unique job ID
#SBATCH --job-name=run_markdown_human
#SBATCH --partition=tinyq # job queue where the job is submitted to
#SBATCH --qos=tinyq # qos must match the partition
#SBATCH --nodes=1 # number of physical nodes
#SBATCH --ntasks=1 # 1 task
#SBATCH --cpus-per-task=1 # 1 task on 1 CPU
#SBATCH --time=02:00:00 # Job time is max 2 hours (it's cancelled if it's still executing after 2 hours)
# Optional parameters
#SBATCH --mem=40000 # using 40gb of memory
#SBATCH --error /nobackup/lab_villunger/bsilke/logs/run_markdown_%j.err # error log file location (stderr), %j stands for unique job ID
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
module load bio/R-bundle-Bioconductor/
module load Pandoc/3.1.5
# *** run the job ***
date
Rscript ../r/src/star_analysis/local_files/Etop_r1to6_create_stardata.R
Rscript ../r/src/star_analysis/local_files/ZM_r1to6_create_stardata.R
Rscript ../r/src/star_analysis/local_files/Noc_r1to6_create_stardata.R
Rscript ../r/src/star_analysis/local_files/Nutl_r1to6_create_stardata.R
Rscript ../r/src/star_analysis/local_files/DHCB_r1to6_create_stardata.R

# sleep 1000

# Rscript ../r/src/star_analysis/local_files/knit_all_files.R
date
seff $SLURM_JOB_ID


