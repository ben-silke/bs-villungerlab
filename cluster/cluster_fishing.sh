#!/bin/bash
#SBATCH --output /nobackup/lab_villunger/bsilke/logs/salmon_quant_%j.log # log file location (stdout), %j stands for unique job ID
#SBATCH --job-name=salmon_quant
#SBATCH --partition=longq # job queue where the job is submitted to
#SBATCH --qos=longq # qos must match the partition
#SBATCH --nodes=1 # number of physical nodes
#SBATCH --ntasks=1 # 1 task
#SBATCH --cpus-per-task=8 # 1 task on 1 CPU
#SBATCH --time=12:00:00 #
# Optional parameters
#SBATCH --mem=40000 # using 40gb of memory
#SBATCH --error /nobackup/lab_villunger/bsilke/logs/salmon_quant_%j.err # error log file location (stderr), %j stands for unique job ID
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
module load Salmon

# *** run the job ***
date
./fishing.sh /nobackup/lab_villunger/bsilke
date

seff $SLURM_JOB_ID