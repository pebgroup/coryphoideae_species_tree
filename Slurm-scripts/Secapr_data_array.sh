#!/bin/bash

#SBATCH --job-name=Secapr_data_array
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --account=Coryphoideae
#              D-HH:MM:SS This is the max running time allowed on GDK
#SBATCH --time=0-01:00:00
#SBATCH --array 0-205%20

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda Secapr environment 
conda activate secapr_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/03_hybpiper

#uncomment to analyse all *.fq.gz files in directory
inputfiles=( $( ls /home/owrisberg/Coryphoideae/work_flow/01_data | grep 3*.fastq ) )

#Program to run in array:
secapr quality_check --input ${inputfiles[$SLURM_ARRAY_TASK_ID]} --output ${inputfiles[$SLURM_ARRAY_TASK_ID]}