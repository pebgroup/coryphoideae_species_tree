#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Secapr
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=16

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate secapr_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/01_data

# Running secapr
secapr quality_check --input . --output .

