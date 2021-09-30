#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Secapr_data
#SBATCH --partition normal
#SBATCH --mem-per-cpu=40g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS
#SBATCH --time=4-00:00:00


#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate secapr_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/01_data

# Running secapr
secapr quality_check --input . --output .

#Moving the results to the correct folder
mv *.zip ../00_secapr/0_data
mv *.html ../00_secapr/0_data
mv *.txt ../00_secapr/0_data
mv *.pdf ../00_secapr/0_data



