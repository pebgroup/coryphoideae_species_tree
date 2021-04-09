#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Coverage
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=16

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/3_hybpiper

# Running Wolfs Coverage tester
while read name; do python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/coverage.py $name; done < namelist.txt

