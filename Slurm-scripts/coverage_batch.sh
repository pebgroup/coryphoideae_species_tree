#!/bin/bash

#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1
#SBATCH --account Coryphoideae

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Going to folder with data
cd /home/owrisberg/Coryphoideae/test_data/3_hybpiper

# Running Wolfs Coverage tester
while read name; do /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/coverage.py $name; done < namelist.txt