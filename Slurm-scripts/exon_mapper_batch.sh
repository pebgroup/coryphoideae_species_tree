#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Exon_mapper
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=12
#              D-HH:MM:SS This is the max running time allowed on GDK
#SBATCH --time=0-06:00:00


#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/07_alignment

# Running Wolfs Coverage tester
python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/exon_mapper.py
