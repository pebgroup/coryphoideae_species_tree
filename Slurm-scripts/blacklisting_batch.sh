#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Coverage
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS 
#SBATCH --time=1-00:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

# navigating to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/05_coverage

#Creating a list of the files
ls *trimmed.fasta > filelist.txt

#Running the samples to genes program
python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/samples2genes.py > outstats.csv

#Navigating to folder for blacklisting
cd /home/owrisberg/Coryphoideae/work_flow/06_blacklisting

