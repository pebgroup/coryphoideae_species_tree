#!/bin/bash

#SBATCH --job-name=Trimming
#SBATCH --account Coryphoideae
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1


#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda trimmomatic environment 
conda activate trimmomatic_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed

# Running the trimmomatic
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/trimmomatic_genomedk_sh.sh

