#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Secapr_trimmed
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=16


#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate secapr_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed

#Combining the reads for secapr evalutation
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/comb_postrim_secapr.sh

#
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed/secapr_postrim

# Running secapr
secapr quality_check --input . --output .

#Moving the results to the correct folder
mv *.zip ../../00_secapr/1_trimmed
mv *.html ../../00_secapr/1_trimmed
mv *.txt ../../00_secapr/1_trimmed
mv *.pdf ../../00_secapr/1_trimmed


