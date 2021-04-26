#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Intronerate
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS
#SBATCH --time=1-00:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/03_hybpiper

# Running Intronerate
while read name; do python3 /home/owrisberg/Coryphoideae/github_code/HybPiper/intronerate.py --prefix $name &>> intronerate_out.txt; done < namelist.txt
