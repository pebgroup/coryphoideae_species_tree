#!/bin/bash

#SBATCH --job-name=Trimming
#SBATCH --account Coryphoideae
#              D-HH:MM:SS
#SBATCH --time 0-05:00:00
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=24

#With these settings it takes approximately 1.2 minutes to run a single species
#OBS this was with 1 cpu and with 8 mem per cpu. 

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda trimmomatic environment 
conda activate trimmomatic_env

#going to data folder
cd /home/owrisberg/Coryphoideae/work_flow/01_data


# Creating namelist in data folder
ls *R1.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt

# Running the trimmomatic
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/trimmomatic_genomedk_sh.sh

mv *_1*.fastq ../02_trimmed
mv *_2*.fastq ../02_trimmed
