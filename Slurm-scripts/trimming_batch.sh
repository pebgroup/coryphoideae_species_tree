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

#going to data folder
cd /home/owrisberg/Coryphoideae/work_flow/01_data


# Creating namelist in data folder
ls *R1.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed

# Running the trimmomatic
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/trimmomatic_genomedk_sh.sh

#going to data folder
cd /home/owrisberg/Coryphoideae/work_flow/01_data

rm namelist.txt