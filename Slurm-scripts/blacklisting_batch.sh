#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Blacklisting
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS 
#SBATCH --time=1-00:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate amas_env

# navigating to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/05_coverage

#Creating a list of the files
ls *trimmed.fasta > filelist.txt

#Running the samples to genes program
python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/samples2genes.py > outstats.csv

#Navigating to folder for blacklisting
cd /home/owrisberg/Coryphoideae/work_flow/06_blacklisting

#Cleaning up names of the files in the folder
for f in *.FNA; do (sed -i'.old' -e $'s/-[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*_[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*//g' $f); done
rm *.old 

# Removing blacklisted Taxa from all sequence sets and tidy up file names.
#AMAS remove -x  taxa taxa taxa taxa taxa 