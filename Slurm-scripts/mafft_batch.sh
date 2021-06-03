#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Alignment_mafft
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS 
#SBATCH --time=2-00:00:00

#This line enables the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate mafft_env

# navigating to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/05_blacklisting


## This script runs mafft on the fasta files contained within the directory from which it is executed. Remember to adjust number of threads. 
## From the manual: "*L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information): mafft --localpair --maxiterate 1000 input [> output] linsi input [> output]". These settings are recommended by Matt Johnson in the KewHybSeqWorkshop.

# Pirada's method of doing mafft
for f in *.FNA; do (echo ${f/.FNA} >> ./genenames.txt); done # making a list of the genes

parallel --eta "mafft --localpair --adjustdirectionaccurately --maxiterate 1000 {}.FNA > /home/owrisberg/Coryphoideae/work_flow/07_alignment/${f}_aligned.fasta" :::: genenames.txt