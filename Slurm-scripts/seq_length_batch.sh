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
cd /home/owrisberg/Coryphoideae/work_flow/3_hybpiper

#Getting the sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt