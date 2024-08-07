#!/bin/bash

#SBATCH --job-name=Seq_lengths_and_stats
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4
#SBATCH --account=Coryphoideae
#SBATCH --time=0-23:59:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Generating a new namelist
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed
ls *1P.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > ../03_hybpiper/namelist.txt; rm namelist_temp.txt

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/03_hybpiper

#Getting the sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt

#Calculating stats based on sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt
