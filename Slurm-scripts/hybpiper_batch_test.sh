#!/bin/bash

#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1
#SBATCH --account Coryphoideae
#              D-HH:MM:SS
#SBATCH --time 0-12:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

# Combining unpaired reads to a single file
cd /home/owrisberg/Coryphoideae/work_flow/02_trimmed
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/comb_u_trim_reads.sh

#Creatint a namelist for hybpiper based on the reads in 02_trimmed
ls *1P.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > ../3_hybpiper/namelist.txt; rm namelist_temp.txt

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/03_hybpiper

#Running the hybpiper main functionality
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/hybpiper_gdk.sh 

#Getting the sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt

#Calculating stats based on sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt


