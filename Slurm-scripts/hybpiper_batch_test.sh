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

#Running the hybpiper main functionality
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/hybpiper_gdk.sh 

#Getting the sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt

#Calculating stats based on sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt

#Investigating potential paralog genes
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/paralog2_gdk.sh

#Sorting paralog genes
sort paralog.txt | uniq | sed 's/^.......................//' > gene_paralogs.txt

#Moving sequences to next folder 
mv *.FNA ../4_seqs

