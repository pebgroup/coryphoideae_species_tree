#!/bin/bash

#SBATCH --job-name=Secapr_data_array
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --account=Coryphoideae
#              D-HH:MM:SS This is the max running time allowed on GDK
#SBATCH --time=0-01:00:00
#SBATCH --array 0-205%20

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda Secapr environment 
conda activate secapr_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/03_hybpiper

#uncomment to analyse all *.fq.gz files in directory
inputfiles=( $( ls /home/owrisberg/Coryphoideae/work_flow/01_data | grep *.fastq ) )

#Program to run in array:
secapr quality_check --${inputfiles[$SLURM_ARRAY_TASK_ID]} . --${inputfiles[$SLURM_ARRAY_TASK_ID]} .

#Getting the sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt

#Calculating stats based on sequence lengths
python /home/owrisberg/Coryphoideae/github_code/HybPiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt


secapr quality_check --input . --output .