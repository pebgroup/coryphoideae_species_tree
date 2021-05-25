#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Gap_trimming
#SBATCH --partition normal
#SBATCH --mem-per-cpu=15g
#SBATCH --cpus-per-task=8
#              D-HH:MM:SS 
#SBATCH --time=0-12:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate trimal_env

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/09_optrimal

# replace n's with gaps in alignmenets - this will otherwise trip up TrimAl
for f in *.fasta; do (sed -i'.old' -e 's/n/-/g' $f); done

# change back "exo" to "exon"
for f in *.fasta; do (sed -i'.old' -e 's/exo-/exon/g' $f); done

# create summary tables for all thresholds specified
bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/pasta_taster.sh

# create summary table for the raw alignments
AMAS.py summary -f fasta -d dna -i *.fasta
mv summary.txt summary_0.txt
rm *.fasta

Rscript --vanilla /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/optrimAl.R

#Remove alignments with empty sequences
for f in *.fasta;do(python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/noempty.py $f);done
