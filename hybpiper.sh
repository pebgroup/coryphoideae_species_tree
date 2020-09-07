#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 16

while read name
do
~/github/HybPiper/reads_first.py --cpu 16 -r ../2_trimmed/"$name"_1P.fastq ../2_trimmed/"$name"_2P.fastq --unpaired ../2_trimmed/"$name"_UN.fastq -b /home/petoe/faststorage/target_directory/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix $name --bwa
done < namelist.txt

