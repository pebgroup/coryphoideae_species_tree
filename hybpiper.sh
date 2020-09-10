#!/bin/bash

while read name
do
~/github/hybpiper/reads_first.py --cpu 7 -r ../2_trimmed/"$name"_1P.fastq ../2_trimmed/"$name"_2P.fastq --unpaired ../2_trimmed/"$name"_UN.fastq -b /home/au297565/Documents/temp/sidonie_target_file/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix $name --bwa
done < namelist.txt

