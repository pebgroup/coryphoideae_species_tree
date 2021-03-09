#!/bin/bash

while read name
do
/home/owrisberg/Coryphoideae/github_code/HybPiper/reads_first.py --cpu 10 --readfiles ../2_trimmed/"$name"_1P.fastq ../2_trimmed/"$name"_2P.fastq --unpaired ../2_trimmed/"$name"_UN.fastq -b /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix $name --bwa
done < namelist.txt

