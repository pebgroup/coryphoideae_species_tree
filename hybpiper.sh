#!/bin/bash

target_file=PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta
cpu=16

while read name
do
~/github/hybpiper/reads_first.py --cpu "$cpu" -r ../2_trimmed/"$name"_1P.fastq ../2_trimmed/"$name"_2P.fastq --unpaired ../2_trimmed/"$name"_UN.fastq -b data_vol/peter/target/"$target_file" --prefix $name --bwa
done < namelist.txt

