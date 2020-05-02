#!/bin/bash

while read name
do
~/github/hybpiper/reads_first.py -r ../2_trim/4trim_loop/"$name"_1P.fastq ../2_trim/4trim_loop/"$name"_2P.fastq --unpaired ../2_trim/4trim_loop/"$name"_UN.fastq -b ~/Documents/temp/sidonie_target_file/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix "$name" --bwa
done < namelist.txt

