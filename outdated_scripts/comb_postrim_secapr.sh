#!/bin/bash

## This script combines paired and unpaired reads for forward and reverse reads respectively for each species 
## Run this program before post-trimming secapr quality_check for comparability before and after trimming
## Run from within directory containing trimmed fastq files

mkdir secapr_postrim

for file in *_1P.fastq;
do cat $file ${file/_1P/_1U} > secapr_postrim/${file/_1P/_1PU};
done

for file in *_2P.fastq;
do cat $file ${file/_2P/_2U} > secapr_postrim/${file/_2P/_2PU};
done
