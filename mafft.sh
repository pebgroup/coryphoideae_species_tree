#!/bin/bash

## This script runs mafft on the fasta files contained within the directory from which it is executed. Remember to adjust number of threads. 
## From the manual: "*L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information): mafft --localpair --maxiterate 1000 input [> output] linsi input [> output]". These settings are recommended by Matt Johnson in the KewHybSeqWorkshop.

for f in *
do 
linsi --thread 16 $f > ../5_alignments/${f/.FNA}_aligned.fasta
done

## Modified from Wolf's repo
# `for f in *; do (linsi --thread 16 $f > ../5_alignments/${f/.FNA}_aligned.fasta); done`

