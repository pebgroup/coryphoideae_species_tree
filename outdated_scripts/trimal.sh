#!/bin/bash

## This script runs trimal on the fasta files contained within the directory from which it is executed. 

for f in *; do
	~/apps/trimAl/source/trimal -in $f -out ../6_trimal/${f/_aligned.fasta}_trimmed.fasta -gt 0.8
done

## The -gt flag specifies the minimun number of species that need to have data for a given column in order for it to be preserved. As an alternative conider using the flag -automated1 (like this: trimal -in -out -automated1).

