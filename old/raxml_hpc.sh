#!/bin/bash

## This script runs RAxML

for f in *.fasta; do (raxmlHPC -m GTRGAMMA -f a -p 12345 -x 12345 -N 100 -s $f -n ${f}_tree); done

## "-m GTRGAMMA[X]": GTR + Optimization of substitution rates + GAMMA model of rate heterogeneity (alpha parameter will be estimated). With the optional "X" appendix you can specify a ML estimate of base frequencies.  

## "-f a": rapid Bootstrap analysis and search for best-scoring ML tree in one program run

## The -p and -x options are important for reproducibility, the number does not matter but you should take note of it (see the manual)

## "-#|-N": Specify the number of alternative runs on distinct starting trees. In combination with the "-b" option [or the -x option] this will invoke a multiple boostrap analysis. 

## The -s and -w options specify input and output names

## Be careful with the -T option, which controls the number of threads to use! According to the manual the benfit of using -T will depend on the number of site patterns, so for an average gene it is not worth setting -T to more than 2 or at most 4, although this will depend on the model of evolution and if the sequences are nucleotides or amino-acids.

## The -k option specifies that bootstrapped trees should be printed with branch lengths.

## RAxML help: ~/github/standard-RAxML/raxmlHPC -h | less

## RAxML manual: evince ~/github/standard-RAxML/manual/NewManual.pdf &
