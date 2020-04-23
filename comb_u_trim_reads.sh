#!/bin/bash

## This script combines subsets of unpaired reads
## Run this script before running HybPiper

for file in *_1U.fastq;
do cat $file ${file/_1U/_2U} > ${file/_1U/_UN};
done

