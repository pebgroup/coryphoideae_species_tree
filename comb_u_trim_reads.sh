#!/bin/bash

## This script combines subsets of unpaired reads. Then it deletes the old unpaired reads.
## WARNING! Must be run within the right directory at the right time because of the `rm` command

for file in *_1U.fastq;
do cat $file ${file/_1U/_2U} > ${file/_1U/_UN};
done

rm *U.fastq
