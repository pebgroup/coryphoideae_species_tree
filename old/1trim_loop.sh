#!/bin/bash

## Run Trimmomatic for a list of fastq files contained in namelist.txt
## Adapter trimming
## No quality trimming
## Length trimming
## TEST OUTCOME: This script adequately addresses adapter problem

while read name;
do java -jar /home/au297565/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33\
 ../../1_data/${name}_R1.fastq ../../1_data/${name}_R2.fastq\
 -baseout ${name}.fastq\
 ILLUMINACLIP:/home/au297565/apps/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
 LEADING:20\
 TRAILING:20\
 MINLEN:36\
 2>> shell_output.txt
done < ../../1_data/namelist.txt
