#!/bin/bash

## Run Trimmomatic for a list of fastq files contained in namelist.txt 
## Removes adapters and performs quality trimming using MAXINFO with high strictness parameter values
## CONCLUSION: This script retains 25% more reads than 3trim_loop.sh despite high MAXINFO strictness and similar fastqc results.

while read name;
do java -jar /home/au297565/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33\
 ../../1_data/${name}_R1.fastq ../../1_data/${name}_R2.fastq\
 -baseout ${name}.fastq\
 ILLUMINACLIP:/home/au297565/apps/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
 LEADING:30\
 TRAILING:30\
 MAXINFO:40:0.8\
 MINLEN:36\
 2>> stderr_trim_loop_output.txt
done < ../../1_data/namelist.txt
