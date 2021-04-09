#!/bin/bash

## Designed to run on GenomeDK servers within an environment with trimmomatic installed
## Before running this script make sure the custom adapter (TruSeq3-PE-2.fa) is located in the folder
## Run Trimmomatic for a list of fastq files contained in namelist.txt 
## Removes adapters and performs quality trimming using MAXINFO with high strictness parameter values
## CONCLUSION: This script retains 25% more reads than 3trim_loop.sh despite high MAXINFO strictness and similar fastqc results.

while read name;
do trimmomatic PE -threads 16 -phred33\
 /home/owrisberg/Coryphoideae/work_flow/01_data/${name}_R1.fastq /home/owrisberg/Coryphoideae/work_flow/1_data/${name}_R2.fastq\
 -baseout ${name}.fastq\
 ILLUMINACLIP:/home/owrisberg/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
 LEADING:3\
 TRAILING:3\
 MAXINFO:40:0.8\
 MINLEN:36\
 2>> stderr_trim_loop_output.txt
done < ../01_data/namelist.txt
