##!/bin/bash

#This script is meant to produce the necessary file structure for the workflow within a directory

#Building main folders
mkdir 00_secapr
mkdir 01_data
mkdir 02_trimmed
mkdir 03_hybpiper
mkdir 04_coverage
mkdir 05_blacklisting
mkdir 06_alignment
mkdir 07_optrimal
mkdir 08_cialign
mkdir 09_mapping
mkdir 10_tree_building
mkdir 11_dating_the_tree

#Building sub folders
#00
mkdir 00_secapr/0_data
mkdir 00_secapr/1_trimmed
mkdir 08_cialign/TAPER


# Creating the required environments based on the files in the environment folder
for ENV_FILE in ./environments/*
do
	name_var="${ENV_FILE#./environments/}"
	conda env create --name "${name_var%_env.txt}"  --file $ENV_FILE
done
