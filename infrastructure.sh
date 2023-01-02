##!/bin/bash

#This script is meant to produce the necessary file structure for the workflow within a directory

#Building main folders
mkdir 00_secapr_TEST
mkdir 01_data_TEST
mkdir 02_trimmed_TEST
mkdir 03_hybpiper_TEST
mkdir 04_coverage_TEST
mkdir 05_blacklisting_TEST
mkdir 06_alignment_TEST
mkdir 07_mapping_TEST
mkdir 08_optrimal_TEST
mkdir 09_manual-edit_TEST
mkdir 10_tree_building_TEST

#Building sub folders
#00
mkdir 00_secapr/0_data_TEST
mkdir 00_secapr/1_trimmed_TEST

mkdir10_manual_edit/01_alignments_for_editing_TEST
mkdir10_manual_edit/02_edited_alignments_TEST
mkdir10_manual_edit/03_bad_alignments_TEST
mkdir10_manual_edit/04_alignments_for_trees_TEST


# Creating the required environments based on the files in the environment folder
for ENV_FILE in ./environments/*
do
	conda env create --name "TEST" + ${ENV_FILE%_env.txt}  --file $ENV_FILE
done
