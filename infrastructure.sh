##!/bin/bash

#This script is meant to produce the necessary file structure for the workflow within a directory

#Building main folders
echo mkdir 00_secapr_TEST
echo mkdir 01_data_TEST
echo mkdir 02_trimmed_TEST
echo mkdir 03_hybpiper_TEST
echo mkdir 04_coverage_TEST
echo mkdir 05_blacklisting_TEST
echo mkdir 06_alignment_TEST
echo mkdir 07_mapping_TEST
echo mkdir 08_optrimal_TEST
echo mkdir 09_manual-edit_TEST
echo mkdir 10_tree_building_TEST

#Building sub folders
#00
echo mkdir 00_secapr/0_data_TEST
echo mkdir 00_secapr/1_trimmed_TEST

echo mkdir10_manual_edit/01_alignments_for_editing_TEST
echo mkdir10_manual_edit/02_edited_alignments_TEST
echo mkdir10_manual_edit/03_bad_alignments_TEST
echo mkdir10_manual_edit/04_alignments_for_trees_TEST


# Creating the required environments based on the files in the environment folder
for ENV_FILE in ./environments/*
do
	echo conda env create --name "TEST"+${ENV_FILE%_env.txt}  --file $ENV_FILE
done
