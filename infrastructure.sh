##!/bin/bash

#This script is meant to produce the necessary file structure for the workflow within a directory

#Building main folders
mkdir 00_secapr
mkdir 01_data
mkdir 02_trimmed
mkdir 03_hybpiper
mkdir 04_seqs # This folder is empty?
mkdir 05_coverage
mkdir 06_blacklisting
mkdir 07_alignment
mkdir 08_mapping
mkdir 09_optrimal
mkdir 10_manual-edit
mkdir 11_tree_building

#Building sub folders
#00
mkdir 00_secapr/0_data
mkdir 00_secapr/1_trimmed

mkdir10_manual_edit/01_alignments_for_editing
mkdir10_manual_edit/02_edited_alignments
mkdir10_manual_edit/03_bad_alignments
mkdir10_manual_edit/04_alignments_for_trees