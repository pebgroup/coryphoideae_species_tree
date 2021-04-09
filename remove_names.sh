#!/bin/bash

#Script for removing inferior quality files

rm *.clean_*

for f in $(cat /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_1_data_2_rm_seqs.csv)
do
rm "$f"
done

rm *.sh
