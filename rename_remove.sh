#!/bin/bash

#Script for renaming species in 1_data based on names_1_data_rename_seqs.csv

#Removing unnecessary files 
rm *.xml *.csv

#Unpacking zip files
gunzip *gz

#Copying renaming list into folder
cp /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_1_data_1_rename_seqs.csv .

#Renaming the files
cat names_1_data_1_rename_seqs.csv | while IFS=, read orig new; do mv "$orig" "$new"; done

#Removing the renaming list
rm names_1_data_1_rename_seqs.csv

#Removing hidden characters which were created by the renaming script
for f in *; do echo mv "$f" "$(sed 's/[^0-9A-Za-z_.]/_/g' <<< "$f")"; done

rename -n 's/(.*).{1}/$1/' *

# Removing files with bad clean in name
rm *.clean_*

# Removing files with in names_1_data_2_rm_seqs.csv
for f in $(cat /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_1_data_2_rm_seqs.csv)
do
rm "$f"
done

#Removing leftover shell scripts within datafolder 
rm *.sh