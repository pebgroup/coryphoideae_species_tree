# Coryphoideae targeted capture bioinformatics work-flow
Peter Petoe (peter.petoe@bio.au.dk), 7 April 2020

## 0. Workspace
Data directory on cluster: `/data_vol/peter/coryphoideae_species_tree` 
- `fastqc`
    - `secapr_1_data`
    - `secapr_2_trimmed`
- `1_data`
- `2_trimmed`
- `3_hybpiper`
- `4_seqs`
- `5_alignments`
- `6_trimal` 
- `7_raxml`   
GitHub project repo clone on cluster: `~/github/coryphoideae_species_tree`   
GitHub HybPiper repo clone on cluster: `~/github/hybpiper`   
Directory containing target file on cluster: `/data_vol/peter/target`   
 
## 0. Downloading and renaming data 
Transfer the first set of sequence files from storage device to `1_data`. These files need renaming among other things. Within `1_data` run    
`rm *.xml *.csv Undtermined*`   
`gunzip *.gz`   
`cp ~/github/coryphoideae_species_tree/names_1_data_1_rename_seqs.csv .`   
Rename files   
`cat names_1_data_1_rename_seqs.csv | while IFS=, read orig new; do echo mv "$orig" "$new"; done`   
The above line introduces hidden characters ($'\r). To remove these run the following   
`rm names_1_data_1_rename_seqs.csv`   
`for f in *; do echo mv "$f" "$(sed 's/[^0-9A-Za-z_.]/_/g' <<< "$f")"; done`   
`rename -n 's/(.*).{1}/$1/' *`  
The last line removes an underscore, which was introduced when hidden characters were removed (for changes to take effect remove the -n flag).    

Transfer the remaining sequences to `1_data`. They will already have been renamed.   
Remove files which are of inferior quality   
`rm *.clean_*`   
`for f in $(cat ~/github/coryphoideae_species_tree/names_1_data_2_rm_seqs.csv); do rm "$f"; done`   
`rm *.sh`   
 
## 1. Before commencing analysis
SECAPR quality check is run on the raw data in `1_data`.   
In order to invoke `secapr` first activate the environment   
`conda activate secapr_env`  
Then run SECAPR from within the directory   
`secapr quality_check --input . --output .`  
Move the secapr files to `secapr_1_data`   

A list of the fastq files is created for the directory `1_data` by running the following script from within it   
`ls *R1.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt`  
The second command removes the last 9 characters. This list is needed for running Trimmomatics on all the names in the list.   

## 2. Trimming
From within the `2_trimmed` directory run trimmomatics on the sequences in `1_data`   
`bash ~/github/coryphoideae_species_tree/trimmomatic_sh`   
All trimmed files are found in `2_trimmed` along with the stderr output in the file `stderr_trim_loop_output.txt`  

SECAPR quality check is now run on the trimmed data in the directory `2_trimmed` but first, for comparability with the raw SECAPR quality check, paired and unpaired reads are combined for each sample. Within `2_trimmed` run   
`bash ~/github/coryphoideae_species_tree/comb_postrim_secapr.sh`   
This script should be run from within `2_trimmed` and results in the creation of a subdirectory `secapr_postrim` that contains the combined files.   
Now run SECAPR from within `secapr_postrim`   
`secapr quality_check --input . --output .`  
Transfer the FastQC files to `fastqc/secapr_2_trimmed` and delete `secapr_postrim` with all the combined files to save storage.   
 
## 3. Assembly (HybPiper)
### Combine unpaired reads into a single file for each sample
Run the following script within `2_trimmed`   
`bash ~/github/coryphoideae_species_tree/comb_u_trim_reads.sh`  
This merges `####_1U.fastq` and `####_2U.fastq` into `####_UN.fastq`   
 
### Generate namelist
A list of the fastq files within the directory `2_trimmed` is created in the directory `3_hybpiper` by running the following script from within `2_trimmed`   
`ls *1P.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > ../3_hybpiper/namelist.txt; rm namelist_temp.txt`
The second command removes the last 9 characters. This list is needed for running hybpiper on all the listed names.   
If all trimmed data are to go into hybpiper the namelist within `1_data` can alternatively be copied to `3_hybpiper` and reused. 

### Execute HybPiper
From within `3_hybpiper` run   
`bash ~/github/coryphoideae_species_tree/hybpiper.sh` 
 
### Get assembly stats 
From within `3_hybpiper` run   
`python ~/github/hybpiper/get_seq_lengths.py /data_vol/peter/target/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt`   
and   
`python ~/github/hybpiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt`   
The file `seq_lengths.txt` can be used to make a heatmap in R by running the script `gene_recovery_heatmap_ggplot.R`

### Paralogs
From within `3_hybpiper` run   
`bash ~/github/coryphoideae_species_tree/paralog.sh`   
The output from this program is saved in the file `paralog.txt`. This file lists putatively paralogous genes.   

## 4. Retrieve sequences (HybPiper) 
From within `3_hybpiper` run   
`python ~/github/hybpiper/retrieve_sequences.py /data_vol/peter/target/sidonie_target_file/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta . dna > stats_seq_retr.txt`  
and   
`mv *.FNA ../4_seqs`   
The recovered, unaligned multi-FASTA files for each gene can now be found in `4_seqs`. The text file `stats_seq_retr.txt` which contains the stdout for the programme `retrieve_sequences.py` can be used for the crude exclusion of genes based on number of retrieved sequences   

## 5. Alignment
From within `4_seqs` run MAFFT on all genes 
`bash ~/github/coryphoideae_species_tree/mafft.sh`   
Aligned genes are found within `5_alignments`

Visualize single gene alignments with AliView. Launch the program with command: `aliview`   
Visualize multiple gene alignments with Geneious. Launch from Nautilus.  

## 6. Gap trimming
From within `5_alignments` run 
`bash ~/github/coryphoideae_species_tree/trimal.sh`   
The trimmed alignments can now be found in `6_trimal`   

## 7. RAxML gene tree inference 
From within `6_trimal` run   
`bash ~/github/coryphoideae_species_tree/raxml_ng.sh`  
Move RAxML output to `7_raxml`   


