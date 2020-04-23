# Coryphoideae targeted capture bioinformatics work-flow
Peter Petoe (peter.petoe@bio.au.dk), 7 April 2020

## 0. Workspace
Data folder on GIS07: `data_vol/peter/github/coryphoideae`
- `1_data`
- `secapr_quality_check`
    - `raw`
    - `trimmed`
- `2_trimmed`

## 1. Before commencing analysis
SECAPR quality_check is run on the raw data in the directory `1_data`.  
In order to invoke `secapr` first activate the environment   
`conda activate secapr_env`  
Then run SECAPR from within the `raw` directory which is also where the fastqc output will be stored  
`secapr quality_check --input ../../1_data --output .`   

A list of the fastq files is created for the directory `1_data` by running the following script from within it  
`ls *.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt`  
The second command removes the last 9 characters
 
## 2. Trimming
From within the `2_trimmed` directory run trimmomatics using the following shell script  
`bash ~/github/coryphoideae/trim_loop.sh`  
All trimmed files are found in `2_trimmed` along with the stderr output in the file `stderr_trim_loop_output.txt`  

SECAPR quality_check is now run on the trimmed data in the directory `trimmed` but first, for comparability with the first SECAPR quality_check, paired and unpaired reads are combined for each sample   
`bash ~/github/coryphoideae/comb_postrim_secapr.sh`   
This script should be run from within `trimmed` and results in the creation of a subfolder `secapr_postrim` which contains the combined files.   
Now run SECAPR from within the directory `secapr_postrim`  
`secapr quality_check --input . --ouput ../../secapr_quality_check/trimmed`  
The temporary directory `secapr_postrim` can now be deleted to save storage space as the combined paired and unpaired reads are no longer need. The new fastqc output is stored within `secapr_quality_check/trimmed`  
 
## 3. Assembly (HybPiper)
### Combine unpaired reads into a single file for each sample
Run within `2_trimmed`   
`bash ~/github/coryphoideae/comb_u_trim_reads.sh`   

 
