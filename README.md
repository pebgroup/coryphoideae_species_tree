# Coryphoideae targeted capture bioinformatics work-flow
Peter Petoe (peter.petoe@bio.au.dk), 7 April 2020

## 0. Workspace
Data folder on GIS07: `data_vol/peter/coryphoideae`
- `1_data`
- `secapr_quality_check`
    - `raw`
    - `trimmed`
- `2_trimmed`

## 1. Preparing data for analysis
SECAPR quality_check is run on the raw data in the directory `1_data`. 
In order to invoke `secapr` first activate the environment 
`conda activate secapr_env`
Then run SECAPR from within the `raw` directory which is also where the fastqc output will be stored
`secapr quality_check --input ../../1_data --output .` 

A list of the fastq files is created for the directory `1_data` by running the following script from within it 
`ls *.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt`
The second command removes the last 9 characters
 
## 2. Trimming
From within the `2_trimmed` directory run trimmomatics using the following script
`bash ~/github/coryphoideae/trim_loop.sh`
All trimmed files are found in `2_trimmed` along with the stderr output in the file `stderr_trim_loop_output.txt`

