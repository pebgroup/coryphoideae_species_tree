# Coryphoideae targeted capture bioinformatics work-flow

Inherited by Oscar Wrisberg ([Oscar.wrisberg@bio.au.dk](mailto:Oscar.wrisberg@bio.au.dk)), 11th January 2021

* * *

## 0\. Workspace
**Required file structure and environments**  
In order to run this pipeline you need a directory with the following folders along with a set of conda environments.
These can be quickly produced by running `infrastructure.sh` within the desired directory.
**OBS** one of the environments is called **base** and running the infrastructure script might **erase** your own base environment, proceed at your own risk. 

- `00_secapr`
  - `0_data`
  - `1_trimmed`
- `01_data`
- `02_trimmed`
- `03_hybpiper`
- `04_coverage`
- `05_blacklisting`
- `06_alignment`
- `07_mapping`
- `08_optrimal`
- `09_manual_edit`
- `10_tree_building`

## 00\. Downloading and renaming data

Transfer all sequence files from storage device to `01_data`.  
These files need renaming and unpacking among other things. Within `01_data` run  
`rm *.xml *.csv Undtermined*`  
`gunzip *.gz`  
`cp /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_1_data_1_rename_seqs.csv .`

**Rename files**  
You can test if this script is printing the correct things before running it by including an echo between do mv, like so `do echo mv`

`cat names_1_data_1_rename_seqs.csv | while IFS=, read orig new; do mv "$orig" "$new"; done`  
The above line introduces hidden characters ($'  
). To remove these run the following.

`rm names_1_data_1_rename_seqs.csv`  
`for f in *; do echo mv "$f" "$(sed 's/[^0-9A-Za-z_.]/_/g' <<< "$f")"; done`  
`rename -n 's/(.*).{1}/$1/' *`

The last line removes an underscore, which was introduced when hidden characters were removed (for changes to take effect remove the -n flag).

Transfer the remaining sequences to `01_data`. They will already have been renamed. (This is probably the files from Angela as these already have names following the formula of xxxx_R1.fastq)

Remove files which are of inferior quality  
run `bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/rename_remove.sh`

* * *

## Pipeline
This analysis of the Coryphoideae subfamily utilizes a pipeline / workflow approach using the GWF workflow manager.
The entire pipeline consists of three subsequent scripts that need to be run in the correct order. 
Before commencing the first workflow, make sure that GWF is installed and configured to work on the cluster that you are using.
See https://gwf.app/ for information on how to do this. 


## 01\. species_workflow.py
  Trims the reads from each species using a custom adapter and then combines the paired and unpaired reads for forward and reverse reads respectively in order to enable post trimming fastqc quality check for comparability before and after trimming.

  All the trimmed reads are then given to Hybpiper [link](https://github.com/mossmatters/HybPiper/wiki/) which takes the trimmed reads and builds supercontigs for each target in the target file.  If hybpiper is detecting potential paralogs these are then investigated and removed using the parallel scribt from Hybpiper.

  After Hybpiper all supercontigs for each specimen is collected into a single fasta file. The paired and unpaired reads are then mapped to that fastafile, the reads are then deduplicated and the depth of each base is calculated. Any base with a depth less than 2 are then trimmed and a new fasta file is then created for the supercontigs. 

## 02\. genes_workflow.py
  The genes workflow starts off by distributing all the different supercontigs into different files so that each file contains all the versions of a specific target.

  All of these versions of a supercontig is then aligned using mafft.
  The aligned supercontigs are then trimmed using a list of gap trimming thresholds ranging from 0.1 up to 0.95 in increments of 0.05.

  For each of these trimmed alignments we calculate various summary statistics using AMAS and use these summary statistics to find the gap trimmin threshold which gives us the highest proportion of parsimony informative characters while not removing more data than one median abselute deviation above the median data loss across the entire range of trimming thresholds being tested.

  Additional trimming is then carried out by CIAlign which removes divergent sequences from the multiple sequence alignment and TAPER which removes outlier stretches from each sequence. 

  For each multiple sequence alignment, the two sequences of exons with the highest recovery were added to the multiple sequence alignment.

## 03\. trees_workflow.py
  The two exons added to the alignment are then used to find the best substitution model for the exons and the introns of the multiple sequence alignment and create a RAxML style partition file for the multiple sequence alignment.

  Each of these multiple sequence alignments are then given to IQtree which creates a gene tree for each gene. These genetrees are then all collected in a single newick file. This newick file is then used to create a species tree using ASTRAL.

  This speciestree is then evaluated using the ASTRAl annotation feature.

 

* * *
## Important Notes
***Make sure to download the developmental version of intronerate from Github, as the standard one causes errors when run***

* * *
