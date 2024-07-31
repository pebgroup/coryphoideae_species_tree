# Coryphoideae targeted capture bioinformatics work-flow

Inherited by Oscar Wrisberg ([Oscar.wrisberg@bio.au.dk](mailto:Oscar.wrisberg@bio.au.dk)), 11th January 2021

* * *

## 0\. Workspace

**Required file structure and environments**  
In order to run this pipeline you need a directory with the following folders along with a set of conda environments.
These can be quickly produced by running `infrastructure.sh` within the desired directory.
**OBS!!** one of the environments is called **base** and running the infrastructure script might **erase** your own base environment, proceed at your own risk.

- `00_secapr`
  - `0_data`
  - `1_trimmed`
- `01_data`
- `02_trimmed`
- `03_hybpiper`
- `04_coverage`
- `05_blacklisting`
- `06_alignment`
- `07_optrimal`
- `08_cialign`
- `09_mapping`
- `10_tree_building`
- `11_dating_the_tree`

## Downloading and renaming data

Transfer all sequence files into `01_data` from where the analysis begins.  
In order for the pipeline to run, the files need to follow a strict naming regime, see [names_number_species.txt](./names_number_species.txt)

* * *

## Pipeline

This analysis of the Coryphoideae subfamily utilizes a pipeline / workflow approach using the GWF workflow manager.
The entire pipeline consists of three subsequent scripts that need to be run in the correct order.
Before commencing the first workflow, make sure that GWF is installed and configured to work on the cluster that you are using.
See [GWF](https://gwf.app/) for information on how to do this.

## 01\. species_workflow.py

  Trims the reads from each species using [trimmomatic](https://github.com/usadellab/Trimmomatic) and a custom adapter and then combines the paired and unpaired reads for forward and reverse reads respectively in order to enable post trimming [fastqc](https://github.com/s-andrews/FastQC) quality check for comparability before and after trimming.

  All the trimmed reads are then given to [Hybpiper](https://github.com/mossmatters/HybPiper/wiki/) which takes the trimmed reads and builds supercontigs for each target in the target file.  If hybpiper is detecting potential paralogs these are then investigated and removed using the parallel scribt from Hybpiper.

  After Hybpiper all supercontigs for each specimen is collected into a single fasta file. The paired and unpaired reads are then mapped to that fastafile, the reads are then deduplicated and the depth of each base is calculated. Any base with a depth less than 2 are then trimmed and a new fasta file is then created for the supercontigs.

## 02\. genes_workflow.py

  The genes workflow starts off by distributing all the different supercontigs into different files so that each file contains all the versions of a specific target.

  All of these versions of a supercontig is then aligned using [MAFFT](https://mafft.cbrc.jp/alignment/software/).
  The aligned supercontigs are then trimmed using a list of gap trimming thresholds ranging from 0.1 up to 0.95 in increments of 0.05 using [trimAl](http://trimal.cgenomics.org/)

  For each of these trimmed alignments we calculate various summary statistics using [AMAS](https://github.com/marekborowiec/AMAS) and use these summary statistics to find the gap trimmin threshold which gives us the highest proportion of parsimony informative characters while not removing more data than one median abselute deviation above the median data loss across the entire range of trimming thresholds being tested. This step is done using the R-script [optrimal](https://github.com/baileyp1/PhylogenomicsPipelines)

  Additional trimming is then carried out by [CIAlign](https://github.com/KatyBrown/CIAlign) which removes divergent sequences from the multiple sequence alignment and [TAPER](https://github.com/chaoszhang/TAPER) which removes outlier stretches from each sequence.

  For each multiple sequence alignment, the two sequences of exons with the highest recovery were added to the multiple sequence alignment.

## 03\. trees_workflow.py

  The two exons added to the alignment are then used to find the best substitution model for the exons and the introns of the multiple sequence alignment and create a RAxML style partition file for the multiple sequence alignment.

  Each of these multiple sequence alignments are then given to [IQtree](http://www.iqtree.org/) which creates a gene tree for each gene. These genetrees are then all collected in a single newick file. This newick file is then used to create a species tree using [ASTRAL](https://github.com/smirarab/ASTRAL).

  This speciestree is then evaluated using the [ASTRAl](https://github.com/smirarab/ASTRAL) annotation feature.

* * *

## 04\. Final scripts
  There are some final scripts which produce the figures and stats provided in the article.
  These scripts are small enough to not require a cluster/workflow manager and are therefore not included in the pipeline.
  Run these scripts on your local machine with the output from the pipeline as input.
  *OBS* None of these scripts have relative paths and you need to check that you have the correct paths to the input files. 

  The scripts are:
   - 'load_trees.r' which loads all the output files from the pipeline and produces a list of all the trees.
   - 'functions.r' contains all the custom plotting functions used in main_figure.r
   - 'main_figure.r' which produces the main figure of the article along with the supplementary figures.
   - 'Finding_wcvp_coryphoideae.r' which produces a subset of the World checklist of vascular plants for the Coryphoideae subfamily.
   - 'counting_genera_coverage.r' which produces the coverage of the different genera in the dataset.

## Important Notes

***Make sure to download the developmental version of intronerate from Github, as the standard one causes errors when run.***

* * *
