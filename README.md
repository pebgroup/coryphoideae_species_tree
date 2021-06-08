# Coryphoideae targeted capture bioinformatics work-flow

Inherited by Oscar Wrisberg ([Oscar.wrisberg@bio.au.dk](mailto:Oscar.wrisberg@bio.au.dk)), 11th January 2021

* * *

## 0\. Workspace

All raw sequences are found on the Ecoinf drive:  
`/home/au543206/Ecoinf/C_Write/_Proj/@PEB/lab/sequence_data`  

Copy of raw sequences on GenomeDK are found at:  
`/home/owrisberg/Coryphoideae/sequence_data/Coryphoideae`

GitHub project repo clone on GenomeDK:  
`/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree`

GitHub HybPiper repo clone on GenomeDKr:  
`/home/owrisberg/Coryphoideae/github_code/HybPiper`

Directory containing target file on GenomeDK:  
`/home/owrisberg/Coryphoideae/target_sequence`

**Required file structure**  
In order to run this pipeline you need a directory with the following folders, this can be quickly produced by running `infrastructure.sh` within the desired directory.

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

## 01\. data

SECAPR quality check is run on the raw data in `01_data`.  
In order to invoke `secapr` first activate the environment  
`conda activate secapr_env`  
Then run SECAPR from within the directory  
`secapr quality_check --input . --output .`  
Move the secapr files to `0_secapr/0_data`

A list of the fastq files is created for the directory `01_data` by running the following script from within it  
`ls *R1.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > namelist.txt; rm namelist_temp.txt`  
The second command removes the last 9 characters. This list is needed for running Trimmomatics on all the names in the list.

* * *

## 02\. Trimmed

In order to trim the files, `run trimming_batch.sh`.
**OBS** remember to check if the TruSeq3-PE-2.fa file is uploaded to the folder located at `/home/owrisberg/miniconda3/pkgs/trimmomatic-0.39-1/adapters/`

In order to produce the quality check on the trimmed files run the `secapr_batch_trimmed.sh`.  

* * *

## 03\. hybpiper

### Combine unpaired reads into a single file for each sample

Run the following script within `02_trimmed`  
`/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/comb_u_trim_reads.sh`  
This merges `####_1U.fastq` and `####_2U.fastq` into `####_UN.fastq`

### Generate namelist

A list of the fastq files within the directory `02_trimmed` is created in the directory `03_hybpiper` by running the following script from within `02_trimmed`  
`ls *1P.fastq > namelist_temp.txt; sed 's/.........$//' namelist_temp.txt > ../3_hybpiper/namelist.txt; rm namelist_temp.txt`  
The second command removes the last 9 characters. This list is needed for running hybpiper on all the listed names.  
If all trimmed data are to go into hybpiper the namelist within `01_data` can alternatively be copied to `03_hybpiper` and reused.

### Execute HybPiper

From within `03_hybpiper` run  
`bash /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/hybpiper_gdk.sh`

### Get assembly stats

From within `03_hybpiper` run  
`python /home/owrisberg/Coryphoideae/github_code/HybPiper/get_seq_lengths.py /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist.txt dna > seq_lengths.txt`  
and  
`python ~/Coryphoideae/github_code/HybPiper/hybpiper_stats.py seq_lengths.txt namelist.txt > stats.txt`  
The file `seq_lengths.txt` can be used to make a heatmap in R by running the script `gene_recovery_heatmap_ggplot.R`

### Paralogs

From within `03_hybpiper` run  
`bash ~/Coryphoideae/github_code/coryphoideae_species_tree/paralog2.sh`  
The output from this program is saved in the file `paralog.txt`. This file lists putatively paralogous genes. In order to have the name of every gene only listed once only run the following  
`sort paralog.txt | uniq | sed 's/^.......................//' > gene_paralogs.txt`

### Intronerate

In order to generate the super contigs we need to run intronerate.
Run the `intronerate_batch.sh`
***Make sure to download the developmental version of intronerate from Github, as the standard one causes errors when run***

* * *

## 04\. Coverage trimming and Length filtering

Run the coverage program by running the `coverage_batch.sh`.

Ensure that "supercontig" is chosen in the coverage.py script. This is currently done by commenting two lines of code.

The Coverage.py script does the following:

- Gather all contigs from each sample in one fasta file
- map paired and unpaired reads to that fasta using BWA mem
- Deduplicate reads using Picard
- Calculate depth using samtools
- Mask/strip any bases with coverage less than 2
- Generate a new trimmed sample-level fasta.

Then, in `04_coverage`, run:

`ls *trimmed.fasta > filelist.txt`
`/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/samples2genes.py > outstats.csv`

* * *

## 05\. Blacklisting

This step is kept for removing troublesome species.
These species usually show up further downstream in the analysis or as a part of the tree building process.
In order to remove the troublesome species run the `blacklisting_batch.sh`. ***OBS species are hardcoded into the script***

***

# From here on out, everything needs to be done together

## 06\. Alignment

From within `05_blacklisting` run MAFFT on all genes  
`mafft_batch.sh`
Aligned genes are found within `06_alignments`

***OBS! There are currently several mafft_batch.sh scripts in the folder as the different genes were run in parallel due to urgency***

Visualize single gene alignments with AliView. Launch the program with command: `aliview`  
Visualize multiple gene alignments with Geneious. Launch from Nautilus.

* * *

## 07\. mapping

Run the `exon_mapper_batch` script.
This creates new alignments in `07_mapping` that contain the original alignments plus the exon sequences of the two species that had the highest recovery success at each locus.

then run `cp *.fasta ../08_optrimal` in order to copy the alignments to the optrimal folder.

* * *

## 08\. Optrimal

Create a file called cutoff_trim.txt with the -gt values which should be tested.

Run the `gaptrimming_batch.sh` script.

This will a folder for each of the -gt values.
In each of these folders there will be created a .fasta file, a .htm file and a noempty.fasta file for each gene.
in addition a summary file for the trimal process will also be created in the folder.

* * *

## 09\. Manual editing

Move all the *aligned.fasta files to folders in manual alignment using the command
`mv *aligned.fasta /home/owrisberg/Coryphoideae/work_flow/09_manual-edit/01_alignments_for_editing`

Manually edit sequences to ensure proper alignment.
When manually editing an alignment, move it to `02_edited_alignments` and edit it.
* * *

## 10\. Tree building

Run the `treebuilder_batch.sh` script.

The treebuilder script it will copy the files from `09_manual_edit/02_edited_alignments` into `09_manual_edit/04_alignments_for_trees` and perform the following on these alignments.

This batch file will first run the `partitioner.py` with a smoothing parameter of 10bp (i.e. ignoring any mini-partitions <10bp long) to generate RAxML-style partition files called *_part.txt, and remove the exon sequences from the alignment (new alignment file saved as*_clean.fasta)
If the exons of a specific gene are of unequal length the gene is added to a text file called badgenes.text. remember to have a look through this file.

it will then run IQtree on each gene within the directory, and add the genetrees to the genetrees.tre file in the `/home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree` folder.

When the script is done running you can run `java -jar /home/owrisberg/Coryphoideae/github_code/ASTRAL/astral.5.7.7.jar -i genetrees.tre -o astral_tree.tre  2> astral.log` in order to produce the species tree from the genetrees.

You can rename the astral tree by running.
`python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv astral_tree.tre astral_tree_renamed.tre`

When astral is done, the astral tree needs to be evaluated using QuartetScores
`/home/owrisberg/Coryphoideae/github_code/QuartetScores -o astral_tree_QS.tre -e genetrees.tre -r astral_tree.tre -v`

Run these commands in order to apply the correct labels on the Quartet scores
`sed astral_tree_QS.tre -i'.old' -e s/[0-9]\.*[0-9]*\(:[0-9]\.*[0-9]*\)\[qp-ic:-*[0-9]\.[0-9]*;lq-ic:-*[0-9]\.[0-9]*;eqp-ic:\(-*[0-9]\.[0-9]*\)\]/\2\1/g``
`sed astral_tree_QS.tre -i'.old' -e 's/\[eqp-ic:-*[0-9]\.*[0-9]*\]//g``

Finally in order to rename the tips in the Quartet_scored tree run
`python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv astral_tree_QS.tre astral_tree_QS_renamed.tre --bs 1`
* * *
