#!/bin/bash

#SBATCH --job-name=Treebuilder 
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --account=Coryphoideae
#              D-HH:MM:SS
#SBATCH --time=00-00:30:00

echo "Starting Script"

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate treebuilder_env

#Cleaning folder with data
cd /home/owrisberg/Coryphoideae/work_flow/09_manual_edit/04_alignments_for_trees
rm *

#Copying data from manual alignment folder
cd /home/owrisberg/Coryphoideae/work_flow/09_manual_edit/02_edited_alignments
cp *fasta ../04_alignments_for_trees

#Going to folder with data
cd /home/owrisberg/Coryphoideae/work_flow/09_manual_edit/04_alignments_for_trees

#The partitioner should produce 2 files for each gene
#one file called {gene}_aligned_part.txt which is the partitioning file
#another called {gene}_aligned_clean.fasta which are just the sequences without the exons
python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/partitioner.py --smoother 10 

echo "Beginning IQtree genetree search"

for f in *_part.txt; do (cp $f ${f/_part.txt}_clean.part); done

for f in *_aligned_clean.fasta
do
	if [[$f -a ]] && [["${f/_aligned_clean_fasta}_aligned_part.txt" -a ]]
	then 
		ls *clean.fasta | parallel -j 6 iqtree2 -s {} -T AUTO -ntmax 4 -p {.}.part -B 1000
	else
		echo "Skipping gene" $f
		continue
	fi
done

echo "Done with IQtree"

echo "Beginning cleanup and file transfer"

for f in *_aligned_clean.fasta 
do
	mv ${f/clean.fasta}clean.part.treefile /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees/${f/clean.fasta}part.txt.tre
	mv ${f/clean.fasta}part.txt* /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees #This works
 	mv ${f/_aligned_clean.fasta}.fasta 05_alignments_used_in_trees
 	rm ${f}
done


cd /home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees
rm -f /home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree

echo "Done with cleanup"
echo "Starting rerooting of genetrees"

#Some rerooting, removing of gene in species name and general cleanup
for f in *.tre
do 
	echo $f ${f/_aligned_part.txt.tre}
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/genenameremover.py $f ${f/_aligned_part.txt.tre} #Removes gene-name
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/rooter.py $f #Roots the phylogeny with outgroup
	nw_ed temp.tre 'i & (b<30)' o >> /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree/genetrees.tre #Moves trees used in treebuilding
	rm temp.tre
done

echo "Done with reroot of genetrees"

# Species Tree using Astral
cd /home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree

rm -f astral*