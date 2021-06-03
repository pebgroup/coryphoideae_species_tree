#!/bin/bash

#SBATCH --job-name=Treebuilder 
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --account=Coryphoideae
#              D-HH:MM:SS
#SBATCH --time=00-03:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate treebuilder_env

#Cleaning folder with data
#cd /home/owrisberg/Coryphoideae/work_flow/10_manual-edit/04_alignments_for_trees
#rm *

#Copying data from manual alignment folder
#cd /home/owrisberg/Coryphoideae/work_flow/10_manual-edit/02_edited_alignments
#cp *fasta ../04_alignments_for_trees

#Going to folder with data
#cd /home/owrisberg/Coryphoideae/work_flow/10_manual-edit/04_alignments_for_trees

#The partitioner should produce 2 files for each gene
#one file called {gene}_aligned_part.txt which is the partitioning file
#another called {gene}_aligned_clean.fasta which are just the sequences without the exons
#python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/partitioner.py --smoother 10 


# echo "Beginning IQtree genetree search"

# for f in *_part.txt; do (cp $f ${f/_part.txt}_clean.part); done
# ls *clean.fasta | parallel -j 6 iqtree2 -s {} -T AUTO -ntmax 4 -p {.}.part -B 1000

# echo "Done with IQtree"

# echo "Beginning cleanup and file transfer"

# for f in *_aligned_clean.fasta 
# do
# 	mv ${f/clean.fasta}clean.part.treefile /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees/${f/clean.fasta}part.txt.tre
# 	mv ${f/clean.fasta}part.txt* /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees #This works
# 	mv ${f/_aligned_clean.fasta}.fasta 05_alignments_used_in_trees
# 	rm ${f}
# done


cd /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees
# rm -f /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree

# echo "Done with cleanup"
echo "Starting rerooting of genetrees"

#Some rerooting, renaming and general cleanup
for f in *.tre
do 
	echo $f ${f/_aligned_part.txt.tre}
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/genenameremover.py $f ${f/_aligned_part.txt.tre}
	python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree/genetrees.tre 
	rm temp.tre
done

echo "Done with reroot of genetrees"
echo "Starting Astral species tree search"

# Species Tree using Astral
cd /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree

rm -f astral*
java -jar /home/owrisberg/Coryphoideae/github_code/ASTRAL/astral.5.7.7.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv astral_tree.tre astral_tree_renamed.tre

/home/owrisberg/Coryphoideae/github_code/QuartetScores -o astral_tree_QS.tre -e genetrees.tre -r astral_tree.tre -v
sed astral_tree_QS.tre -i'.old' -e 's/[0-9]\.*[0-9]*\(:[0-9]\.*[0-9]*\)\[qp-ic:-*[0-9]\.[0-9]*;lq-ic:-*[0-9]\.[0-9]*;eqp-ic:\(-*[0-9]\.[0-9]*\)\]/\2\1/g'
sed astral_tree_QS.tre -i'.old' -e 's/\[eqp-ic:-*[0-9]\.*[0-9]*\]//g'

echo "Done with Astral"
echo "starting renaming"

python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv astral_tree_QS.tre astral_tree_QS_renamed.tre --bs 1

echo "Done with renaming"