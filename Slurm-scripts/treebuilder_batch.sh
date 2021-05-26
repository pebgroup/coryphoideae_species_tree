#!/bin/bash

#SBATCH --job-name=Hybpiper 
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=16
#SBATCH --account=Coryphoideae
#              D-HH:MM:SS
#SBATCH --time=00-02:00:00

#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate treebuilder_env

#Going to folder with data
/home/owrisberg/Coryphoideae/work_flow/10_manual-edit/02_edited_alignments

python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/partitioner.py --smoother 10 

#Genetrees using IQtree
for f in *_clean.fasta
do
	iqtree2 -s $f -T AUTO -ntmax 16 -p ${f/clean.fasta}part.txt -B 1000 # 1000 bootstrap replicates and 16 cores
	mv ${f/clean.fasta}part.txt.treefile /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees/${f/clean.fasta}part.txt.tre
	mv ${f/clean.fasta}part.txt* /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees
	mv ${f/_clean.fasta}.fasta done
	rm $f
done		

cd /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees
rm -f /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree

#Some rerooting, renaming and general cleanup
for f in *.tre
do 
	/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree/genetrees.tre 
	rm temp.tre
done

# Species Tree using Astral
cd /home/owrisberg/Coryphoideae/work_flow/11_tree_building/02_speciestree

rm -f astral*
java -jar /home/owrisberg/Coryphoideae/github_code/ASTRAL/astral.5.7.7.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
 ../rename.csv astral_tree.tre astral_tree_renamed.tre

/home/owrisberg/Coryphoideae/github_code/QuartetScores -o astral_tree_QS.tre -e genetrees.tre -r astral_tree.tre -v
sed astral_tree_QS.tre -i'.old' -e 's/[0-9]\.*[0-9]*\(:[0-9]\.*[0-9]*\)\[qp-ic:-*[0-9]\.[0-9]*;lq-ic:-*[0-9]\.[0-9]*;eqp-ic:\(-*[0-9]\.[0-9]*\)\]/\2\1/g'
sed astral_tree_QS.tre -i'.old' -e 's/\[eqp-ic:-*[0-9]\.*[0-9]*\]//g'

#/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/renamer.py ../rename.csv astral_tree_QS.tre astral_tree_QS_renamed.tre --bs 1