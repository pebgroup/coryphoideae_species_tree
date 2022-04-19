#!/bin/bash

#SBATCH --account Coryphoideae
#SBATCH --job-name=Intronerate
#SBATCH --partition normal
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=24
#              D-HH:MM:SS
#SBATCH --time=0-10:00:00

# This File runs the sortadate functions on the gene trees and gives you a list of stats which can be used to determine which genes should be used in dating
#


#This line should enable the activation of specific conda environments
source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#Activating conda base environment 
conda activate base

#Get the root-to-tip variance with
python /home/owrisberg/Coryphoideae/github_code/SortaDate/src/get_var_length.py /home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/ --flend _part.txt.tre --outf /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/var --outg 1079,1080,1081,1082

#Get the bipartition support with
#What kind of values is the values shown in the example file Chrono_Tent_Bird_study.new, and do i need to care?
python /home/owrisberg/Coryphoideae/github_code/SortaDate/src/get_bp_genetrees.py /home/owrisberg/Coryphoideae/work_flow/10_tree_building/01_genetrees/ /home/owrisberg/Coryphoideae/work_flow/10_tree_building/02_speciestree/astral_treee.tre --flend _part.txt.tre --outf /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/bp

#Combine the results from these two runs with
python /home/owrisberg/Coryphoideae/github_code/SortaDate/src/combine_results.py /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/var /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/bp --outf /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/comb

#Sort and get the list of the good genes with
python /home/owrisberg/Coryphoideae/github_code/SortaDate/src/get_good_genes.py /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/comb --max 3 --order 3,1,2 --outf /home/owrisberg/Coryphoideae/work_flow/11_dating_the_tree/00_sortadate/gg






 # Outgroup species
 #if "1079" in tips and "1080" in tips and "1081" in tips and "1082" in tips:


 # Sortadate example functions
 # python src/get_var_length.py examples/genes_trees/ --flend .tre.rr --outf examples/var --outg Struthio_camelus,Tinamou_guttatus

 # python src/get_bp_genetrees.py examples/genes_trees/ examples/Chrono_Tent_Bird_study.new --flend .tre.rr --outf examples/bp

 # python src/combine_results.py examples/var examples/bp --outf examples/comb

 #  python src/get_good_genes.py examples/comb --max 3 --order 3,1,2 --outf examples/gg