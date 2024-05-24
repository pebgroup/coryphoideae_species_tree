# This script is used to count the number of species, subspecies and variants in the final tree and the ortholog tree.
# This is done using the latest version of the World Checklist of Vascular Plants (WCVP).

# Loading some command line elements
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("ape","phytools","data.table")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Loading the functions
repo <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/Cory_tree_gene_support"
source(paste(repo,"/load_trees.R", sep=""))

# Loading the arguments manually
wcvp <- fread("/home/au543206/Documents/world_checklist/World_checklist_downloads/05_09_2024/wcvp_names.csv") # here you define the name of your input wcvp file


astral_tree_for_figure_orthologs <- astral_tree_orthologs_EN
#astral_tree_for_figure_orthologs <- drop.tip(astral_tree_for_figure_orthologs,"1362")
astral_tree_for_figure_orthologs$tip.label = figurename_idx[astral_tree_for_figure_orthologs$tip.label]
astral_tree_for_figure_orthologs$edge.length[which(astral_tree_for_figure_orthologs$edge.length == "NaN")] <- 0.001

astral_tree_for_figure_orthologs$tip.label

# split all the species names so i can get the genera
unlist_sp_names <- unlist(strsplit(astral_tree_for_figure_orthologs$tip.label, " "))
unlist_sp_names
