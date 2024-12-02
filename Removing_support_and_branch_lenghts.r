# Install and load ape if not already installed
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
library(ape)

# Setting the working directory
setwd("/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/02_speciestree")

# Read the phylogenetic tree from a Newick file
tree <- read.tree("astral_tree.tre")

# Set branch lengths to NULL to remove them
tree$edge.length <- NULL

# Set node labels to NULL to remove node support values
tree$node.label <- NULL


# Write the modified tree to a new Newick file
write.tree(tree, file = "astral_tree_cleaned.tre")

# Now do the same for the orthologs tree
tree_orthologs <- read.tree("astral_tree_orthologs.tre")

# Set Branch lenghts to NULL to remove them
tree_orthologs$edge.length <- NULL

# Set node labels to be NULL to remove them
tree_orthologs$node.label <- NULL

# Write modified tree
write.tree(tree_orthologs, file = "astral_tree_orthologs_cleaned.tre")


# Removing branch lenghts and support for gene trees

# set folder 
setwd("/home/au543206/GenomeDK/coryphoideae/work_flow/10_tree_building/01_genetrees/copy_of_ortholog_gene_trees/without_support")

# Get all files in folder
files <- list.files()

# List of genetrees
genetrees <- files

for (i in 1:length(genetrees)){

  genetree <- genetrees[i]

  # Read the phylogenetic tree from a Newick file
  tree <- read.tree(genetree)

  # Set branch lengths to NULL to remove them
  tree$edge.length <- NULL

  # Set node labels to NULL to remove node support values
  tree$node.label <- NULL

  # Setting the output filename

  file_name <- sub(".tre$", "", genetree)
  file_name_clean <- paste0(file_name, "_cleaned.tre")

  # Write the modified tree to a new Newick file
  write.tree(tree, file = file_name_clean)

}
