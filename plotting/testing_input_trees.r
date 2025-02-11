# This script is meant to test whether all the species found in the genetrees are also found in the final astral tree.
library(phytools)

# Astral tree based on all genes
# We start of by loading the trees
astral_all_gene_tree <- read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree.tre")

# We then load the gene trees
astral_all_gene_genetree <- read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/genetrees.tre")


# Function to get unique tip labels from a multiPhylo object
get_unique_tip_labels <- function(multi_phylo) {
  unique(unlist(lapply(multi_phylo, function(tree) tree$tip.label)))
}

# Get unique tip labels from the trees
unique_labels <- get_unique_tip_labels(astral_all_gene_genetree)
unique_labels

all(astral_all_gene_tree$tip.label %in% unique_labels) # True
all(unique_labels %in% astral_all_gene_tree$tip.label) # True

# Now we do it on all the ortholog genetrees
astral_ortholog_tree <- read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree_orthologs.tre")

# We then load the gene trees
astral_orthologs_genetree <- read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/genetrees_orthologs.tre")

unique_labels_orthologs <- get_unique_tip_labels(astral_orthologs_genetree)

all(astral_ortholog_tree$tip.label %in% unique_labels_orthologs) # True
all(unique_labels_orthologs %in% astral_ortholog_tree$tip.label) # True

all(astral_ortholog_tree$tip.label %in% astral_all_gene_tree$tip.label)
all(astral_all_gene_tree$tip.label %in% astral_ortholog_tree$tip.label) # False

# Function to find species in one tree but not in another
find_missing_species <- function(tree1, tree2) {
  setdiff(tree1$tip.label, tree2$tip.label)
}

# Example usage (replace with actual tree objects)
missing_species <- find_missing_species(astral_all_gene_tree, astral_ortholog_tree)
print(missing_species)

# 4021,Coccothrinax yuraguana Bro. León 18468
# 4100,Pritchardia martii
# 3326,Saribus chocolatinus
# 3224,Licuala ruthiae


# I also need to test the folder which contains all the "raw" genetrees
file_folder <- "/home/au543206/GenomeDK/Coryphoideae/work_flow/10_tree_building/01_genetrees/" # Folder on GenomeDK

# Function to load multiple .tre files into a multiPhylo object
load_trees <- function(folder, pattern = "_rooted.tre$") {
  # List all files in the folder that match the pattern
  tree_files <- list.files(folder, pattern = pattern, full.names = TRUE)
  
  # Read all tree files into a list
  trees <- lapply(tree_files, read.tree)
  
  # Convert list to multiPhylo object
  multi_trees <- do.call(c, trees)
  class(multi_trees) <- "multiPhylo"
  
  return(multi_trees)
}

# Finding all genetrees
genetrees_raw <- load_trees(file_folder, "_rooted.tre$")

# Finding the unique labels from the raw genetrees.
genetrees_raw_unique_labels <- get_unique_tip_labels(genetrees_raw)

# Are all the unique labels from the genetrees found in the species tree
all(genetrees_raw_unique_labels %in% astral_all_gene_tree$tip.label) # TRUE

# Are all the unique labels from the raw genetrees foun in the genetrees that are concatenated
all(genetrees_raw_unique_labels %in% unique_labels) # TRUE

# Which are not found in the 
which(genetrees_raw_unique_labels %in% astral_all_gene_genetree)


