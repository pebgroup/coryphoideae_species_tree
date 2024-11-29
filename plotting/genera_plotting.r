# In this script I will make a figure which shows the relation between all the genera in the Coryphoideae tree

# Load the packages
library(phytools)

# I will start by loading the tree

astral_tree <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree.tre")

# Then We will grap all the tip labels
tibs <- astral_tree$tip.label

# We then split all the tip labels in order to get all the genera
split_labels <- strsplit(tibs, " ")

# Extracting Genera
genera <- sapply(split_labels, `[`, 1)

# Unique genera
unique_genera <- unique(genera)

