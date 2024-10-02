# In this paper I want to plot the coryphoid phylogeny using the daylight method as a beautiful figure for my thesis

library(ggtree)
library(phytools)

# Load the data
setwd("/home/au543206/Documents/Coryphoideae/Figures_in_r/data")

# Load the tree
tree <- read.tree("astral_tree_annotated_EN.tre")

# Dropping the outgroup
outgroup <- c("1079","1080","1081","1082","1002","1003")

tree <- drop.tip(tree, outgroup)

# plot the tree
daylight <- ggtree(tree, layout="daylight", branch.length = 'none')


circular <- ggtree(tree, branch.length='none', layout='circular', size = 1.2)
circular

circular_tip_labels <- ggtree(tree, branch.length='none', layout='circular', size = 0.9) + geom_tiplab(aes(angle = angle), size = 2, hjust = 0, vjust = 0.5)
circular_tip_labels

tree$tip.label
