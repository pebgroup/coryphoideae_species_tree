# In this script I will make a figure which shows the relation between all the genera in the Coryphoideae tree

# Load the packages
library(phytools)

# I will start by loading the tree
astral_tree <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree.tre")

# Renaming the tips in the tree based on the names for tips 
# read figurename translation table (SECAPR No. to figure name)
figurename <- read.table("/home/au543206/GenomeDK/Coryphoideae/github_code/coryphoideae_species_tree/names_for_tips.csv", sep=",", colClasses = "character")
#figurename <- figurename[1:518,]
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1
figurename_idx

# Renaming tips in tree
astral_tree$tip.label = figurename_idx[astral_tree$tip.label]

# Then We will grap all the tip labels
tips <- astral_tree$tip.label
tips

# We then split all the tip labels in order to get all the genera
split_labels <- strsplit(tibs, " ")
split_labels

# Extracting Genera
genera <- sapply(split_labels, `[`, 1)

# Unique genera
unique_genera <- unique(genera)
unique_genera

# Rename all tips to just their genus name. 
# Modify the tip labels to only include the genus name
tree <- astral_tree
tree$tip.label <- sapply(strsplit(tree$tip.label, " "), `[`, 1)

# For each genera, check if it is monophyletic in the tree
# Function to test monophyly of a genus
test_monophyly <- function(genus, tree, genera) {
  # Get the tips corresponding to the genus
  tips <- which(genera == genus)
  # Check if the tips form a monophyletic group
  is.monophyletic(tree, tips)
}

# Apply the function to all unique genera
results <- sapply(unique_genera, test_monophyly, tree = tree, genera = genera)

# Print the results
cat("Monophyly test results:\n")
print(results)

results[which(results == FALSE)]

# Now for all the genera which are monophyletic, drop all but one of the tips
# Identify monophyletic genera
monophyletic_genera <- unique_genera[sapply(unique_genera, test_monophyly, tree = tree, genera = genera)]
monophyletic_genera

# Drop all but one tip for each monophyletic genus
keep_tips <- unlist(lapply(monophyletic_genera, function(genus) {
  which(genera == genus)[1]  # Keep the first tip of the monophyletic genus
}))

# Add the non monophyletic genera to the keep tips


tree_test <- keep.tip(tree, keep_tips)
plot(tree_test)


# Identify monophyletic and non-monophyletic genera
genera_unique <- unique(genera)
monophyletic_genera <- genera_unique[sapply(genera_unique, test_monophyly, tree = tree, genera = genera)]
non_monophyletic_genera <- setdiff(genera_unique, monophyletic_genera)

# Tips to keep
keep_tips <- c(
  # One representative from each monophyletic genus
  unlist(lapply(monophyletic_genera, function(genus) {
    which(genera == genus)[1]  # Keep the first tip
  })),
  # All tips from non-monophyletic genera
  unlist(lapply(non_monophyletic_genera, function(genus) {
    which(genera == genus)  # Keep all tips from non-monophyletic genera
  }))
)

# I Need to Keep Livistona exigua

tree_test_1 <- keep.tip(tree, keep_tips)
plot(tree_test_1)

# Root the tree using the outgroup
rooted_tree <- root(tree_test_1, outgroup = "Nypa", resolve.root = TRUE)
plot(rooted_tree, show.node.label = TRUE)

tips_for_mrca <- which(tree$tip.label %in% c("Caryota","Wallichia","Arenga"))
tips_for_mrca
mrca <- findMRCA(tree_test_1, tree_test_1$tip.label[tips_for_mrca])
mrca

tree_test_1_collapsed <- collapse.to.star(rooted_tree, 139)

plot(tree_test_1_collapsed) 


tree_test_1$tip.label[tips_for_mrca]

plot(rooted_tree,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:rooted_tree$Nnode,node=1:rooted_tree$Nnode+Ntip(rooted_tree))


par(mar=c(1,1,1,1)) # set margins
plot(rooted_tree,font=1); nodelabels(bg="white")
