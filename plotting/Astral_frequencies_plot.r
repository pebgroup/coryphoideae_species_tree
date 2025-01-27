library(phytools)

# read tree
tree <- read.tree("/Users/au265104/Library/CloudStorage/OneDrive-Aarhusuniversitet/PROJECTS/95 Coryphoideae phylogeny/gene tree discordance/Wrisberg_et_al_data/05_species_trees/all_genes/astral_tree_annotated.tre")

# set all edges to 1
tree$edge.length <- rep(1, length(tree$edge.length))

# translate tip numbers to tip names 
transl <- read.csv("/Users/au265104/Library/CloudStorage/OneDrive-Aarhusuniversitet/PROJECTS/95 Coryphoideae phylogeny/gene tree discordance/Wrisberg_et_al_data/06_Supporting_Information/names_for_tips.csv", header=F)
transl <- transl[-c(473, 497, 424, 472),] #remove duplicates
rownames(transl) <- transl$V1
tree$tip.label <- transl[tree$tip.label,"V2"]

# tidy up outgroup
tree <- drop.tip(tree, c("Nypa fruticans", "Eugeissona tristis"))

# root tree
tree <- root.phylo(tree, c("Chrysalidocarpus ambositrae", "Aphandra natalia"))

# ladderize tree
tree <- ladderize(tree, right=F)

# extract quartet frequenies (q1, q2, q3) for all nodes
qs <- matrix(nrow=tree$Nnode,ncol=3)

for(i in 1:tree$Nnode){
  qs[i,1] <- as.numeric(substring(strsplit(tree$node.label[i],";")[[1]][1],6))
  qs[i,2]  <- as.numeric(substring(strsplit(tree$node.label[i],";")[[1]][2],4))
  qs[i,3]  <- as.numeric(substring(strsplit(tree$node.label[i],";")[[1]][3],4))
}

# plot tree
pdf("quartet_scores.pdf", width = 8.3, height = 6*11.7)
plot(tree)
nodelabels(pie = qs, piecol=c("blue","darksalmon","grey"), cex=.5)
dev.off()