wd <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/Cory_tree_gene_support"
setwd(wd)

figurepath <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/"


library(tidytree)
library(ape)
library(phytools)
library(classInt)
library(colorspace)
library(grDevices)
library(stringr)
library(adephylo)
library(MetBrewer)



##### Testing
# astral_tree$tip.label[duplicated(astral_tree$tip.label)]
# astral_tree$tip.label
# which(duplicated(astral_tree$tip.label)== TRUE)
# which(duplicated(astral_tree_EN$tip.label)==TRUE)
# which(duplicated(astral_tree_for_figure$tip.label)==TRUE)

repo <- "/home/au543206/Documents/Coryphoideae/coryphoideae_species_tree/plotting/"
source(paste(repo,"/functions.r", sep=""))
source(paste(repo,"/load_trees.r", sep=""))

astral_tree <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree.tre")
astral_tree <- drop.tip(astral_tree_orthologs, "4039")
#astral_tree <- ape::read.nexus("/home/au543206/Documents/Coryphoideae/Figures_in_r/astral_tree_renamed.nexus")
#astral_tree <- root(astral_tree, outgroup=c("1081"))
astral_tree <- ladderize(astral_tree, right=FALSE)

astral_tree$edge.length[which(astral_tree$edge.length == "NaN")] <- 0.00001


genetrees <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/genetrees.tre")

# create list of representation across gene trees: tipname, locusname
locusname <- c()
for(i in 1:length(gtnames)){ # gtnames is from load_tree.R
  locusname <- c(locusname, strsplit(gtnames[i],"_")[[1]][1])
}
rm(i)

data = matrix(nrow=0,ncol=2)
colnames(data) <- c("tipname","locusname")

#This for loop creates a dataframe which contains all the tipnames in each genetree.
#This takes a really long time to run through 970 genes X 400 species = 388.000 rows
#The dataframe is therefore saved as a RDS at the end for easier loading afterwards

# for(i in 1:length(gts)){
#   #i = 1
#   for(tip in gts[[i]]$tip.label){
#     data <- rbind(data,c(tip,locusname[i]))
#     print(locusname[i])
#   }
# }

# data <- as.data.frame(data)
# saveRDS(data,file="gts_data")

data <- readRDS("gts_data")

# Get support values as edge labels
astral_tree$node.label -> supportvals
names(supportvals) <- 1:length(astral_tree$node.label) + Ntip(astral_tree)

# Translate node labels into edge labels
edgelabs <- c()
for(i in 1:nrow(astral_tree$edge)){ # loops through all edges
  dscndnt <- as.character(astral_tree$edge[i,2]) # selects the crown node of each edge
  if(dscndnt %in% names(supportvals)){ # if this is an internal branch.... 
    edgelabs <- c(edgelabs, supportvals[dscndnt]) # assign the support value of the crown node
  } else {
    edgelabs <- c(edgelabs, "") # if not, assign nothing. 
  }
}

edgelabs[edgelabs=="1"] <- "" # suppress support values == 1
# reformat numbers (no leading zero, two decimals)
edgelabs <- rapply(as.list(as.numeric(edgelabs)), sprintf, fmt = "%0.2f", how = "replace") 
edgelabs <- unlist(edgelabs)
edgelabs[edgelabs == "NA"] <- ""
edgelabs <- str_replace(edgelabs, "0.", ".")

# Get edge colours

#Rooting tree on outgroup species. 
# 1081, Eugeissona tristis
# 1082, Nypa fructicans

astral_tree_EN <- root(astral_tree_EN, outgroup=c( "1081", "1082"))
astral_tree_EN <- drop.tip(astral_tree_EN, "4039")
astral_tree_EN <- drop.tip(astral_tree_EN, "4091")
astral_tree_EN <- ladderize(astral_tree_EN, right=F)
astral_tree_EN$edge.length[is.na(astral_tree_EN$edge.length)] <- 0.001

EN <- as.numeric(astral_tree_EN$node.label)
names(EN) <- 1:length(astral_tree_EN$node.label) + Ntip(astral_tree_EN)

nodeclass <- classIntervals(EN[!is.na(EN)], 6, style = "jenks")
pal <- c(met.brewer("Peru1",6, type = "discrete"))
as.vector(findColours(nodeclass, pal)) -> ENcolour

names(ENcolour) <- names(EN[!is.na(EN)])

edgecols <- rep("black",nrow(astral_tree_EN$edge))

# loop through internal nodes and colour branches
for(nl in (1:length(astral_tree_EN$node.label) + Ntip(astral_tree_EN))){
  if(nl %in% names(ENcolour)){
    print(nl)
    edgecols[astral_tree_EN$edge[,2]==nl] <- ENcolour[as.character(nl)]
  }
}

tipres <- table(data$tipname)
tipclass <- nodeclass
tipclass$var <- tipres
tipclass$brks[7] <- max(tipres)
as.vector(findColours(tipclass, pal)) -> tipcols
names(tipcols) <- names(tipres)

# loop through terminal nodes and colour branches
i = 1
for(nl in astral_tree_EN$tip.label){
  if(nl %in% names(tipcols)){
    print(tipcols[nl])
    edgecols[astral_tree_EN$edge[,2]==i] <- tipcols[nl]
    i = i +1
  }
}


tipcols <- tipcols[astral_tree_EN$tip.label]

# Creating tree for Figure
astral_tree_for_figure <- astral_tree_EN

# Renaming tips in tree
astral_tree_for_figure$tip.label = figurename_idx[astral_tree_for_figure$tip.label]


#pdf(paste(figurepath,"main_tree_narrow.pdf",sep=""), height=15.3, width=8.25)
# #pdf(paste(figurepath,"main_tree_broad.pdf",sep=""), height=15.3, width=11)
# pdf(paste(figurepath, "Coryphoideae_tree_genetree_support_final.pdf", sep = ""), height = (15.3*5), width = 27)
# plotwe(astral_tree_for_figure,
#  direction = "rightwards",
#  cex = 0.9,
#  align.tip.label = T,
#  edge.color = edgecols,
#  link.color = tipcols,
#  label.offset = 0.08,
#  x.lim = c(-1.2,max(distRoot(astral_tree_for_figure))+3),
#  edge.width = 2.5),
# #nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
#  edgelabels(text = edgelabs, frame="none", cex=1, adj=c(0.5,-0.30)),
#  legend("topleft", legend=c("0-86","87-269","270-506","507-735","736-882","883-960"),fill=pal,title="No. gene trees")
#  add.scale.bar(x=0, y=-10)
# dev.off()

# pdf(paste(figurepath, "Coryphoideae_tree_ortholog_genetrees_support_leftwards.pdf", sep = ""), height = (15.3*5), width = 27) # height was 15.3 * 3 and width = 11
# plotwe(astral_tree_for_figure_orthologs,
# direction = "rightwards",
# cex = 0.9,
# align.tip.label = T,
# edge.color = edgecols_orthologs,
# link.color = tipcols_orthologs,
# label.offset = 0.08,
# x.lim = c(-1.2,max(distRoot(astral_tree_for_figure_orthologs))+3),
# edge.width = 2.5) # was 1.2
# #nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
# edgelabels(text = edgelabs_orthologs, frame="none", cex=1, adj=c(0.5,-0.30)) # cex was 0.45
# legend("topleft", legend=c("0-34","35-79","80-129","139-177","178-214","215-231"),fill=pal_orthologs,title="No. gene trees")
# add.scale.bar(x=0, y=-10)
# dev.off()




# pdf(paste(figurepath, "Coryphoideae_tree_genetree_support_big_branch.pdf", sep = ""), height = (15.3*3), width = 11)
# plotwe(astral_tree_for_figure, cex = 0.35, align.tip.label = T, edge.color = edgecols, link.color = tipcols, label.offset = 0.04, x.lim = c(0,max(distRoot(astral_tree_for_figure))+1), edge.width = 3)
# #nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
# edgelabels(text = edgelabs, frame="none", cex=0.45, adj=c(0.5,-0.30))
# legend("topleft", legend=c("0-86","87-269","270-506","507-735","736-882","883-960"),fill=pal,title="No. gene trees")
# add.scale.bar(x=0, y=20)
# dev.off()


# nodeclass$brks
# WOLF Data
# 274        215        298        233        236        329        249 
# 1.992441  24.974194  65.057874  91.352708 120.897633 141.353996 160.930571 
# c("1-24","25-65","66-91","92-120","121-141","142-161")

# Oscar Data
# 646       668       544       507       562       810       524 
# 0.00000  86.86602 269.43718 506.52628 735.51260 882.38263 960.89080 
# c("0-86","87-269","270-506","507-735","736-882","883-960")


#######################################################################
############################ Orthologs#################################
#######################################################################

astral_tree_orthologs <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree_orthologs.tre")
astral_tree_orthologs <- drop.tip(astral_tree_orthologs, "4039")

astral_tree_orthologs <- root(astral_tree_orthologs, outgroup=c("1081"))
astral_tree_orthologs <- ladderize(astral_tree_orthologs, right=TRUE)

astral_tree_orthologs$edge.length
astral_tree_orthologs$edge.length[which(astral_tree_orthologs$edge.length == "NaN")] <- 0.00001


genetrees <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/genetrees_orthologs.tre")

# create list of representation across gene trees: tipname, locusname
locusname <- c()
for(i in 1:length(gtnames)){
  locusname <- c(locusname, strsplit(gtnames[i],"_")[[1]][1])
}
rm(i)

##########################################################################################
# data_orthologs = matrix(nrow=0,ncol=2)
# colnames(data_orthologs) <- c("tipname","locusname")
# 
# # This for loop creates a dataframe which contains all the tipnames in each genetree.
# # This takes a really long time to run through 300 genes X 400 species = 120.000 rows
# # The dataframe is therefore saved as a RDS at the end for easier loading afterwards
# for(i in 1:length(gts_orthologs)){
#   for(tip in gts_orthologs[[i]]$tip.label){
#     data_orthologs <- rbind(data_orthologs,c(tip,locusname[i]))
#     print(locusname[i])
#   }
# }
# data_orthologs <- as.data.frame(data_orthologs)
# saveRDS(data_orthologs,file="gts_orthologs_data.Rds")
#########################################################################################

data_orthologs <- readRDS("gts_orthologs_data.Rds")
data_orthologs <- as.data.frame(data_orthologs)

# Get support values as edge labels
astral_tree_orthologs$node.label -> supportvals_orthologs
names(supportvals_orthologs) <- 1:length(astral_tree_orthologs$node.label) + Ntip(astral_tree_orthologs)

# Translate node labels into edge labels
edgelabs_orthologs <- c()
for(i in 1:nrow(astral_tree_orthologs$edge)){ # loops through all edges
  dscndnt <- as.character(astral_tree_orthologs$edge[i,2]) # selects the crown node of each edge
  if(dscndnt %in% names(supportvals_orthologs)){ # if this is an internal branch.... 
    edgelabs_orthologs <- c(edgelabs_orthologs, supportvals_orthologs[dscndnt]) # assign the support value of the crown node
  } else {
    edgelabs_orthologs <- c(edgelabs_orthologs, "") # if not, assign nothing. 
  }
}

edgelabs_orthologs[edgelabs_orthologs=="1"] <- "" # suppress support values == 1

# reformat numbers (no leading zero, two decimals)
edgelabs_orthologs <- rapply(as.list(as.numeric(edgelabs_orthologs)), sprintf, fmt = "%0.2f", how = "replace") 
edgelabs_orthologs <- unlist(edgelabs_orthologs)
edgelabs_orthologs[edgelabs_orthologs == "NA"] <- ""
edgelabs_orthologs <- str_replace(edgelabs_orthologs, "0.", ".")

# Get edge colours

#Rooting tree on outgroup species. 
# 1081, Eugeissona tristis
# 1082, Nypa fructicans

astral_tree_orthologs_EN <- root(astral_tree_orthologs_EN, outgroup=c( "1081", "1082"))
astral_tree_orthologs_EN <- ladderize(astral_tree_orthologs_EN, right=FALSE)
astral_tree_orthologs_EN$edge.length[is.na(astral_tree_orthologs_EN$edge.length)] <- 0.001

EN_orthologs <- as.numeric(astral_tree_orthologs_EN$node.label)
names(EN_orthologs) <- 1:length(astral_tree_orthologs_EN$node.label) + Ntip(astral_tree_orthologs_EN)

# Creating jenks of intervals for the number of genes
nodeclass_orthologs <- classIntervals(EN_orthologs[!is.na(EN_orthologs)], 6, style = "jenks") 
pal_orthologs <- c(met.brewer("Peru1",6, type = "discrete"))
as.vector(findColours(nodeclass_orthologs, pal_orthologs)) -> ENcolour_orthologs

names(ENcolour_orthologs) <- names(EN_orthologs[!is.na(EN_orthologs)]) 
edgecols_orthologs <- rep("black",nrow(astral_tree_orthologs_EN$edge))

# loop through internal nodes and colour branches
# we do this by looping through the internal nodes and grapping the number of genes used for each node.
# We then check if the node is in 
for(nl_orthologs in (1:length(astral_tree_orthologs_EN$node.label) + Ntip(astral_tree_orthologs_EN))){ # looping through all internal nodes
  print(nl_orthologs)
  if(nl_orthologs %in% names(ENcolour_orthologs)){
    edgecols_orthologs[astral_tree_orthologs_EN$edge[,2]==nl_orthologs] <- ENcolour_orthologs[as.character(nl_orthologs)]
    print(edgecols_orthologs[astral_tree_orthologs_EN$edge[,2]==nl_orthologs])
  }
}

tipres_orthologs <- table(data_orthologs$tipname) # Nr of genes found for each specimen
tipclass_orthologs <- nodeclass_orthologs
tipclass_orthologs$var <- tipres_orthologs
tipclass_orthologs$brks[7] <- max(tipres_orthologs)
as.vector(findColours(tipclass_orthologs, pal_orthologs)) -> tipcols_orthologs
names(tipcols_orthologs) <- names(tipres_orthologs)

# loop through terminal nodes and colour branches
i = 1
for(nl_orthologs in astral_tree_orthologs_EN$tip.label){
  print(nl_orthologs)
  if(nl_orthologs %in% names(tipcols_orthologs)){
    print(tipcols_orthologs[nl_orthologs])
    edgecols_orthologs[astral_tree_orthologs_EN$edge[,2]==i] <- tipcols_orthologs[nl_orthologs]
    i = i +1
  }
}

tipcols_orthologs <- tipcols_orthologs[astral_tree_orthologs_EN$tip.label]

astral_tree_for_figure_orthologs <- astral_tree_orthologs_EN
#astral_tree_for_figure_orthologs <- drop.tip(astral_tree_for_figure_orthologs,"1362")
astral_tree_for_figure_orthologs$tip.label = figurename_idx[astral_tree_for_figure_orthologs$tip.label]
astral_tree_for_figure_orthologs$edge.length[which(astral_tree_for_figure_orthologs$edge.length == "NaN")] <- 0.001
#pdf(paste(figurepath,"main_tree_narrow.pdf",sep=""), height=15.3, width=8.25)
#pdf(paste(figurepath,"main_tree_broad.pdf",sep=""), height=15.3, width=11)
pdf(paste(figurepath, "Coryphoideae_tree_ortholog_genetrees_support_leftwards.pdf", sep = ""), height = (15.3*5), width = 27) # height was 15.3 * 3 and width = 11
plotwe(astral_tree_for_figure_orthologs, direction = "rightwards", cex = 0.8, align.tip.label = T,
       edge.color = edgecols_orthologs, link.color = tipcols_orthologs, label.offset = 0.08,
       x.lim = c(-1.2,max(distRoot(astral_tree_for_figure_orthologs))+1.5), edge.width = 2.5) # was 1.2
#nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
edgelabels(text = edgelabs_orthologs, frame="none", cex=1, adj=c(0.5,-0.30)) # cex was 0.45
legend("topleft", legend=c("0-34","35-79","80-129","139-177","178-214","215-231"),fill=pal_orthologs,title="No. gene trees")
#add.scale.bar(x=0, y=-1)
dev.off()


###
astral_tree_for_figure_orthologs_uniform <- astral_tree_for_figure_orthologs
astral_tree_for_figure_orthologs_uniform$edge.length <- rep(2,959)
pdf(paste(figurepath, "Coryphoideae_tree_ortholog_genetrees_support_leftwards.pdf", sep = ""), height = (15.3*5), width = 27) # height was 15.3 * 3 and width = 11
plotwe(astral_tree_for_figure_orthologs_uniform, direction = "rightwards", cex = 0.8, align.tip.label = T,
       edge.color = edgecols_orthologs, link.color = tipcols_orthologs, label.offset = 0.08, edge.width = 5, use.edge.length = TRUE) # was 1.2
#nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
edgelabels(text = edgelabs_orthologs, frame="none", cex=1, adj=c(0.5,-0.30)) # cex was 0.45
legend("topleft", legend=c("0-34","35-79","80-129","139-177","178-214","215-231"),fill=pal_orthologs,title="No. gene trees")
#add.scale.bar(x=0, y=-1)
dev.off()

#########################################################################################
######################-- Creating sub plots of the phylogeny--###########################
#########################################################################################

# Finding the MRCA node for the following species pairs
# 5 sub phylogenies 1 for each page
# Licuala c("Licuala lauterbachii var. bougainvillensis","Licuala petiolulata")
# Rest of Livistoninae c("Johannesteijsmannia perakensis","Acoelorraphe wrightii")
# Rhapidinae c("Rhapis humilis","Washingtonia filifera")
# CSP c("Coccothrinax salvatoris","Phoenix paludosa 1")
# Syncarpous clades c("Arenga ryukyuensis","Nannorrhops ritchieana")

# For each of these nodes I want to create a phylogeny which fits well on a page using the plotwe function.

# #########################################################################################
# ###################################--- Licuala ---#######################################
# #########################################################################################

# # Find MRCA node
# MRCA_licuala <- getMRCA(astral_tree_for_figure_orthologs, c("Licuala lauterbachii var. bougainvillensis","Licuala petiolulata")) 

# # processing phylogeny
# licuala_results <- process_phylogeny(MRCA_licuala, astral_tree_orthologs, astral_tree_orthologs_EN, tipcols_orthologs, data_orthologs,"Peru1",6)

# # Saving the results from process_phylogeny
# licuala_tree_figure <- licuala_results[[1]]
# licuala_edgecols_orthologs <- licuala_results[[2]]
# licuala_edgelabs_orthologs <- licuala_results[[3]]
# licuala_tipcols_orthologs <- licuala_results[[4]]

# # Renaming tips in tree
# licuala_tree_figure$tip.label = figurename_idx[licuala_tree_figure$tip.label]
# licuala_tree_figure$edge.length[which(licuala_tree_figure$edge.length == "NaN")] <- 0.00001

# #Plotting the sub phylogeny 
# pdf(paste(figurepath, "Licuala_tree_ortholog_genetrees_support_leftwards.pdf", sep = ""), paper = "a4", height = 11.7, width = 8.3)
# plotwe(licuala_tree_figure, direction = "rightwards", cex = 0.45, align.tip.label = T,
#        edge.color = licuala_edgecols_orthologs, link.color = licuala_tipcols_orthologs, label.offset = 0.08,
#        x.lim = c(0,4), edge.width = 2.5)
# edgelabels(text = licuala_edgelabs_orthologs, frame="none", cex=0.3, adj=c(0.5,-0.30)) # cex was 0.45
# #legend("bottomleft", legend=c("0-34","35-79","80-129","139-177","178-214","215-231"),fill=pal_orthologs,title="No. gene trees")
# #add.scale.bar(x=2, y=-1)
# dev.off()

# #########################################################################################
# #############################--- Rest of Livistoninae ---################################
# #########################################################################################

# # Find MRCA node
# MRCA_livistoninae <- getMRCA(astral_tree_for_figure_orthologs, c("Johannesteijsmannia perakensis","Acoelorraphe wrightii")) 

# # processing phylogeny
# livistoninae_results <- process_phylogeny(MRCA_livistoninae, astral_tree_orthologs, astral_tree_orthologs_EN, tipcols_orthologs, data_orthologs, "Peru1",6)

# # Saving the results from process_phylogeny
# livistoninae_tree_figure <- livistoninae_results[[1]]
# livistoninae_edgecols_orthologs <- livistoninae_results[[2]]
# livistoninae_edgecols_orthologs
# livistoninae_edgelabs_orthologs <- livistoninae_results[[3]]
# livistoninae_edgelabs_orthologs
# livistoninae_tipcols_orthologs <- livistoninae_results[[4]]
# livistoninae_tipcols_orthologs

# # Renaming tips in tree
# livistoninae_tree_figure$tip.label = figurename_idx[livistoninae_tree_figure$tip.label]
# livistoninae_tree_figure$edge.length[which(livistoninae_tree_figure$edge.length == "NaN")] <- 0.00001

# # Converting the tips from previous phylogeny to a triangle

# #Plotting the sub phylogeny 
# pdf(paste(figurepath, "livistoninae_tree_ortholog_genetrees_support_leftwards.pdf", sep = ""), paper = "a4", height = 11.7, width = 8.3)
# plotwe(livistoninae_tree_figure, direction = "rightwards", cex = 0.45, align.tip.label = T,
#        edge.color = livistoninae_edgecols_orthologs, link.color = livistoninae_tipcols_orthologs, label.offset = 0.08,
#        x.lim = c(0,4), edge.width = 2.5)
# edgelabels(text = livistoninae_edgelabs_orthologs, frame="none", cex=0.3, adj=c(0.5,-0.30)) # cex was 0.45
# #legend("bottomleft", legend=c("0-34","35-79","80-129","139-177","178-214","215-231"),fill=pal_orthologs,title="No. gene trees")
# #add.scale.bar(x=2, y=-1)
# dev.off()
# #########################################################################################
# #########################Co-phyloplot of the 2 previous trees############################
# #########################################################################################

cphyl <- cophylo(astral_tree_for_figure, astral_tree_for_figure_orthologs, fsize = 10, cex = 9)

pdf(paste(figurepath, "Coryphoideae_cophylo_wide.pdf", sep = ""), height = (70), width = (11*2))

plot(cphyl)

dev.off()

cphyl <- cophylo(astral_tree_for_figure, astral_tree_for_figure_orthologs, cex = 0.01)

cphyl <- cophylo(astral_tree_for_figure, astral_tree_for_figure_orthologs, fsize = 0.00000001)

pdf(paste(figurepath, "Coryphoideae_cophylo_small.pdf", sep = ""), height = (11.7*2), width = (8.3))

plot(cphyl, fsize = 0.3)

dev.off()

tips_per_page <- 97
total_tips <- length(astral_tree_for_figure$tip.label) # Assuming both trees have the same number of tips
pages <- ceiling(total_tips / tips_per_page)

for (i in 1:pages) {
  # Determine the range of tips for the current page
  start_tip <- (i - 1) * tips_per_page + 1
  end_tip <- min(i * tips_per_page, total_tips)
  
  # Subset the tips for both trees
  tips_subset <- astral_tree_for_figure$tip.label[start_tip:end_tip]

  tree1 <- keep.tip(astral_tree_for_figure, tips_subset)

  if (all(tips_subset %in% astral_tree_for_figure_orthologs$tip.label)) {
    tree2 <- keep.tip(astral_tree_for_figure_orthologs, tips_subset)
  } else {
    # If not all tips are in the second tree, subset the tips that are
    tips_subset <- tips_subset[tips_subset %in% astral_tree_for_figure_orthologs$tip.label]
    tree1 <- keep.tip(astral_tree_for_figure, tips_subset) # Update tree1 to match the new subset
    tree2 <- keep.tip(astral_tree_for_figure_orthologs, tips_subset)
  }
  
  # Create a cophylo object for the subset of tips
  cophy_subset <- cophylo(tree1, tree2)
  
  # Save each plot to a separate file
  pdf(paste0("cophylo_page_", i, ".pdf"), height = (11.7), width = (8.3))
  plot(cophy_subset, fsize = 0.5, main = paste("Page", i))
  dev.off()
}


# ##########################################################################################
# ##################################### Circle Tree ########################################
# ##########################################################################################
# pdf(paste(figurepath, "Coryphoideae_tree_ortholog_genetrees_support_leftwards_big_branches.pdf", sep = ""), height = (15.3*3), width = 11)
# plotwe(astral_tree_for_figure_orthologs, direction = "rightwards", cex = 0.35, align.tip.label = T, edge.color = edgecols_orthologs, link.color = tipcols_orthologs, label.offset = 0.08, x.lim = c(-1.2,max(distRoot(astral_tree_for_figure_orthologs))+1), edge.width = 3)
# #nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
# edgelabels(text = edgelabs_orthologs, frame="none", cex=0.45, adj=c(0.5,-0.30))
# legend("topleft", legend=c("0-25","26-64","65-111","112-168","169-206","207-230"),fill=pal_orthologs,title="No. gene trees")
# add.scale.bar(x=0, y=20)
# dev.off()

# nodeclass_orthologs$brks

# figurename_idx
############################################################################################
###################################-- Subfamilies --########################################
############################################################################################
#Sabaleae
  #Sabal 
  # Sabal Yapa Sabal Etonio
MRCA_sabal <- MRCA(astral_tree_for_figure_orthologs, c("Sabal yapa", "Sabal etonia")) #502

# Cryosophileae
  # Schippia, Itaya, Cryosophila, Trithrinax, Chelyocarpus, Zombia, Coccothrinax, Thrinax, Hemithrinax, Leucothrinax
  # Cryosophila nana, Coccothrinax macroglossa
MRCA_cryosophileae <- MRCA(astral_tree_for_figure_orthologs, c("Sabinaria magnifica", "Coccothrinax macroglossa")) #495

#Phoeniceae
  # Phoenix
  # Phoenix rupicola, Phoenix dactylifera
MRCA_phoeniceae <- MRCA(astral_tree_for_figure_orthologs, c("Phoenix rupicola", "Phoenix dactylifera")) #447

# Brahea
  # Brahea
  # Brahea aculeata , Brahea calcarea
MRCA_brahea <- MRCA(astral_tree_for_figure_orthologs,c("Brahea aculeata","Brahea calcarea"))

#Trachycarpeae
# Subtribe Rhapidinae
  # Chamaerops, Guihaia, Trachycarpus, Rhapidophyllum, Maxburretia, Rhapis, (Brahea, Colpothrinax)
  # Colpothrinax aphanopetala, Rhapis puhuongensis
MRCA_rhapidinae <- MRCA(astral_tree_for_figure_orthologs, c("Rhapidophyllum hystrix", "Rhapis humilis")) #636

#Subtribe livistoninae 
  #Livistona, Licuala, johannesteijsmannia, Pholidocarpus, Saribus, Acoelorrhaphe, Serenoa
  # Serenoa repens, Licuala simplex
MRCA_livistoninae <-MRCA(astral_tree_for_figure_orthologs, c("Livistona carinensis", "Licuala simplex")) #636

# Chuniophoeniceae
  # Chuniophoenix, Kerriodoxa, Nannorrhops, Tahina
  # Nannorrhops ritchieana, chuniophoenix hainanensis
MRCA_chuniophoeniceae <- MRCA(astral_tree_for_figure_orthologs, c("Nannorrhops ritchieana", "Chuniophoenix hainanensis")) #409

# Caryoteae
  # Caryota, Wallichia, Arenga
  # Caryota obtusa, Wallichia gracilis
MRCA_caryoteae <- MRCA(astral_tree_for_figure_orthologs, c("Caryota obtusa", "Wallichia gracilis")) #436

#Corypheae
  # Corypha
  # Corypha lecomtei, Corypha taliera
MRCA_corypheae <- MRCA(astral_tree_for_figure_orthologs, c("Corypha lecomtei 1", "Corypha taliera")) #415

# Borasseae
  # Subtribe Hyphaeninae
  # Bismarckia, Satranala, Hyphaene, Medemia
  # Bismarkia nobilis, Hyphaene thebaica
MRCA_borasseae <- MRCA(astral_tree_for_figure_orthologs, c("Bismarckia nobilis", "Hyphaene thebaica")) #419

  # Subtribe Lataniinae
  # Borassodendron, Latania, Borassus, Lodoicea
  # Borassus aethiopum, Latania lontaroides
MRCA_lataniinae <- MRCA(astral_tree_for_figure_orthologs, c("Borassus aethiopum", "Latania lontaroides")) #427



dat_subfam <- data.frame(
           node = c(MRCA_sabal,MRCA_cryosophileae,MRCA_phoeniceae,MRCA_rhapidinae,MRCA_livistoninae,MRCA_chuniophoeniceae,MRCA_caryoteae,MRCA_corypheae,MRCA_borasseae,MRCA_lataniinae),
           name = c("Sabaleae","Cryosophileae","Phoeniceae","Rhapidineae","Livistoninae","Chuniophoeniceae","Caryoteae","Corypheae", "Hyphaeninae","Lataniinae"))








##########################################################################################
################################### split up plot ########################################
###################################-- Orthologs --########################################
##########################################################################################

test_tree <- ladderize(astral_tree_for_figure_orthologs)

# Creating a small vector of colours and states
pal_orthologs_t <- pal_orthologs
pal_orthologs_t
pal_orthologs_t <- setNames(pal_orthologs,1:6)

# and the inverse vector with states and colours
pal_orthologs_t1 <- (1:6)
names(pal_orthologs_t1) <- pal_orthologs

# Painting the branches based on the number of genes in the orthologs gene tree
for (i in (2:length(astral_tree_orthologs_EN$node.label) + Ntip(astral_tree_orthologs_EN))){
	if(i %in% names(ENcolour_orthologs)){
		print(i)
		col <- ENcolour_orthologs[which(names(ENcolour_orthologs) == i)] # finding the colour for the node in ENcolour_orthologs
		nr <- which(toupper(names(pal_orthologs_t1)) == col[[1]][1]) # converting that colour to a number based on pal_orthologs_t1
		test_tree <- paintBranches(test_tree, edge = i, state = nr,) # painting the state on the branches of the tree
	}
}

# Painting the tips based on the number of genes in the orthologs gene tree
for (i in 1:length(test_tree$tip.label)){
	tip_name <- names(test_tree$tip.label[i])
	if(tip_name %in% names(tipcols_orthologs)){
		col_tip <- tipcols_orthologs[which(names(tipcols_orthologs) == tip_name)] # finding the colour for the tip in tipscols_orthologs
		nr_tip <- which(toupper(names(pal_orthologs_t1)) == col_tip[[1]][1]) # converting that colour to a number based on pal_orthologs_t1
		cat(tip_name," ",col_tip," ",nr_tip,"\n")
		test_tree <- paintBranches(test_tree, edge = i, state = nr_tip,) # painting the state on the branches of the tree
	}
}

#Splitting the tree into just 1 part
plotTree_clades(test_tree, ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 2, cex = 0.70, file = "split_plot_cory_orthologs_1_part.pdf", width = 15, height = 60, tipcols = tipcols_orthologs, type = "cladogram")

# Splitting the tree into 5 parts
split.plotTree_clades(test_tree_no_branch_lenghts, splits = c(0.2819,0.398,0.6081,0.845), ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory_orthologs.pdf", x_lim = c(0,23), tipcols = tipcols_orthologs)

# Splitting the tree into 5 parts
split.plotTree_clades(test_tree, splits = c(0.2819,0.398,0.6081,0.845), ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory_orthologs.pdf", x_lim = c(0,23), tipcols = tipcols_orthologs)

 # Splitting the tree into 8 equal parts
split.plotTree_clades(test_tree, splits = c(0.125,0.25,0.375,0.5,0.625,0.75,0.875), ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory_orthologs_8_parts.pdf", x_lim = c(0,23), tipcols = tipcols_orthologs)

# Adding clade labels
split.plotTree_clades(test_tree, splits = c(0.2819,0.398,0.6081,0.845), ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory_orthologs_clades.pdf", x_lim = c(0,23), tipcols = tipcols_orthologs, clades = dat_subfam)

c(0.2819,0.398,0.6081,0.845) # Ladderized
c(0.281,0.397,0.6065,0.8427) # Not Ladderized

# Editing the edge lengths so that the ones that are 0.001000000 are now 0.1
test_tree$edge.length
#test_tree$edge.length[test_tree$edge.length == 0.100000000] <- 0.001
test_tree$edge.length[test_tree$edge.length == 0.000000000] <- 0.001
test_tree$edge.length[order(test_tree$edge.length)]
# Assigning the legend to a variable
# Plotting the tree

plot(test_tree, colors	= pal_orthologs_t)
View(test_tree)
legend_plot <- legend("bottomleft", legend=c("0-25","26-64","65-111","112-168","169-206","207-230"), fill=pal_orthologs, title="No. gene trees", horiz=TRUE)
legend_plot

# When you need to convert the pages from the pdf to indidual jpeg files use the following command in the terminal
# mkdir -p images && pdftoppm -jpeg -r 300 split_plot_cory_orthologs.pdf images/pg


##########################################################################################
################################### split up plot ########################################
###################################-- All Genes --########################################
##########################################################################################

test_tree_ag <- ladderize(astral_tree_for_figure)

"Maxburretia gracilis" %in% test_tree_ag$tip.label
"Arenga australasica" %in% test_tree_ag$tip.label


# Creating a small vector of colours and states
pal_t <- pal
pal_t <- setNames(pal,1:6)

# and the inverse vector with states and colours
pal_t1 <- (1:6)
names(pal_t1) <- pal

# Painting the branches based on the number of genes in the orthologs gene tree
for (i in (2:length(astral_tree_EN$node.label) + Ntip(astral_tree_EN))){
	if(i %in% names(ENcolour)){
		#print(i)
		col <- ENcolour[which(names(ENcolour) == i)] # finding the colour for the node in ENcolour_orthologs
		nr <- which(toupper(names(pal_t1)) == col[[1]][1]) # converting that colour to a number based on pal_orthologs_t1
		test_tree_ag <- paintBranches(test_tree_ag, edge = i, state = nr,) # painting the state on the branches of the tree
	}
}

# Painting the tips based on the number of genes in the orthologs gene tree
for (i in 1:length(test_tree_ag$tip.label)){
	tip_name <- names(test_tree_ag$tip.label[i])
	if(tip_name %in% names(tipcols)){
		col_tip <- tipcols[which(names(tipcols) == tip_name)] # finding the colour for the tip in tipscols_orthologs
		nr_tip <- which(toupper(names(pal_t1)) == col_tip[[1]][1]) # converting that colour to a number based on pal_orthologs_t1
		cat(tip_name," ",col_tip," ",nr_tip,"\n")
		test_tree_ag <- paintBranches(test_tree_ag, edge = i, state = nr_tip,) # painting the state on the branches of the tree
	}
}


# Splitting the tree into 5 parts
split.plotTree_clades(test_tree_ag, splits = c(0.2818,0.399,0.6092,0.8463), ftype = "i", mar = c(3,1,1,1), color = pal_t,
 fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory.pdf", x_lim = c(0,16.5), tipcols = tipcols)

# # Adding clade labels
# split.plotTree_clades(test_tree_ag, splits = c(0.2819,0.398,0.6081,0.845), ftype = "i", mar = c(3,1,1,1), color = pal_t,
#  fsize = 1, type = "phylogram", cex = 0.5, file = "split_plot_cory_clades.pdf", x_lim = c(0,23), tipcols = tipcols, clades = dat_subfam)

# c(0.2819,0.398,0.6081,0.845) # Ladderized
# c(0.281,0.397,0.6065,0.8427) # Not Ladderized

# # Editing the edge lengths so that the ones that are 0.001000000 are now 0.1
# test_tree_ag$edge.length
# #test_tree$edge.length[test_tree$edge.length == 0.100000000] <- 0.001
# test_tree_ag$edge.length[test_tree_ag$edge.length == 0.000000000] <- 0.001
# test_tree_ag$edge.length[order(test_tree_ag$edge.length)]
# # Assigning the legend to a variable
# # Plotting the tree
# plot(test_tree_ag, colors	= pal_t)

# legend_plot <- legend("bottomleft", legend=c("0-25","26-64","65-111","112-168","169-206","207-230"), fill=pal_orthologs, title="No. gene trees", horiz=TRUE)
# legend_plot


# length(test_tree_ag$tip.label)



# # Creating a small phylogeny which I can add to the pictures
# plot(astral_tree_orthologs, show.tip.label = FALSE, use.edge.length = FALSE, edge.width = 10)

# library(ggtree)
# p <- ggtree(astral_tree_orthologs, ladderize = FALSE, branch.length = "none", size = 3)

# p + theme(legend.position = "none")

# # plotting this picture as a png with no background and the same dimensions as the tree used in the article.
# png(paste(figurepath, "phylogeny_orthologs.png", sep = ""), width = 1000, height = 4000, units = "px", bg = "transparent")
# p + theme(legend.position = "none")
# dev.off()

# Make a new figure which shows the phylogeny of the genera



####################################################################################################################################
#Testing area
test_tree_no_branch_lenghts <- test_tree
length(test_tree_no_branch_lenghts$edge.length)
test_tree_no_branch_lenghts$edge.length <- rep(0.25,959)

test_tree_no_branch_lenghts$edge.length[which(test_tree_no_branch_lenghts$edge.length != 0)] <- 0.25

depth <- max(node.depth.edgelength(test_tree))
new_edge_lengths <- rep(depth / (Ntip(test_tree) + test_tree$Nnode), length(test_tree$edge.length))

tree_ultrametric_test <- force.ultrametric(test_tree)


# Splitting the tree into 5 parts
split.plotTree_clades(test_tree_no_branch_lenghts, splits = c(0.2819,0.398,0.6081,0.845), ftype = "i", mar = c(3,1,1,1), color = pal_orthologs_t,
 fsize = 1, cex = 0.5, file = "split_plot_cory_orthologs.pdf", x_lim = c(0,23), tipcols = tipcols_orthologs)

par(fg = "transparent")
plotTree(test_tree_no_branch_lenghts, align.tip.label = TRUE)

plotSimmap(test_tree_no_branch_lenghts)

plotSimmap(test_tree)

