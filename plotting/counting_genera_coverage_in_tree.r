# This script is used to count the number of species, subspecies and variants in the final tree and the ortholog tree.
# This is done using the latest version of the World Checklist of Vascular Plants (WCVP).

# Loading some command line elements
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("ape","phytools","data.table","tidyr")

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

astral_tree_orthologs <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree_orthologs.tre")
astral_tree_orthologs <- drop.tip(astral_tree_orthologs, "4039")

astral_tree_orthologs <- root(astral_tree_orthologs, outgroup=c("1081"))
astral_tree_orthologs <- ladderize(astral_tree_orthologs, right=TRUE)
astral_tree_orthologs$tip.label = figurename_idx[astral_tree_orthologs$tip.label]


test_tree <- astral_tree_orthologs


# List of Genera in Coryphoideae
genera <- c("Sabal","Schippia","Trithrinax","Zombia","Coccothrinax","Hemithrinax","Leucothrinax","Thrinax","Chelyocarpus",
            "Cryosophila","Itaya","Phoenix","Chamaerops","Guihaia","Trachycarpus","Rhapidophyllum","Maxburretia","Rhapis",
            "Livistona","Licuala","Johannesteijsmannia","Pholidocarpus","Acoelorraphe","Serenoa","Brahea",
            "Colpothrinax","Copernicia","Pritchardia","Washingtonia","Chuniophoenix","Kerriodoxa","Nannorrhops","Tahina",
            "Caryota","Arenga","Wallichia","Corypha","Bismarckia","Satranala","Hyphaene","Medemia","Latania","Lodoicea",
            "Borassodendron","Borassus","Lanonia","Saribus","Sabinaria")

monotypic_genera <- c("Schippia", "Zombia", "Leucothrinax","Itaya","Chamaerops","Rhapidophyllum", "Acoelorraphe","Serenoa","Washingtonia","Kerriodoxa","Nannorrhops","Tahina",
                      "Bismarckia","Satranala","Medemia","Lodoicea","Sabinaria")

# Remove monotypic genera from genera.
genera <- genera[!genera %in% monotypic_genera]

# Create subset of figurename_idx which only contains the tips which are in the tree
figurename_idx_sub <- figurename_idx[figurename_idx %in% test_tree$tip.label]

# Now I also need to remove samples from duplicate species.
# I will start by finding all the tip names which includes either a number, var. or subsp.
duplicate_species <- figurename_idx_sub[grepl("[0-9]|var.|subsp.", figurename_idx_sub)]

# Now I group them based on the first 2 elements of their name and create a list of lists where each list within the list contains the names of duplicate samples of a species.
# Splitting them based on genus names
duplicate_species_split <- split(duplicate_species, substr(duplicate_species, 1, 5))
duplicate_species_split

list_of_samples <- strsplit(duplicate_species, names(duplicate_species),split = " ")
list_of_samples

unique_names_duplicates <- c()
# I think ill do it with a for loop
for ( i in 1:length(list_of_samples)){
  print(paste("Genus:", list_of_samples[[i]][[1]]))
  print(paste("Species:", paste(list_of_samples[[i]][[2]], collapse=", ")))

  # Combine them
  sp_name <- paste(list_of_samples[[i]][[1]], list_of_samples[[i]][[2]], sep=" ")

  # find out if the sp_name is already in unique_names_duplicates.
  if (sp_name %in% unique_names_duplicates){
    #print(paste("Error: sp_name", sp_name, "is already in unique_names_duplicates"))
  } else {
    unique_names_duplicates <- c(unique_names_duplicates, sp_name)
  }
}
unique_names_duplicates

# Now for each of these unique names duplicates I need to add the species to a list.
# I will do this by looping through the unique_names_duplicates and finding the samples in list_of_samples which contains the unique name.
# the tips which are found will be added to a list. which will then be added to a list of lists.

list_of_species <- list()
for (i in 1:length(unique_names_duplicates)){
#   print(paste("Genus:", unique_names_duplicates[i]))
  species <- duplicate_species[grepl(unique_names_duplicates[[i]], duplicate_species)]
  list_of_species[[i]] <- species
  names(list_of_species)[i] <- unique_names_duplicates[i]
}

# I will remove species where the name contains sp.
not_identified_species <- list_of_species[grepl("sp.", names(list_of_species))] # Remove from the tree and not from the list

# Remove the sp. species from the tree
test_tree <- drop.tip(test_tree, not_identified_species[[1]])

# I need to do some manual edits of this list
#I just want to check if there are any tips in the tree which have the same name as the name of each of the lists
# I will do this by looping through the list_of_species and checking if the names of the list is in the tree.
# if they are I will add them to the list of tips with that name.
for ( i in 1:length(list_of_species)){
  test_name <- names(list_of_species)[i]
  #print(test_name)
  # Find the tips which are in the tree
  if ( test_name %in% tree$tip.label){
    print(paste("Adding tip ", test_name, " to the list"))
    list_of_species[[i]] <- c(list_of_species[[i]], test_name)
  }
} 
list_of_species

# will remove the lists which contain only 1 element
list_of_species <- list_of_species[lengths(list_of_species) > 1]

# I will remove species where the name contains sp. from the list aswell as we have removed all of these.
list_of_species <- list_of_species[!grepl("sp.", names(list_of_species))] # Remove from the tree and not from the list

# Now I need to randomly remove the duplicate samples from the tree.
# I will do this 100 times in order to create 100 different trees.
# I will then calculate the average local posterior probability for each genus in each of the 100 trees.

list_of_trees <-  as.list(NA, 100)
removed_tips <- as.list(NA, 100)
kept_tips <- as.list(NA, 100)

for (i in 1:100){ # 100 replicates
  loop_tree <- test_tree

  internal_removed_tips <- c()
  internal_kept_tips <- c()

  for (j in 1:length(list_of_species)){ # each tip which there are duplicates of.
    # Select 1 tip at random from list j of list_of_species
    tip_keep <- sample(list_of_species[[j]], 1)

    # Select the other tips from the list
    tips_remove <- list_of_species[[j]][!list_of_species[[j]] %in% tip_keep]

    # Remove the tips from the tree
    loop_tree <- drop.tip(loop_tree, tips_remove)

    # Add the removed tips to the list
    internal_removed_tips <- c(internal_removed_tips, tips_remove)

    # Add the kept tips to the list
    internal_kept_tips <- c(internal_kept_tips, tip_keep)

  }
  # Adding the loop tree to the list.
  list_of_trees[[i]] <- loop_tree

  # Adding the removed tips to the list
  removed_tips[[i]] <- internal_removed_tips

  # Adding the kept tips to the list
  kept_tips[[i]] <- internal_kept_tips
}

# Testing if all the elements of kept tips are the same length
all(lengths(kept_tips) == lengths(kept_tips)[1]) 

# Testing if all the elements of removed tips are the same length
all(lengths(removed_tips) == lengths(removed_tips)[1])

# Testing if there are the same number of tips in all the trees.
all(lengths(list_of_trees$tip.label) == lengths(list_of_trees$tip.label)[1])



# For each of these genera I want to find the MRCA of all tips in that genus
# I want to check if all the descendants of that MRCA are in the same genus
# and then I want to calculate the average local posterior probability for that genus
list_of_results <- as.list(NA, 100)
for (q in 1:length(list_of_trees)){
  results <- data.frame()

  for (genus in genera){
    print(paste("Genus:", genus))

    tips <- unname(figurename_idx[which(grepl(genus, figurename_idx))])
    tips <- as.vector(tips)
    tips

    # are all the tips found in the tree?
    if(!all(tips %in% list_of_trees[[q]]$tip.label)){
      print(paste("Error: not all tips found in the tree for genus", genus))
    
      # Which tips are not found in the tree?
      print(paste("Tips not found in the tree for genus", genus, ":", paste(tips[!tips %in% list_of_trees[[q]]$tip.label], collapse=", ")))

      # Removing these tips from the tips vector
      tips <- tips[tips %in% list_of_trees[[q]]$tip.label]
    } else {
      print(paste("All tips found in the tree for genus", genus))
    }


    # I need to remove genera with one 1 tip as I cannot find their MRCA
    if (length(tips) == 1){
      print(paste("Error: only 1 tip found in the tree for genus", genus))
      next
    }
    
    # Find the MRCA of all tips in the genus
    mrca <- getMRCA(list_of_trees[[q]], tips)
    print(paste("MRCA for genus", genus, ":", mrca))
    
    if (is.na(mrca)) {
      print(paste("Error: MRCA is NA for genus", genus))
      next
    } else {
      "MRCA is not NA"
    }
    
    # Check if all descendants of the MRCA are in the same genus
    descendants <- getDescendants(list_of_trees[[q]], node = mrca)
    #print(paste("Descendants for genus", genus, ":", descendants))
    tip_labels <- unname(list_of_trees[[q]]$tip.label[descendants[descendants <= length(list_of_trees[[q]]$tip.label)]])

    
    if (any(is.na(descendants))) {
      print(paste("Error: NA in descendants for genus", genus))
      next
    } else {
      "No descendants are NA"
    }
    
    # Find all descendant tips
    # Split the labels by space and extract the genus names
    descendants_genera <- unique(sapply(tip_labels, function(x) strsplit(x, " ")[[1]][1]))
    if (length(unique(descendants_genera)) > 1){
      print(paste("Error: Descendants of MRCA of", genus, "are not all in the same genus"))
    } else {
      print(paste("Descendants of MRCA of", genus, "are all in the same genus"))
    }
    
    # Create a subtree which contains only the descendants of the MRCA
    subtree <- extract.clade(list_of_trees[[q]], node = mrca)

    # Calculate the average local posterior probability for the genus
    local_posterior_probabilities <- subtree$node.label
    print(paste("Local posterior probabilities for genus", genus, ":", paste(local_posterior_probabilities, collapse=", ")))
    
    if (any(is.na(local_posterior_probabilities))) {
      print(paste("Error: NA in local posterior probabilities for genus", genus))
      next
    }
    
    average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
    
    # Print the results
    print(paste("Genus:", genus))
    print(paste("Average local posterior probability:", average_local_posterior_probability))
    print(paste("Number of species:", length(tip_labels)))
    #print(paste("Species:", paste(figurename_idx[descendants], collapse=", ")))
    
    # Save the results
    output <- data.frame(genus = genus, average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))

    results <- rbind(results, output)
  }

  # Can I order the list aphabetically?
  results <- results[order(results$genus),]
  results

  # Calculate support values for Arenga and Livistona without the following species 
  # Livistona exigua
  # Arenga distincta & Arenga hastata

  # Arenga
  tips <- unname(figurename_idx[which(grepl("Arenga", figurename_idx))])
  tips <- as.vector(tips)
  tips

  # Are all the tips found in the tree?
  if(!all(tips %in% list_of_trees[[q]]$tip.label)){
    print("Error: not all tips found in the tree for genus Arenga")
    
    # Which tips are not found in the tree?
    print(paste("Tips not found in the tree for genus Arenga:", paste(tips[!tips %in% list_of_trees[[q]]$tip.label], collapse=", ")))

    # Removing these tips from the tips vector
    tips <- tips[tips %in% list_of_trees[[q]]$tip.label]
  } else {
    print("All tips found in the tree for genus Arenga")
  }

  # I need to remove genera with one 1 tip as I cannot find their MRCA
  if (length(tips) == 1){
    print("Error: only 1 tip found in the tree for genus Arenga")
    next
  }

  # Remove Arenga distincta and Arenga hastata
  tips <- tips[!grepl("distincta", tips)]
  tips <- tips[!grepl("hastata", tips)]
  tips

  # Find the MRCA of all tips in the genus
  mrca <- getMRCA(list_of_trees[[q]], tips)
  print(paste("MRCA for genus Arenga:", mrca))

  if (is.na(mrca)) {
    print("Error: MRCA is NA for genus Arenga")
  } else {
    "MRCA is not NA"
  }

  # Check if all descendants of the MRCA are in the same genus
  descendants <- getDescendants(list_of_trees[[q]], node = mrca)

  tip_labels <- unname(list_of_trees[[q]]$tip.label[descendants[descendants <= length(list_of_trees[[q]]$tip.label)]])

  # Create subtree
  subtree <- extract.clade(list_of_trees[[q]], node = mrca)

  # Calculate the average local posterior probability for the genus
  local_posterior_probabilities <- subtree$node.label

  # find mean
  average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
  print(paste("Average local posterior probability for Arenga:", average_local_posterior_probability))

  # and the number of species
  print(paste("Number of species for Arenga:", length(descendants)))
  length(tip_labels)

  # add the new results for Arenga to the results dataframe
  output <- data.frame(genus = "Arenga", average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))
  # Update the results dataframe
  results[results$genus == "Arenga",] <- output
  results

  # Livistona
  # Livistona exigua
  tips <- unname(figurename_idx[which(grepl("Livistona", figurename_idx))])
  tips <- as.vector(tips)
  tips

  # Are all the tips found in the tree?
  if(!all(tips %in% list_of_trees[[q]]$tip.label)){
    print("Error: not all tips found in the tree for genus Arenga")
    
    # Which tips are not found in the tree?
    print(paste("Tips not found in the tree for genus Arenga:", paste(tips[!tips %in% list_of_trees[[q]]$tip.label], collapse=", ")))

    # Removing these tips from the tips vector
    tips <- tips[tips %in% list_of_trees[[q]]$tip.label]
  } else {
    print("All tips found in the tree for genus Arenga")
  }

  # I need to remove genera with one 1 tip as I cannot find their MRCA
  if (length(tips) == 1){
    print("Error: only 1 tip found in the tree for genus Arenga")
    next
  }

  # Remove Arenga distincta and Arenga hastata
  tips <- tips[!grepl("exigua", tips)]
  tips

  # Find the MRCA of all tips in the genus
  mrca <- getMRCA(list_of_trees[[q]], tips)
  print(paste("MRCA for genus Livistona:", mrca))

  if (is.na(mrca)) {
    print("Error: MRCA is NA for genus Livistona")
  } else {
    "MRCA is not NA"
  }

  # Check if all descendants of the MRCA are in the same genus
  descendants <- getDescendants(list_of_trees[[q]], node = mrca)

  tip_labels <- unname(list_of_trees[[q]]$tip.label[descendants[descendants <= length(list_of_trees[[q]]$tip.label)]])

  # Create subtree
  subtree <- extract.clade(list_of_trees[[q]], node = mrca)

  # Calculate the average local posterior probability for the genus
  local_posterior_probabilities <- subtree$node.label

  # find mean
  average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
  print(paste("Average local posterior probability for Livistona:", average_local_posterior_probability))

  # and the number of species
  print(paste("Number of species for Livistona:", length(tip_labels)))
  length(tip_labels)
  tip_labels

  # add the new results for Arenga to the results dataframe
  output <- data.frame(genus = "Livistona", average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))
  # Update the results dataframe
  results[results$genus == "Livistona",] <- output
  results



  # Can I make a for loop for the genera which finds the support value for the missing genera

  # Finding the genera which are NOT in results$genus
  miss_genera <- genera[!genera %in% results$genus]
  miss_genera


  for (genus in miss_genera){
    print(paste("Genus:", genus))

    tips <- unname(figurename_idx[which(grepl(genus, figurename_idx))])
    tips <- as.vector(tips)
    tips

    # are all the tips found in the tree?
    if(!all(tips %in% list_of_trees[[q]]$tip.label)){
      print(paste("Error: not all tips found in the tree for genus", genus))
    
      # Which tips are not found in the tree?
      print(paste("Tips not found in the tree for genus", genus, ":", paste(tips[!tips %in% list_of_trees[[q]]$tip.label], collapse=", ")))

      # Removing these tips from the tips vector
      tips <- tips[tips %in% list_of_trees[[q]]$tip.label]
    } else {
      print(paste("All tips found in the tree for genus", genus))
    }
    
    # Find the tip index
    tip_index <- unname(which(list_of_trees[[q]]$tip == tips))
    tip_index

    # Find the parent node
    parent_node <- list_of_trees[[q]]$edge[list_of_trees[[q]]$edge[,2] == tip_index, 1]
    parent_node <- parent_node - length(list_of_trees[[q]]$tip.label)
    parent_node

    lpp_val <- list_of_trees[[q]]$node.label[parent_node]

    
    # Print the results
    print(paste("Genus:", genus))
    print(paste("Average local posterior probability:", lpp_val))
    print(paste("Number of species:", length(tips)))
    #print(paste("Species:", paste(figurename_idx[descendants], collapse=", ")))
    
    # Save the results
    output <- data.frame(genus = genus, average_local_posterior_probability = lpp_val, number_of_species = length(tips))

    results <- rbind(results, output)
  }

  # # Adding the missing species
  # miss_genera <- data.frame(genus = miss_genera, average_local_posterior_probability = NA, number_of_species = 1)
  # results <- rbind(results, miss_genera)

  results <- results[order(results$genus),]
  list_of_results[[q]] <- results
}

# Now I need to calculate the average support values for each genus in each of the 100 trees.
# I will do this by calculating the mean of the average local posterior probabilities for each genus in each of the 100 trees.
# I will also calculate the standard deviation of the average local posterior probabilities for each genus in each of the 100 trees.
average_local_across_100_trees <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(average_local_across_100_trees) <- c("genus", "mean_average_lpp", "sd_average_lpp", "number_of_species")

for (i in 1:length(genera)) {
  genus <- genera[i]
  print(paste("Genus:", genus))
  mean_average_lpp <- mean(sapply(list_of_results, function(x) x$average_local_posterior_probability[x$genus == genus]))
  sd_average_lpp <- sd(sapply(list_of_results, function(x) x$average_local_posterior_probability[x$genus == genus]))
  number_of_species <- mean(sapply(list_of_results, function(x) x$number_of_species[x$genus == genus]))
  average_local_across_100_trees <- rbind(average_local_across_100_trees, data.frame(genus = genus, mean_average_lpp = mean_average_lpp, sd_average_lpp = sd_average_lpp, number_of_species = number_of_species))
}

average_local_across_100_trees

# write csv file
write.csv(average_local_across_100_trees, file = "/home/au543206/Documents/Coryphoideae/Figures_in_r/data/average_lpp_per_genera_orthologs_100_trees.csv", row.names = FALSE)


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

# Doing the same thing but in the all genes tree.
astral_tree <- ape::read.tree("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/astral_tree.tre")
astral_tree <- drop.tip(astral_tree, "4039") # Maxburretia gracilis
astral_tree <- drop.tip(astral_tree, "4091") # Arenga australasica

astral_tree <- root(astral_tree, outgroup=c("1081"))
astral_tree <- ladderize(astral_tree, right=TRUE)
astral_tree$tip.label = figurename_idx[astral_tree$tip.label]

# List of Genera in Coryphoideae
genera <- c("Sabal","Schippia","Trithrinax","Zombia","Coccothrinax","Hemithrinax","Leucothrinax","Thrinax","Chelyocarpus",
            "Cryosophila","Itaya","Phoenix","Chamaerops","Guihaia","Trachycarpus","Rhapidophyllum","Maxburretia","Rhapis",
            "Livistona","Licuala","Johannesteijsmannia","Pholidocarpus","Acoelorraphe","Serenoa","Brahea",
            "Colpothrinax","Copernicia","Pritchardia","Washingtonia","Chuniophoenix","Kerriodoxa","Nannorrhops","Tahina",
            "Caryota","Arenga","Wallichia","Corypha","Bismarckia","Satranala","Hyphaene","Medemia","Latania","Lodoicea",
            "Borassodendron","Borassus","Lanonia","Saribus","Sabinaria")

monotypic_genera <- c("Schippia", "Zombia", "Leucothrinax","Itaya","Chamaerops","Rhapidophyllum", "Acoelorraphe","Serenoa","Washingtonia","Kerriodoxa","Nannorrhops","Tahina",
                      "Bismarckia","Satranala","Medemia","Lodoicea","Sabinaria")

# Remove monotypic genera from genera.
genera <- genera[!genera %in% monotypic_genera]     

# Create subset of figurename_idx which only contains the tips which are in the tree
figurename_idx_sub <- figurename_idx[figurename_idx %in% astral_tree$tip.label]

# Now I also need to remove samples from duplicate species.
# I will start by finding all the tip names which includes either a number, var. or subsp.
duplicate_species <- figurename_idx_sub[grepl("[0-9]|var.|subsp.", figurename_idx_sub)]

# Now I group them based on the first 2 elements of their name and create a list of lists where each list within the list contains the names of duplicate samples of a species.
# Splitting them based on genus names
duplicate_species_split <- split(duplicate_species, substr(duplicate_species, 1, 5))
duplicate_species_split

list_of_samples <- strsplit(duplicate_species, names(duplicate_species),split = " ")
list_of_samples

unique_names_duplicates <- c()
# I think ill do it with a for loop
for ( i in 1:length(list_of_samples)){
  #print(paste("Genus:", list_of_samples[[i]][[1]]))
  #print(paste("Species:", paste(list_of_samples[[i]][[2]], collapse=", ")))

  # Combine them
  sp_name <- paste(list_of_samples[[i]][[1]], list_of_samples[[i]][[2]], sep=" ")

  # find out if the sp_name is already in unique_names_duplicates.
  if (sp_name %in% unique_names_duplicates){
    #print(paste("Error: sp_name", sp_name, "is already in unique_names_duplicates"))
  } else {
    unique_names_duplicates <- c(unique_names_duplicates, sp_name)
  }
}
unique_names_duplicates

# Now for each of these unique names duplicates I need to add the species to a list.
# I will do this by looping through the unique_names_duplicates and finding the samples in list_of_samples which contains the unique name.
# the tips which are found will be added to a list. which will then be added to a list of lists.

list_of_species <- list()
for (i in 1:length(unique_names_duplicates)){
#   print(paste("Genus:", unique_names_duplicates[i]))
  species <- duplicate_species[grepl(unique_names_duplicates[[i]], duplicate_species)]
  list_of_species[[i]] <- species
  names(list_of_species)[i] <- unique_names_duplicates[i]
}

# I will remove species where the name contains sp.
not_identified_species <- list_of_species[grepl("sp.", names(list_of_species))] # Remove from the tree and not from the list

# Remove the sp. species from the tree
astral_tree <- drop.tip(astral_tree, not_identified_species[[1]])

# I need to do some manual edits of this list
#I just want to check if there are any tips in the tree which have the same name as the name of each of the lists
# I will do this by looping through the list_of_species and checking if the names of the list is in the tree.
# if they are I will add them to the list of tips with that name.
for ( i in 1:length(list_of_species)){
  test_name <- names(list_of_species)[i]
  #print(test_name)
  # Find the tips which are in the tree
  if ( test_name %in% tree$tip.label){
    print(paste("Adding tip ", test_name, " to the list"))
    list_of_species[[i]] <- c(list_of_species[[i]], test_name)
  }
} 
list_of_species

# will remove the lists which contain only 1 element
list_of_species <- list_of_species[lengths(list_of_species) > 1]

# I will remove species where the name contains sp. from the list aswell as we have removed all of these.
list_of_species <- list_of_species[!grepl("sp.", names(list_of_species))] # Remove from the tree and not from the list

# Now I need to randomly remove the duplicate samples from the tree.
# I will do this 100 times in order to create 100 different trees.
# I will then calculate the average local posterior probability for each genus in each of the 100 trees.

list_of_trees_all_genes <-  as.list(NA, 100)
removed_tips <- as.list(NA, 100)
kept_tips <- as.list(NA, 100)

for (i in 1:100){ # 100 replicates
  loop_tree <- astral_tree

  internal_removed_tips <- c()
  internal_kept_tips <- c()

  for (j in 1:length(list_of_species)){ # each tip which there are duplicates of.
    # Select 1 tip at random from list j of list_of_species
    tip_keep <- sample(list_of_species[[j]], 1)

    # Select the other tips from the list
    tips_remove <- list_of_species[[j]][!list_of_species[[j]] %in% tip_keep]

    # Remove the tips from the tree
    loop_tree <- drop.tip(loop_tree, tips_remove)

    # Add the removed tips to the list
    internal_removed_tips <- c(internal_removed_tips, tips_remove)

    # Add the kept tips to the list
    internal_kept_tips <- c(internal_kept_tips, tip_keep)

  }
  # Adding the loop tree to the list.
  list_of_trees_all_genes[[i]] <- loop_tree

  # Adding the removed tips to the list
  removed_tips[[i]] <- internal_removed_tips

  # Adding the kept tips to the list
  kept_tips[[i]] <- internal_kept_tips
}

# Testing if all the elements of kept tips are the same length
all(lengths(kept_tips) == lengths(kept_tips)[1]) 

# Testing if all the elements of removed tips are the same length
all(lengths(removed_tips) == lengths(removed_tips)[1])

# Testing if there are the same number of tips in all the trees.
all(lengths(list_of_trees_all_genes$tip.label) == lengths(list_of_trees$tip.label)[1])


# For each of these genera I want to find the MRCA of all tips in that genus
# I want to check if all the descendants of that MRCA are in the same genus
# and then I want to calculate the average local posterior probability for that genus
list_of_results_all_genes <- as.list(NA, 100)

for (q in 1:length(list_of_trees_all_genes)){
  results <- data.frame()
  for (genus in genera){
    print(paste("Genus:", genus))

    tips <- unname(figurename_idx[which(grepl(genus, figurename_idx))])
    tips <- as.vector(tips)
    tips

    # are all the tips found in the tree?
    if(!all(tips %in% list_of_trees_all_genes[[q]]$tip.label)){
      print(paste("Error: not all tips found in the tree for genus", genus))

      # Which tips are not found in the tree?
      print(paste("Tips not found in the tree for genus", genus, ":", paste(tips[!tips %in% list_of_trees_all_genes[[q]]$tip.label], collapse=", ")))

      # Removing these tips from the tips vector
      tips <- tips[tips %in% list_of_trees_all_genes[[q]]$tip.label]
    } else {
      print(paste("All tips found in the tree for genus", genus))
    }


    # I need to remove genera with one 1 tip as I cannot find their MRCA
    if (length(tips) == 1){
      print(paste("Error: only 1 tip found in the tree for genus", genus))
      next
    }

    # Find the MRCA of all tips in the genus
    mrca <- getMRCA(list_of_trees_all_genes[[q]], tips)
    print(paste("MRCA for genus", genus, ":", mrca))

    if (is.na(mrca)) {
      print(paste("Error: MRCA is NA for genus", genus))
      next
    } else {
      "MRCA is not NA"
    }

    # Check if all descendants of the MRCA are in the same genus
    descendants <- getDescendants(list_of_trees_all_genes[[q]], node = mrca)
    #print(paste("Descendants for genus", genus, ":", descendants))
    tip_labels <- unname(list_of_trees_all_genes[[q]]$tip.label[descendants[descendants <= length(list_of_trees_all_genes[[q]]$tip.label)]])


    if (any(is.na(descendants))) {
      print(paste("Error: NA in descendants for genus", genus))
      next
    } else {
      "No descendants are NA"
    }

    # Find all descendant tips
    # Split the labels by space and extract the genus names
    descendants_genera <- unique(sapply(tip_labels, function(x) strsplit(x, " ")[[1]][1]))
    if (length(unique(descendants_genera)) > 1){
      print(paste("Error: Descendants of MRCA of", genus, "are not all in the same genus"))
    } else {
      print(paste("Descendants of MRCA of", genus, "are all in the same genus"))
    }

    # Create a subtree which contains only the descendants of the MRCA
    subtree <- extract.clade(list_of_trees_all_genes[[q]], node = mrca)

    # Calculate the average local posterior probability for the genus
    local_posterior_probabilities <- subtree$node.label
    print(paste("Local posterior probabilities for genus", genus, ":", paste(local_posterior_probabilities, collapse=", ")))

    if (any(is.na(local_posterior_probabilities))) {
      print(paste("Error: NA in local posterior probabilities for genus", genus))
      next
    }

    average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)

    # Print the results
    print(paste("Genus:", genus))
    print(paste("Average local posterior probability:", average_local_posterior_probability))
    print(paste("Number of species:", length(tip_labels)))
    #print(paste("Species:", paste(figurename_idx[descendants], collapse=", ")))

    # Save the results
    output <- data.frame(genus = genus, average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))

    results <- rbind(results, output)
  }

  # Can I order the list aphabetically?
  results <- results[order(results$genus),]
  results

  # Calculate support values for Arenga and Livistona without the following species 
  # Arenga distincta & Arenga hastata

  # Arenga
  tips <- unname(figurename_idx[which(grepl("Arenga", figurename_idx))])
  tips <- as.vector(tips)
  tips

  # Are all the tips found in the tree?
  if(!all(tips %in% list_of_trees_all_genes[[q]]$tip.label)){
    print("Error: not all tips found in the tree for genus Arenga")

    # Which tips are not found in the tree?
    print(paste("Tips not found in the tree for genus Arenga:", paste(tips[!tips %in% list_of_trees_all_genes[[q]]$tip.label], collapse=", ")))

    # Removing these tips from the tips vector
    tips <- tips[tips %in% list_of_trees_all_genes[[q]]$tip.label]
  } else {
    print("All tips found in the tree for genus Arenga")
  }

  # I need to remove genera with one 1 tip as I cannot find their MRCA
  if (length(tips) == 1){
    print("Error: only 1 tip found in the tree for genus Arenga")
    next
  }

  # Remove Arenga distincta and Arenga hastata
  tips <- tips[!grepl("distincta", tips)]
  tips <- tips[!grepl("hastata", tips)]
  tips

  # Find the MRCA of all tips in the genus
  mrca <- getMRCA(list_of_trees_all_genes[[q]], tips)
  print(paste("MRCA for genus Arenga:", mrca))

  if (is.na(mrca)) {
    print("Error: MRCA is NA for genus Arenga")
  } else {
    "MRCA is not NA"
  }

  # Check if all descendants of the MRCA are in the same genus
  descendants <- getDescendants(list_of_trees_all_genes[[q]], node = mrca)
  tip_labels <- unname(list_of_trees_all_genes[[q]]$tip.label[descendants[descendants <= length(list_of_trees_all_genes[[q]]$tip.label)]])
  tip_labels

  # Create subtree
  subtree <- extract.clade(list_of_trees_all_genes[[q]], node = mrca)

  # Calculate the average local posterior probability for the genus
  local_posterior_probabilities <- subtree$node.label

  # find mean
  average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
  print(paste("Average local posterior probability for Arenga:", average_local_posterior_probability))

  # and the number of species
  print(paste("Number of species for Arenga:", length(subtree$tip.label)))
  length(tip_labels)

  # add the new results for Arenga to the results dataframe
  output <- data.frame(genus = "Arenga", average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))
  # Update the results dataframe
  results[results$genus == "Arenga",] <- output
  results

  # Livistona
  # Livistona exigua
  tips <- unname(figurename_idx[which(grepl("Livistona", figurename_idx))])
  tips <- as.vector(tips)
  tips

  # Are all the tips found in the tree?
  if(!all(tips %in% list_of_trees_all_genes[[q]]$tip.label)){
    print("Error: not all tips found in the tree for genus Arenga")

    # Which tips are not found in the tree?
    print(paste("Tips not found in the tree for genus Arenga:", paste(tips[!tips %in% list_of_trees_all_genes[[q]]$tip.label], collapse=", ")))

    # Removing these tips from the tips vector
    tips <- tips[tips %in% list_of_trees_all_genes[[q]]$tip.label]
  } else {
    print("All tips found in the tree for genus Arenga")
  }

  # I need to remove genera with one 1 tip as I cannot find their MRCA
  if (length(tips) == 1){
    print("Error: only 1 tip found in the tree for genus Arenga")
    next
  }

  # Remove Arenga distincta and Arenga hastata
  tips <- tips[!grepl("exigua", tips)]
  tips

  # Find the MRCA of all tips in the genus
  mrca <- getMRCA(list_of_trees_all_genes[[q]], tips)
  print(paste("MRCA for genus Livistona:", mrca))

  if (is.na(mrca)) {
    print("Error: MRCA is NA for genus Livistona")
  } else {
    "MRCA is not NA"
  }

  # Check if all descendants of the MRCA are in the same genus
  descendants <- getDescendants(list_of_trees_all_genes[[q]], node = mrca)

  tip_labels <- unname(list_of_trees_all_genes[[q]]$tip.label[descendants[descendants <= length(list_of_trees_all_genes[[q]]$tip.label)]])

  # Create subtree
  subtree <- extract.clade(list_of_trees_all_genes[[q]], node = mrca)

  # Calculate the average local posterior probability for the genus
  local_posterior_probabilities <- subtree$node.label

  # find mean
  average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
  print(paste("Average local posterior probability for Livistona:", average_local_posterior_probability))

  # and the number of species
  print(paste("Number of species for Livistona:", length(tip_labels)))
  length(tip_labels)
  tip_labels

  # add the new results for Livistona to the results dataframe
  output <- data.frame(genus = "Livistona", average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))
  # Update the results dataframe
  results[results$genus == "Livistona",] <- output
  results


  # Saribus
  # Saribus brevifolius
  tips <- unname(figurename_idx[which(grepl("Saribus", figurename_idx))])
  tips <- as.vector(tips)
  tips

  # Are all the tips found in the tree?
  if(!all(tips %in% list_of_trees_all_genes[[q]]$tip.label)){
    print("Error: not all tips found in the tree for genus Saribus")

    # Which tips are not found in the tree?
    print(paste("Tips not found in the tree for genus Saribus:", paste(tips[!tips %in% list_of_trees_all_genes[[q]]$tip.label], collapse=", ")))

    # Removing these tips from the tips vector
    tips <- tips[tips %in% list_of_trees_all_genes[[q]]$tip.label]
  } else {
    print("All tips found in the tree for genus Saribus")
  }

  # I need to remove genera with one 1 tip as I cannot find their MRCA
  if (length(tips) == 1){
    print("Error: only 1 tip found in the tree for genus Saribus")
    next
  }

  # Remove Saribus brevifolius
  tips <- tips[!grepl("brevifolius", tips)]
  tips

  # Find the MRCA of all tips in the genus
  mrca <- getMRCA(list_of_trees_all_genes[[q]], tips)
  print(paste("MRCA for genus Saribus:", mrca))

  if (is.na(mrca)) {
    print("Error: MRCA is NA for genus Saribus")
  } else {
    "MRCA is not NA"
  }

  # Check if all descendants of the MRCA are in the same genus
  descendants <- getDescendants(list_of_trees_all_genes[[q]], node = mrca)

  tip_labels <- unname(list_of_trees_all_genes[[q]]$tip.label[descendants[descendants <= length(list_of_trees_all_genes[[q]]$tip.label)]])

  # Create subtree
  subtree <- extract.clade(list_of_trees_all_genes[[q]], node = mrca)

  # Calculate the average local posterior probability for the genus
  local_posterior_probabilities <- subtree$node.label

  # find mean
  average_local_posterior_probability <- mean(as.numeric(local_posterior_probabilities), na.rm=TRUE)
  print(paste("Average local posterior probability for Saribus:", average_local_posterior_probability))

  # and the number of species
  print(paste("Number of species for Saribus:", length(tip_labels)))
  length(tip_labels)
  tip_labels

  # add the new results for Saribus to the results dataframe
  output <- data.frame(genus = "Saribus", average_local_posterior_probability = average_local_posterior_probability, number_of_species = length(tip_labels))
  # Update the results dataframe
  results[results$genus == "Saribus",] <- output
  results


  # Can I make a for loop for the genera which finds the support value for the missing genera

  # Finding the genera which are NOT in results$genus
  miss_genera <- genera[!genera %in% results$genus]
  miss_genera


  for (genus in miss_genera){
    print(paste("Genus:", genus))

    tips <- unname(figurename_idx[which(grepl(genus, figurename_idx))])
    tips <- as.vector(tips)
    tips

    # are all the tips found in the tree?
    if(!all(tips %in% list_of_trees_all_genes[[q]]$tip.label)){
      print(paste("Error: not all tips found in the tree for genus", genus))

      # Which tips are not found in the tree?
      print(paste("Tips not found in the tree for genus", genus, ":", paste(tips[!tips %in% list_of_trees_all_genes[[q]]$tip.label], collapse=", ")))

      # Removing these tips from the tips vector
      tips <- tips[tips %in% list_of_trees_all_genes[[q]]$tip.label]
    } else {
      print(paste("All tips found in the tree for genus", genus))
    }

    # Find the tip index
    tip_index <- unname(which(list_of_trees_all_genes[[q]]$tip == tips))
    tip_index

    # Find the parent node
    parent_node <- list_of_trees_all_genes[[q]]$edge[list_of_trees_all_genes[[q]]$edge[,2] == tip_index, 1]
    parent_node <- parent_node - length(list_of_trees_all_genes[[q]]$tip.label)
    parent_node

    lpp_val <- list_of_trees_all_genes[[q]]$node.label[parent_node]


    # Print the results
    print(paste("Genus:", genus))
    print(paste("Average local posterior probability:", lpp_val))
    print(paste("Number of species:", length(tips)))
    #print(paste("Species:", paste(figurename_idx[descendants], collapse=", ")))

    # Save the results
    output <- data.frame(genus = genus, average_local_posterior_probability = lpp_val, number_of_species = length(tips))

    results <- rbind(results, output)
  }

  # # Adding the missing species
  # miss_genera <- data.frame(genus = miss_genera, average_local_posterior_probability = NA, number_of_species = 1)
  # results <- rbind(results, miss_genera)

  results <- results[order(results$genus),]
  list_of_results_all_genes[[q]] <- results
}

# Now I need to calculate the average support values for each genus in each of the 100 trees.
# I will do this by calculating the mean of the average local posterior probabilities for each genus in each of the 100 trees.
# I will also calculate the standard deviation of the average local posterior probabilities for each genus in each of the 100 trees.
average_local_across_100_trees_all_genes <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(average_local_across_100_trees_all_genes) <- c("genus", "mean_average_lpp", "sd_average_lpp", "number_of_species")

for (i in 1:length(genera)) {
  genus <- genera[i]
  print(paste("Genus:", genus))
  mean_average_lpp <- mean(sapply(list_of_results_all_genes, function(x) x$average_local_posterior_probability[x$genus == genus]))
  sd_average_lpp <- sd(sapply(list_of_results_all_genes, function(x) x$average_local_posterior_probability[x$genus == genus]))
  number_of_species <- mean(sapply(list_of_results_all_genes, function(x) x$number_of_species[x$genus == genus]))
  average_local_across_100_trees_all_genes <- rbind(average_local_across_100_trees_all_genes, data.frame(genus = genus, mean_average_lpp = mean_average_lpp, sd_average_lpp = sd_average_lpp, number_of_species = number_of_species))
}

average_local_across_100_trees_all_genes

# write csv file
write.csv(average_local_across_100_trees_all_genes, file = "/home/au543206/Documents/Coryphoideae/Figures_in_r/data/average_lpp_per_genera_all_genes_100_trees.csv", row.names = FALSE)



# Can we check wether there is a linear correlation between the number of species in a genus and the average local posterior probability?
# We can do this for both the ortholog tree and the all genes tree

# Ortholog tree
results_orthologs <- read.csv("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/average_lpp_per_genera_orthologs_100_trees.csv")

# All genes tree
results_all_genes <- read.csv("/home/au543206/Documents/Coryphoideae/Figures_in_r/data/average_lpp_per_genera_all_genes_100_trees.csv")


results_orthologs

# Ortholog tree
cor.test(results_orthologs$number_of_species, results_orthologs$mean_average_lpp)
# get estimate  R² value
correlation_coefficient_orthologs <- cor.test(results_orthologs$number_of_species, results_orthologs$mean_average_lpp)$estimate
r_squared_orthologs <- correlation_coefficient_orthologs^2
print(r_squared_orthologs)

# All genes tree
cor.test(results_all_genes$number_of_species, results_all_genes$mean_average_lpp)

# get estimate  R² value
correlation_coefficient_all_genes <- cor.test(results_all_genes$number_of_species, results_all_genes$mean_average_lpp)$estimate
r_squared_all_genes <- correlation_coefficient_all_genes^2
print(r_squared_all_genes)


# Lets also test wether there is a correlation with the total number of species in the genus and the average local posterior probability
# First we need to find the total number of species in each genus, we can do this with the wcvp.

# Load the wcvp
path_to_wcvp_coryphoideae <- "/home/au543206/Documents/world_checklist/World_checklist_downloads/05_09_2024/coryphoideae_accepted.csv" # This file is created by the finding_wcvp_coryphoideae.r script

# read csv using fread from data.table
wcvp <- data.table::fread(path_to_wcvp_coryphoideae)

# Only keep rows where taxon_rank == Species.
wcvp <- wcvp[which(wcvp$taxon_rank == "Species"),]

# Remove rows where hybrid_formula is != from ""
wcvp <- wcvp[which(wcvp$hybrid_formula == ""),]

# Now we count the number of species in each genus in wcvp.
number_of_species_per_genus <- wcvp[, .(number_of_species = .N), by = .(genus)]
sp_per_genus <- as.data.frame(table(wcvp$genus))

# We have to do some manual editing because we have joined wallichia and arenga
sp_per_genus[sp_per_genus$Var1 == "Arenga",2] <- (sp_per_genus[sp_per_genus$Var1 == "Arenga",2] + sp_per_genus[sp_per_genus$Var1 == "Wallichia",2])

# Can we now join these number of sp per genus to the results dataframes?
results_orthologs <- merge(results_orthologs, sp_per_genus, by.x = "genus", by.y = "Var1")
results_all_genes <- merge(results_all_genes, sp_per_genus, by.x = "genus", by.y = "Var1")

# Can we rename the new column to make it easier to understand?
results_orthologs <- results_orthologs %>% rename(total_number_of_species = Freq)
results_all_genes <- results_all_genes %>% rename(total_number_of_species = Freq)

# and now we add a new column to both datasets which is the proportion of species sampled.
results_orthologs$proportion_sampled <- results_orthologs$number_of_species / results_orthologs$total_number_of_species
results_all_genes$proportion_sampled <- results_all_genes$number_of_species / results_all_genes$total_number_of_species

# And now we will check for a significant correlation for both the total number of species and the proportion of species sampled.
# Ortholog tree
cor.test(results_orthologs$total_number_of_species, results_orthologs$mean_average_lpp) # Significant
cor.test(results_orthologs$proportion_sampled, results_orthologs$mean_average_lpp) # Not significant

# All genes tree
cor.test(results_all_genes$total_number_of_species, results_all_genes$mean_average_lpp) # Significant
cor.test(results_all_genes$proportion_sampled, results_all_genes$mean_average_lpp) # Not significant
