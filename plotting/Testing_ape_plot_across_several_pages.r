# Load the ape package
library(ape)

# Example phylogenetic tree (replace with your actual tree)
tree <- rtree(10) 

# Number of pages to split the plot into
num_pages <- 3 

# Calculate the number of trees per page
trees_per_page <- ceiling(length(tree) / num_pages)

# Create a list to store the plots
plot_list <- list()

# Loop through each page
for (i in 1:num_pages) {
  # Determine the start and end indices for the current page
  start_idx <- (i - 1) * trees_per_page + 1
  end_idx <- min(i * trees_per_page, length(tree))

  # Create a new plotting device for each page
  pdf(paste0("phylogeny_page_", i, ".pdf"))

  # Loop through the trees on the current page
  for (j in start_idx:end_idx) {
    # Check for missing values (optional)
    if (any(is.na(tree[[j]]))) {
      warning("Tree", j, "contains missing values. Consider handling them before plotting.")
      next  # Skip plotting this tree if missing values are found (optional)
    }
    
    # Check for infinite values (optional)
    if (any(Inf %in% tree[[j]] | -Inf %in% tree[[j]])) {
      warning("Tree", j, "contains infinite values. Consider correcting them before plotting.")
      next  # Skip plotting this tree if infinite values are found (optional)
    }
    
    # Plot the current tree with (optional) default y-axis limits
    plot(tree[[j]], ylim = c(0, 1))  # Adjust limits as needed
  }

  # Close the plotting device
  dev.off()

  # Add the plot to the list
  plot_list[[i]] <- paste0("phylogeny_page_", i, ".pdf")
}

# Optionally, display the list of generated PDF files
print(plot_list)

plot(tree)


dev.off()
