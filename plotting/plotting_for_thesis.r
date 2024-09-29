library(data.table)

# In this script I want to create some figures for my thesis
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)

  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

# First figure Is a pie chart showing the proportion of species assigned to different biomes according to the WCVP.

# Load the wcvp data
setwd("/home/au543206/Documents/world_checklist/World_checklist_downloads/05_09_2024")

wcvp <- fread("wcvp_names.csv")

head(wcvp)

# Filter the wcvp so we only have rows where taxon_status = Accepted, taxon_rank = Species 
wcvp <- wcvp[taxon_status == "Accepted" & taxon_rank == "Species"]

# Count the number of species in each biome
biome_counts <- wcvp[, .N, by = climate_description]



# Change "" to "Unknown"
biome_counts[climate_description == "", climate_description := "Unknown"]

# Make the first letter of each biome uppercase
biome_counts[, climate_description := tolower(climate_description)]
biome_counts[, climate_description := paste0(toupper(substr(climate_description, 1, 1)), substr(climate_description, 2, nchar(climate_description)))]

biome_counts

# Remove the subtropical or tropical ones as there are only 26
biome_counts <- biome_counts[climate_description != "Subtropical or tropical"]

# plot the pie chart
library(ggplot2)
library(MetBrewer)

# Set the colour palette to be vangogh2
pal <- met.brewer(name="Redon", n=12, type="discrete")

# Select colours 1-6 and 10-12
pal <- pal[c(1:4, 10:12)]
pal

# Adding grey to the palette
pal <- c(pal, "lightgrey")

# Can I set specific colours for each biome?
# print the colours

# Creating a list which assigns each biome to a colour
biome_colours <- list (
              "Temperate" = pal[1],
           "Wet tropical" = pal[4],
            "Subtropical" = pal[3],
"Seasonally dry tropical" = pal[6],
"Desert or dry shrubland" = pal[7],
 "Subalpine or subarctic" = pal[2],
       "Montane tropical" = pal[5],
            "Unknown" = pal[8] 
)

# Plot the pie chart
circle_sp <- ggplot(biome_counts, aes(x = "", y = N, fill = climate_description)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = biome_colours) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(fill = "Biome") +
  guides(fill=guide_legend(ncol=3, position = "bottom")) +
  ggtitle("Species richness")

circle_sp

# Now I want to plot the Koppen-geiger biomes ( I think)
koppen_path <- "/home/au543206/Documents/koppen_biomes/1991_2020/"
koppen_file <- "koppen_geiger_0p1.tif"

# Load the package for working with raster files
library(terra)
library(dplyr)

# Load the raster file
koppen_raster <- rast(file.path(koppen_path, koppen_file))

# Plot the raster file
plot(koppen_raster)

# Define the biome categories
biome_mapping <- list(
  "Temperate" = c(10, 12, 13, 15, 16, 17, 18, 21, 22, 25, 26),
  "Wet tropical" = c(1, 2),
  "Subtropical" = c(8, 9, 11, 14),
  "Seasonally dry tropical" = c(3),
  "Desert or dry shrubland" = c(4, 5, 6, 7),
  "Subalpine or subarctic" = c(19, 20, 23, 24, 27, 28, 29, 30),
  "Montane tropical" = c(),
  "Unknown" = c()
)

# Create a data frame with numbers 1 to 30
biome_df <- data.frame(Number = 1:30)

# Assign biome names based on the mapping
biome_df$Biome <- NA  # Initialize the Biome column with NA

for (biome in names(biome_mapping)) {
  for (number in biome_mapping[[biome]]) {
    biome_df$Biome[biome_df$Number == number] <- biome
  }
}

levels(koppen_raster) <- biome_df

plot(koppen_raster)

# plot the spatial raster in some sort of ggplot way to show the different biomes with the same colours as the pie chart
# I think I need to convert the raster to a data frame first
raster_df <- as.data.frame(koppen_raster, xy = TRUE)
head(raster_df)

# Plot using ggplot
ggplot_obj <- ggplot(raster_df, aes(x = x, y = y, fill = Biome)) +  # Use 'Biome' for fill
  geom_raster() +
  scale_fill_manual(values = biome_colours) +
  theme_void() +
  theme(legend.position = "None") +
  labs(fill = "Biome")

# Print the ggplot object
print(ggplot_obj)

##################################################################################################

# Define your equal-area projection (e.g., cylindrical Equal-Area)
equal_area_proj <- "+proj=cea"

# Reproject the raster to equal-area projection
koppen_raster_equal_area <- project(koppen_raster, equal_area_proj, "near")

# Assuming you have your raster as 'koppen_raster'
# 1. Calculate the cell sizes
cell_sizes <- cellSize(koppen_raster_equal_area)

# 2. Get the values of the raster (biome classifications)
biome_values <- values(koppen_raster_equal_area$Biome)
biome_values

# 3. Create a data frame with biome values and their corresponding cell sizes
biome_area_df <- data.frame(Biome = biome_values, Area = values(cell_sizes))

# Remove rows with NA values (in case there are any)
biome_area_df <- na.omit(biome_area_df)

head(biome_area_df)

# 4. Aggregate the areas by biome
total_area_by_biome <- aggregate(area ~ Biome, data = biome_area_df, FUN = sum)

total_area_by_biome

# Assuming total_area_by_biome and biome_df are already defined

# Make sure total_area_by_biome is a data frame (if it isn't already)
total_area_by_biome <- as.data.frame(total_area_by_biome)

# Rename the columns for clarity
colnames(total_area_by_biome) <- c("Number", "Area")

# Merge the data frames on the biome number
merged_biome_data <- merge(total_area_by_biome, biome_df, by = "Number")

# Reorder the columns for better readability
merged_biome_data <- merged_biome_data[, c("Biome", "Area")]

# Display the final result
print(merged_biome_data)

# Assuming merged_biome_data is the data frame with Biome and Area

# Group by Biome and summarize the total area for each biome
total_area_by_biome_summary <- merged_biome_data %>%
  group_by(Biome) %>%
  summarize(Total_Area = sum(Area, na.rm = TRUE))

# Display the summarized data
print(total_area_by_biome_summary)

# Calculate the proportion of each biome based on the total area
total_area_by_biome_summary$Proportion <- total_area_by_biome_summary$Total_Area / sum(total_area_by_biome_summary$Total_Area)

total_area_by_biome_summary


# Create the pie chart
circle_area <- ggplot(total_area_by_biome_summary, aes(x = "", y = Proportion, fill = Biome)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = biome_colours) +
  theme_void() +
  theme(legend.position = "none") +
  labs(fill = "Biome") +
  guides(fill=guide_legend(ncol=2)) +
  ggtitle("Area proportions")

circle_area

###################################################################################################
library(cowplot)

circle <- plot_grid(circle_area, circle_sp + theme(legend.position="none") , ncol = 1, rel_widths = c(1, 1))
circle


# extract a legend that is laid out horizontally
legend_b <- get_legend(circle_sp)

legend_b

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
circle_leg <- plot_grid(circle, legend_b, ncol = 1, rel_heights = c(2, .4))

circle_leg

# Combine the two plots
combined_plot <- plot_grid(ggplot_obj, circle_leg, align = "h", axis = "bt", ncol = 2, rel_widths = c(3, 2))

# Display the combined plot
print(combined_plot)
