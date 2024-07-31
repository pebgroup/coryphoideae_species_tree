library(data.table)

datadir <- "/home/au543206/Documents/world_checklist/World_checklist_downloads/05_09_2024"

#Creating a list of all accepted names within Coryphoideae
data <- as.data.frame(fread(file.path(datadir,"wcvp_names.csv")))

species <- data[which(data$taxon_rank == "Species" | data$taxon_rank == "Variety" | data$taxon_rank == "Subspecies"),] # Selecting only species and varieties/subspecies
palms <- species[which(species$family=="Arecaceae"),] # Selecting only Palms
apalms <- palms[which(palms$taxon_status=="Accepted"),] # Selecting only Accepted species

# List of Genera in Coryphoideae
genera <- c("Sabal","Schippia","Trithrinax","Zombia","Coccothrinax","Hemithrinax","Leucothrinax","Thrinax","Chelyocarpus",
            "Cryosophila","Itaya","Phoenix","Chamaerops","Guihaia","Trachycarpus","Rhapidophyllum","Maxburretia","Rhapis",
            "Livistona","Licuala","Johannesteijsmannia","Pholidocarpus","Pritchardiopsis","Acoelorraphe","Serenoa","Brahea",
            "Colpothrinax","Copernicia","Pritchardia","Washingtonia","Chuniophoenix","Kerriodoxa","Nannorrhops","Tahina",
            "Caryota","Arenga","Wallichia","Corypha","Bismarckia","Satranala","Hyphaene","Medemia","Latania","Lodoicea",
            "Borassodendron","Borassus","Lanonia","Saribus","Sabinaria")


subfamilies <- c("Phoeniceae","Trachycarpeae","Sabaleae","Cryosophileae","Chuniophoeniceae","Caryoteae","Corypheae","Borasseae")


# Selecting only palms in Coryphoideae
cory_data <- apalms[which(apalms$genus %in% genera),]

#Dropping synonym species names rows in data.frame
cory_data_accepted <- cory_data[which(cory_data$taxon_status =="Accepted"),] # | (cory_data$taxon_rank=="Variety"))

dim(cory_data_accepted)
head(cory_data_accepted)

# Write the filtered data to a CSV file
output_file <- file.path(datadir, "coryphoideae_accepted.csv")
write.csv(cory_data_accepted, output_file, row.names = FALSE)

cat("Filtered data has been written to:", output_file)
