# optrimAl.R - run this second

cutoff_trim <- readLines('cutoff_trim.txt')

# create one multiple tables for each threshold value to store AMAS results
amas_table <- read.table('summary_0.txt', header = TRUE)
sites <- data.frame(row.names = amas_table$Alignment_name)

pct <- data.frame(row.names = amas_table$Alignment_name)
filled <- data.frame(row.names = amas_table$Alignment_name)
lost <- data.frame(row.names = amas_table$Alignment_name)

# Reading all the summary files and storing the data in the tables
for(i in 1:length(cutoff_trim)){
  amas_table <- read.table(paste('summary_', cutoff_trim[i], '.txt', sep = ''), header = TRUE)
  for(j in amas_table$Alignment_name){ # looping through alignment names
    sites[rownames(sites) == j,i] <- amas_table$Parsimony_informative_sites[amas_table$Alignment_name == j]
    pct[rownames(pct) == j,i] <- as.numeric(amas_table$Proportion_parsimony_informative[amas_table$Alignment_name == j])
    filled[rownames(filled) == j,i] <- amas_table$Total_matrix_cells[amas_table$Alignment_name == j] * (1 - amas_table$Missing_percent[amas_table$Alignment_name == j] / 100)
  }
}

# calculate data loss for each trimming threshold

sites[is.na(sites)] <- 0
pct[is.na(pct)] <- 0

for(i in 1:ncol(filled)){
  lost[,i] <- 1 - filled[,i] / filled[,1]
}

lost[is.na(lost)] <- 1

colnames(sites) <- cutoff_trim
colnames(pct) <- cutoff_trim
colnames(filled) <- cutoff_trim
colnames(lost) <- cutoff_trim


# select optimal trimming threshold
# current criterion is maximum proportion of parsimony informative sites where data loss is no more than one median absolute deviation above the median

optrim <- numeric()
optrim_loss <- numeric()

# I is the thrimming treshholds
for(i in rownames(pct)){
  lost_i <- unlist(lost[rownames(lost) == i, ])
  pct_i <- unlist(pct[rownames(pct) == i, ])
  dldp <- data.frame(pct_i, lost_i, row.names = cutoff_trim)
  write.csv(dldp, paste('dldp_', i, '.csv', sep = ''))

  # Finding the loss of the different trimming values 
  real_loss <- dldp$lost_i[dldp$lost_i < 1]
  if(length(real_loss) > 0) {
    diff_loss <- real_loss[2:length(real_loss)] - real_loss[1:(length(real_loss) - 1)]
    median_loss <- median(diff_loss[diff_loss != 0])
  } else {
    median_loss <- 0
  }
  dldp <- subset(dldp, dldp$lost_i <= (median(real_loss) + median_loss))

  if(length(dldp$pct_i) > 0){
    optrim[i] <- rownames(dldp)[dldp$pct_i == max(dldp$pct_i)][[1]]
    optrim_loss[i] <- dldp$lost_i[rownames(dldp) == optrim[i][[1]]]
  } else {
    optrim[i] <- 0
    optrim_loss[i] <- 0
  }
}


# generate graphs to show effect of trimming on informativeness and data loss
for(i in rownames(pct)){
	print(i)
  dldp <- read.csv(paste('dldp_', i, '.csv', sep = ''))

  print(paste("Length of dldp$lost_i:", length(dldp$lost_i)))
  print(paste("Length of cutoff_trim:", length(cutoff_trim)))

  
#   png(paste('dldp_', i, '.png', sep = ''))
#   par(mar = c(5,5,2,5))
#   plot(main = i, dldp$lost_i ~ cutoff_trim, ylim = c(0,1), ylab = 'proportion of data lost', xlab = 'strictness of trimming (trimAl gap threshold)', pch = 18, col = 'red')
#   par(new = T)
#   plot(dldp$pct_i ~ cutoff_trim, xlab = NA, ylab = NA, ylim = c(0,1), axes = F, pch = 16, col = 'blue')
#   axis(side = 4)
#   mtext(side = 4, line = 3, 'proportion parsimony informative')
#   legend(x = 0, y = 1, legend = c('proportion of data lost', 'proportion of parsimony informative sites', 'selected trimming threshold'), pch = c(18, 16, NA), lty = c(NA, NA, 2), col = c('red', 'blue', 'black'), cex = 0.9, bty = 'n')
#   if(is.na(optrim[i]) == FALSE){
#     lines(c(-0.5, optrim[i]), c(optrim_loss[i], optrim_loss[i]), lty = 2)
#     lines(c(-0.5, optrim[i]), c(dldp$pct_i[dldp$X == optrim[i]], dldp$pct_i[dldp$X == optrim[i]]), lty = 2)
#     lines(c(optrim[i], optrim[i]), c(-0.5, max(optrim_loss[i], dldp$pct_i[dldp$X == optrim[i]])), lty = 2)
#   }
#   dev.off()
# }

print("checkpoint")

overlost <- names(optrim_loss[optrim_loss > 0.3])

cat("The following sequences have more than 30% data loss at the optimal trimming threshold:\n")
cat(overlost, sep = '\n')
write(overlost, 'overlost.txt', sep = '\n')

cat("Copying the sequence with the optimal trimming value into the current directory.\n")
cat(paste(optrim, '/', names(optrim), sep = ''), sep = '\n')
#This line of code copies the sequence with the optimal trimming value into the current directory. 
#file.copy(paste(optrim, '/', names(optrim), sep = ''), getwd()) This line is commented out because it didnt work, potentially it is only because overwrite = TRUE isent given as a argument.

# file.remove(paste(overlost, sep = ''))

#### Testing other file copying code
# Loop over the files and move them from the respective subfolder
for (file in names(optrim)) {
  subfolder <- optrim[file]  # Get the corresponding subfolder
  source_path <- file.path(getwd(), paste0(subfolder), file)
  destination_path <- file.path(getwd(), file)
  
  # Copy the file and print the result
  success <- file.copy(source_path, overwrite = TRUE, destination_path)
  
  if (!success) {
    message("Failed to copy: ", source_path)
  } else {
    message("Copied: ", source_path, " -> ", destination_path)
  }
}


print("optrimal is done")