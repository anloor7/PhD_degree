


# Applying k-means with DTW distance to the dataset Character Trajectories 

data('uciCT')
ct <- CharTrajMV # Loading the data 

ground_truth <- numeric()
for (i in  1 : 20) { # Ground truth 
  ground_truth <- c(ground_truth, rep(i, 5))
}

# Removing the rows with all-zero values from each MTS

for (i in 1 : 100) {
  ct[[i]] <- ct[[i]][as.logical(rowSums(ct[[i]] != 0)), ]
  colnames(ct[[i]]) <- NULL
}

# Performing hierarchical clustering

d_matrix <- matrix(0, 100, 100)

for (i in 1 : 100) {
  for (j in 1 : 100) {
    d_matrix[i, j] <- dtw_mts(ct[[i]], ct[[j]])
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 20)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering 
