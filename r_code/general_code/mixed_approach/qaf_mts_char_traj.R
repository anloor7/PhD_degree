

# Applying k-means clustering with QAF distance to the dataset Character Trajectories 

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


# Applying k-means

gamma <- listTomatrix(lapply(ct, qaf_mts_coefs_xy))
clustering <- kmeans(gamma, 20)$cluster

# Average adjusted rand index 

b <- numeric()

for (i in 1 : 3000) {
  clustering <- kmeans(gamma, 20)$cluster
  b <- external_validation(ground_truth, clustering)
}

mean(b)

# Fuzzy c-means

external_validation(fuzzytocrisp(t(fcm(gamma, 20)$u)), ground_truth)

# Applying hierarchical clustering 

d_matrix <- matrix(0, 100, 100)

for (i in 1 : 100) {
  for (j in 1 : 100) {
    a <- qaf_mts_coefs_xy(ct[[i]])
    b <- qaf_mts_coefs_xy(ct[[j]])
    d_matrix[i, j] <- EuclideanDistance(a, b)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 20)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering 


# Hybrid approach with hierarchical clustering 

d_matrix <- matrix(0, 100, 100)

for (i in 1 : 100) {
  for (j in 1 : 100) {
    a <- qaf_mts_coefs_xy(ct[[i]])
    b <- qaf_mts_coefs_xy(ct[[j]])
    c <- EuclideanDistance(a, b)
    d <- dtw_mts(ct[[i]], ct[[j]])
    d_matrix[i, j] <- 0.05*(c/1.31) + 0.95*(d/266.8)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 20)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering
