

# Applying hierarchical clustering with QAF and DTW distances to the dataset Arabic Digits

# Loading the dataset standing for a part of dataset AUSLAN

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/data_sets/raw_datasets/auslan')
temp = list.files(pattern="*.tsd")
auslan = lapply(temp, read.delim) # 285 MTSs

# Ground truth 

ground_truth <- numeric(285)

for (i in 1 : 95) {
  ground_truth[ ((i - 1)*3 + 1) : (3*i)] <- i
}

# As the dataset has series of different length, we can not apply k-means 

# As each MTS has 22 variables, we are going to perform dimensionality reduction

# Dimensionality reduction via CPCA

sigma <- list()

for (i in 1 : 285) {
  sigma[[i]] <- cov(auslan[[i]])
}

s <- cpca(sigma, lambda = 0.70)
auslan_reduced <- list()

for (i in 1 : 285) {
  auslan_reduced[[i]] <- as.matrix(auslan[[i]]) %*% s
}


# Hybrid approach with hierarchical clustering

d_matrix <- matrix(0, 285, 285)

for (i in 1 : 285) {
  for (j in 1 : 285) {
    a <- qaf_mts_coefs_xy_sep(auslan_reduced[[i]])
    b <- qaf_mts_coefs_xy_sep(auslan_reduced[[j]])
    c <- EuclideanDistance(a, b)
    d <- dtw_mts(auslan_reduced[[i]], auslan_reduced[[j]])
    d_matrix[i, j] <- d# 0.5*(d/317) + 0.5*(c/1.3) # 0.5*(c/1.61)  +  0.5*(d/528.28)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 95)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering


# QAF applied to the raw dataset 

gamma <- listTomatrix(lapply(auslan, qaf_mts_coefs_xy_sep))

