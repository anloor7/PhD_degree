

# Applying hierarchical clustering with QAF and DTW distances to the dataset Arabic Digits

# Loading the dataset standing for test set of spoken Arabic Digits

arabic <- read.csv('test_arabic_digits.txt', header = F, sep = ' ')
column1 <- arabic$V1
indexes <- which(column1 == 'a')

# We have to obtain 2200 series with 13 variables each one, and different lengths 

n <- 2200
X <- list(n)

for (i in 1 : (n-1)) {
  X[[i]] <- arabic[indexes[i] : indexes[i + 1],]
}
X[[n]] <- arabic[(indexes[n]) : nrow(arabic),]

# Removing useless lines 

for (i in 1 : n) {
  useless <- which(X[[i]]$V1 == 'a')
  X[[i]] <- X[[i]][-useless,]
  rownames(X[[i]]) <- NULL
  colnames(X[[i]]) <- NULL
  X[[i]] <- apply(as.matrix(X[[i]]), c(1,2), as.numeric)
}

# Ground truth 

ground_truth <- numeric(2200)

for (i in 1 : 10) {
  ground_truth[ ((i - 1)*220 + 1) : (220*i)] <- i
}


# As hierarchical clustering is computationally expensive, we are going to choose 10 elements of each one of 10 clusters

indexes <- numeric()

for (i in 0 : 9) {
  indexes <- c(indexes, (i*230+1) : (i*230+10))
}

M <- X[indexes]

ground_truth <- numeric(100)

for (i in 1 : 10) {
  ground_truth[ ((i - 1)*10 + 1) : (10*i)] <- i
}

# Then, we are going to perform dimensionality reduction via CPCA, as 13 variables are a lot of variables

sigma <- list()

for (i in 1 : 100) {
  sigma[[i]] <- cov(M[[i]])
}

s <- cpca(sigma, lambda = 0.80)
M_reduced <- list()

for (i in 1 : 100) {
 M_reduced[[i]] <- M[[i]] %*% s
}


# Hybrid approach with hierarchical clustering (numbers in the denominator are the 0.95 quantiles of the corresponding
# distances)

d_matrix <- matrix(0, 100, 100)

for (i in 1 : 100) {
  for (j in 1 : 100) {
    a <- qaf_mts_coefs_xy_sep(M_reduced[[i]])
    b <- qaf_mts_coefs_xy_sep(M_reduced[[j]])
    c <- EuclideanDistance(a, b)
    d <- dtw_mts(M_reduced[[i]], M_reduced[[j]])
    d_matrix[i, j] <- 0.5*(d/317) + 0.5*(c/0.93)# 0.5*(d/317) + 0.5*(c/1.3) # 0.5*(c/1.61)  +  0.5*(d/528.28)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 10)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering
