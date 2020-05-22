
# Applying hierarchical clustering with QAF and DTW distances to the dataset Arabic Digits

# Loading the dataset

japonese <- read.csv('jap_vowels.txt', header = F, sep = ' ', 
                     blank.lines.skip = FALSE)
column1 <- japonese$V1
indexes <- which(is.na(column1))
indexes <- indexes[-2]

# We have to obtain 370 series with 12 variables each one, and different lengths 

n <- 370
X <- list(n)

for (i in 1 : (n-1)) {
  X[[i]] <- japonese[indexes[i] : indexes[i + 1],]
}
X[[n]] <- japonese[(indexes[n]) : nrow(japonese),]

# Removing useless lines 

for (i in 1 : n) {
  useless <- which(is.na(X[[i]]$V1))
  X[[i]] <- X[[i]][-useless,]
  X[[i]]$V13 <- NULL
  rownames(X[[i]]) <- NULL
  colnames(X[[i]]) <- NULL
  X[[i]] <- apply(as.matrix(X[[i]]), c(1,2), as.numeric)
}

save(X, file = 'japonese_vowels.RData')

# Ground truth 

ground_truth_japonese_vowels <- c(rep(1, 31), rep(2, 35), rep(3, 88), rep(4, 44), rep(5, 29), rep(6, 24),
                  rep(7, 40), rep(8, 50), rep(9, 29))


# Dimensionaly reduction 

sigma <- list()

for (i in 1 : 370) {
  sigma[[i]] <- cov(X[[i]])
}

s <- cpca(sigma, lambda = 0.85)
japonese_reduced <- list()

for (i in 1 : 370) {
  japonese_reduced[[i]] <- X[[i]] %*% s
}


# Hybrid approach with hierarchical clustering

d_matrix <- matrix(0, 370, 370)

for (i in 1 : 370) {
  for (j in 1 : 370) {
   a <- qaf_mts_coefs_xy_sep(japonese_reduced[[i]])
   b <- qaf_mts_coefs_xy_sep(japonese_reduced[[j]])
   c <- EuclideanDistance(a, b)
   d <- dtw_mts(japonese_reduced[[i]], japonese_reduced[[j]])
    d_matrix[i, j] <- 0.5*(c/1.881) + 0.5*(d/95.31)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 9)
external_validation(ground_truth_japonese_vowels, clustering) # External validation of hierarchical clustering

# K means only with QAF 

gamma <- listTomatrix(lapply(X, qaf_mts_coefs_xy_sep))
clustering <- kmeans(gamma, 9)$cluster
external_validation(ground_truth_japonese_vowels, clustering)
