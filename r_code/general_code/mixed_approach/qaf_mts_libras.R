

# Lets perfom clustering via multivariate QAF with Libras dataset


# Loading and preparing the data 

libras <- read.csv('libras.txt', header = F)
colnames(libras) <- NULL
S <- vector(mode = "list", length = nrow(libras))

for (i in 1 : length(S)) {
  S1 <- libras[i, c(seq(1, 90, by = 2), 91)] 
  S2 <- libras[i, c(seq(2, 90, by = 2), 91)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

M <- vector(mode = "list", length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- t(S[[i]][,seq(1, 45)])
}

# Ground truth 

ground_truth <- numeric(length(M))
for (j in 1 : length(M)) {
  ground_truth[j] <- S[[j]][1, 46]
}

# Applying k-means algorithm 

gamma <- listTomatrix(lapply(M, qaf_mts_coefs_xy))
clustering <- kmeans(gamma, 15)$cluster
external_validation(ground_truth, clustering)

# Hybrid approach with hierarchical clustering 

d_matrix <- matrix(0, 360, 360)

for (i in 1 : 360) {
  for (j in 1 : 360) {
   a <- qaf_mts_coefs_xy(M[[i]])
   b <- qaf_mts_coefs_xy(M[[j]])
   c <- EuclideanDistance(a, b)
   d <- dtw_mts(M[[i]], M[[j]])
   d_matrix[i, j] <- 0.1*(c/2.16) + 0.99*(d/60)
  }
}

hierarchical <- hclust(dist(d_matrix))
plot(hierarchical)
clustering <- cutree(hierarchical, 15)
external_validation(ground_truth, clustering) # External validation of hierarchical clustering


# Hybrid approach with kmeans 

hybrid_distance <- function(X, Y){
  
  a <- as.vector(qaf_mts_coefs_xy(X))
  b <- as.vector(qaf_mts_coefs_xy(Y))
  c <- EuclideanDistance(a, b)
  d <- dtw_mts(X, Y)
  0.1 * (c/2.16) + 0.9 * (d/60)
}

clustering <- km_mts(M, 15, dis = hybrid_distance)
external_validation(ground_truth, clustering)

