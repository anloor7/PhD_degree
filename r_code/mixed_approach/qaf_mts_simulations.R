

# Lets perfom simulations in the second scenario of Baragona's paper, and applying hierarchical clustering with QAF

# Simulations

# Second scenario
# Cluster 1: VAR(1) model, 3 variables, 100 time points, 100 elements
# Cluster 2: VMA(1) model, 3 variables, 100 time points, 100 elements
# Cluster 3: VARMA(1) model, 3 variables, 100 time points, 100 elements

cluster12 <- list()
cluster22 <- list()
cluster32 <- list()

phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c3 <- matrix(c(0.5, -0.3, 0.0, -0.3, 0.5, -0.4, 0.0, 0.4, 0.5), nrow = 3)
theta_c3 <- matrix(c(-0.5, 0.4, 0.0, 0.4, -0.5, -0.3, 0.0, 0.3, -0.5), nrow = 3)
sigma_c3 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)

set.seed(1234)
for (i in 1 : 100) {
  cluster12[[i]] <- VARMAsim(30, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster22[[i]] <- VARMAsim(30, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 100) {
  cluster32[[i]] <- VARMAsim(30, arlags = 1, malags = 1, phi = phi_c3,
                             theta = theta_c3, sigma = sigma_c3)$series
}

experiment2 <- c(cluster12, cluster22, cluster32)
ground_truth2 <- c(rep(1, 100), rep(2, 100), rep(3, 100))


# Computing dissimilarity matrix

gamma <- listTomatrix(lapply(experiment2, qaf_mts_coefs))
dis_matrix <- matrix(0, 300, 300)
for (i in 1 : 300) {
  for (j in 1 : 300) {
    dis_matrix[i, j] <- EuclideanDistance(gamma[i,], gamma[j,])
  }
}

# Performing hierarchical clustering 

hierarchical <- hclust(dist(dis_matrix))
clustering <- cutree(hierarchical, 3)
plot(hierarchical)
external_validation(ground_truth, clustering)

# Performing k-means 

clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth, clustering)


# Computing dissimilarity matrix new approach


gamma <- listTomatrix(lapply(experiment2, qaf_mts_coefs_xy))
dis_matrix <- matrix(0, 300, 300)
for (i in 1 : 300) {
  for (j in 1 : 300) {
    dis_matrix[i, j] <- EuclideanDistance(gamma[i,], gamma[j,])
  }
}

# Performing hierarchical clustering new approach

hierarchical <- hclust(dist(dis_matrix))
clustering <- cutree(hierarchical, 3)
plot(hierarchical)
external_validation(ground_truth, clustering)

# Performing k-means new approach 

clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth, clustering)
