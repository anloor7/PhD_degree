

# Cluster creation

cluster1 <- list()
cluster2 <- list()
phi_c1 <- matrix(c(0.3, 0, 0.2, 0.2, 0.4, 0, 0, 0.3, 0.1), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c2 <- matrix(c(0.4, 0, 0.3, 0.2, 0.4, 0, 0, 0.4, 0.2), nrow = 3)
sigma_c2 <- sigma_c1

set.seed(1234)
for (i in 1 : 100) {
  cluster1[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster2[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

cluster <- c(cluster1, cluster2)

# Ground truth 

ground_truth <- c(rep(1, 100), rep(2, 100))


# Distance matrix

dis_matrix <- matrix(0, 200, 200)
for (i in 1 : 200) {
  for (j in 1 : 200) {
    k <- ncol(cluster1[[1]])
    d_matrix <- matrix(0, k, k)
    for (s in 1 : k)
      for (v in 1 : k){
        {
          d_matrix[s, v] <- EuclideanDistance(cluster[[i]][,s], cluster[[j]][,v]) 
        }
      }
    dis_matrix[i, j] <- dtw(d_matrix, distance.only = TRUE)$normalizedDistance
    # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'sakoechiba')$normalizedDistance
    # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'itakura')$normalizedDistance
  }
}


# PAM 

clustering <- pam(dis_matrix, 2)$clustering
clustering
external_validation(ground_truth, clustering)

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
clustering1 <- pam(gamma, 2, diss = F)$cluster
external_validation(ground_truth, clustering1)


