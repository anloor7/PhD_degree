

# Cluster creation

cluster1 <- list()
cluster2 <- list()
phi_c1 <- matrix(c(0.3, 0, 0.2, 0.2, 0.4, 0, 0, 0.3, 0.1), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c2 <- matrix(c(0.4, 0, 0.3, 0.2, 0.4, 0, 0, 0.4, 0.2), nrow = 3)
sigma_c2 <- sigma_c1


for (i in 1 : 5) {
  cluster1[[i]] <- VARMAsim(2000, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 5) {
  cluster2[[i]] <- VARMAsim(2000, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

cluster <- c(cluster1, cluster2)

# Ground truth 

ground_truth <- c(rep(1, 5), rep(2, 5))


# Distance matrix

dis_matrix <- matrix(0, 10, 10)
for (i in 1 : 10) {
  for (j in 1 : 10) {
    dis_matrix[i, j] <- j_divergence(cluster[[i]], cluster[[j]])
  }
}

# PAM 

clustering <- pam(dis_matrix, 2)$clustering
external_validation(ground_truth, clustering)


gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
clustering1 <- pam(gamma, 2)$cluster
external_validation(ground_truth, clustering1)
