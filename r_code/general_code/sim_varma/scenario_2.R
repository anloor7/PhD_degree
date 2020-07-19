


# Comparing QAF, QC1 and QC2. Scenario 2. 

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()

phi_c1 <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
theta_c2 <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c3 <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
theta_c3 <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
sigma_c3 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)

set.seed(1234)
for (i in 1 : 50) {
  cluster1[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 50) {
  cluster2[[i]] <- VARMAsim(200, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 50) {
  cluster3[[i]] <- VARMAsim(200, arlags = 1, malags = 1, phi = phi_c3,
                             theta = theta_c3, sigma = sigma_c3)$series
}

cluster <- c(cluster1, cluster2, cluster3)
ground_truth <- c(rep(1, 50), rep(2, 50), rep(3, 50))



# Parallelization 

n <- 7 # Number of cores 

# c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('cluster', 'qaf_mts_coefs_xy_sep', 'quantile_coherence', 'quantile_coherence_re_im', 'kmeans_mc_av_ari_scenario1'))


#  Performing k-means QAF 

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
gammal <- list(length = 10) 

for (i in 1 : 10) {
  gammal[[i]] <- gamma
}


b <- numeric()
b <- parLapply(c1, gammal, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))

# Performing k-means quantile coherence (modulus)

coherence1 <- listTomatrix(lapply(cluster, quantile_coherence))
coherence1l <- list(length = 10) 

for (i in 1 : 10) {
  coherence1l[[i]] <- coherence1
}

b <- numeric()
b <- parLapply(c1, coherence1l, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))

# Performing k-means quantile coherence (real-imaginary)

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dist_matrix <- proxy::dist(coherence2, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

coherence2l <- list(length = 10) 

for (i in 1 : 10) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))