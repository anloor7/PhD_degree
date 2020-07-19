

# Comparing QAF, QC1 and QC2. Scenario 8. 

n <- 20 # Number of series per cluster
l <- 400 # Length
K <- 4

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()

theta_c1 <- matrix(c(0.3, 0.3, 0.3, 0.0, 0.3, 0.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0, 0.3), nrow = 4)
sigma_c1 <- diag(4)
theta_c2 <- matrix(c(0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5), nrow = 4)
sigma_c2 <- sigma_c1
phi_c3 <- matrix(c(0.2, 0, 0.2, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0.2, 0, 0.2), nrow = 4)
sigma_c3 <- sigma_c1
phi_c4 <- matrix(c(0.1, 0, 0.1, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0, 0.1, 0, 0.1), nrow = 4)
sigma_c4 <- sigma_c1


set.seed(1234)
for (i in 1 : n) {
  cluster1[[i]] <- VARMAsim(l, malags = 1, theta = theta_c1, sigma = sigma_c1)$series
}

for (i in 1 : n) {
  cluster2[[i]] <- VARMAsim(l, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : n) {
  cluster3[[i]] <- VARMAsim(l, arlags = 1, phi = phi_c3,
                            sigma = sigma_c3)$series
}

for (i in 1 : n) {
  cluster4[[i]] <- VARMAsim(l, arlags = 1, phi = phi_c4,
                            sigma = sigma_c4)$series
}


cluster <- c(cluster1, cluster2, cluster3, cluster4)
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n))


# Parallelization 

n <- 7 # Number of cores 

# c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('cluster', 'qaf_mts_coefs_xy_sep', 'quantile_coherence', 'quantile_coherence_re_im', 'kmeans_mc_av_ari_scenario1'))


#  Performing k-means QAF 

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
gammal <- list(length = 20) 

for (i in 1 : 20) {
  gammal[[i]] <- gamma
}


b <- numeric()
b <- parLapply(c1, gammal, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

# Performing k-means quantile coherence (modulus)

coherence1 <- listTomatrix(lapply(cluster, quantile_coherence))
coherence1l <- list(length = 20) 

for (i in 1 : 20) {
  coherence1l[[i]] <- coherence1
}

b <- numeric()
b <- parLapply(c1, coherence1l, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

# Performing k-means quantile coherence (real-imaginary)

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)
coherence2l <- list(length = 20) 

for (i in 1 : 20) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

  
J <- 6 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
wf <- "d4"
features <- lapply(cluster, wavelet_features, wf = wf, J = J) 
dis_matrix <- proxy::dist(features, wave_dist)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# Alonso 

features <- listTomatrix(lapply(cluster, gcc_features_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)


