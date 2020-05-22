

# Comparing QAF, QC1 and QC2. Scenario 6. 


cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()

phi_c1 <- 0.2 * diag(5)
phi_c1[2,1] <- 0.2
phi_c1[1,4] <- 0.2
phi_c1[3,5] <- 0.2
sigma_c1 <- matrix(c(1, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25,
                     0.25, 0.25, 0.25, 0.25, 1), nrow = 5)
phi_c2 <-phi_c1 
phi_c2[1, 4] <- 0.5
phi_c2[3, 5] <- 0.5
phi_c2[2, 1] <- 0.5
sigma_c2 <- sigma_c1
phi_c3 <-phi_c1
phi_c3[1, 4] <- 0.7
phi_c3[3, 5] <- 0.7
phi_c3[2, 1] <- 0.7
sigma_c3 <- sigma_c1
phi_c4 <- phi_c1
phi_c4[1, 4] <- 0.9
phi_c4[3, 5] <- 0.9
phi_c4[2, 1] <- 0.9
sigma_c4 <- sigma_c1

set.seed(1234)

for (i in 1 : 50) {
  cluster1[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 50) {
  cluster2[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

for (i in 1 : 50) {
  cluster3[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c3, sigma = sigma_c3)$series
}

for (i in 1 : 50) {
  cluster4[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c4, sigma = sigma_c4)$series
}

cluster <- c(cluster1, cluster2, cluster3, cluster4)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))


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
b <- parLapply(c1, gammal, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

# Performing k-means quantile coherence (modulus)

coherence1 <- listTomatrix(lapply(cluster, quantile_coherence))
coherence1l <- list(length = 10) 

for (i in 1 : 10) {
  coherence1l[[i]] <- coherence1
}

b <- numeric()
b <- parLapply(c1, coherence1l, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

# Performing k-means quantile coherence (real-imaginary)

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
coherence2l <- list(length = 10) 

for (i in 1 : 10) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))

