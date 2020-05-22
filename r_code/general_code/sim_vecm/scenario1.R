

# Comparing QAF, QC1 and QC2. Escenario 1

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()

B_c1 <- rbind(c(-0.2, 0.1, 0), c(0.2, 0, 0.1))
beta_c1 <- 1
B_c2 <- rbind(c(-0.2, 0.1, 0), c(0.2, 0, 0.1))
beta_c2 <- 0.8
B_c3 <- rbind(c(-0.2, 0.1, 0), c(0.2, 0, 0.1))
beta_c3 <- 0.6

set.seed(1234)
for (i in 1 : 100) {
  cluster1[[i]] <- VECM.sim(B = B_c1, beta = beta_c1, n = 100, lag = 1, include = "none")
}

for (i in 1 : 100) {
  cluster2[[i]] <- VECM.sim(B = B_c2, beta = beta_c2, n = 100, lag = 1, include = "none")
}

for (i in 1 : 100) {
  cluster3[[i]] <- VECM.sim(B = B_c3, beta = beta_c3, n = 100, lag = 1, include = "none")
}

cluster <- c(cluster1, cluster2, cluster3)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))


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
coherence2l <- list(length = 10) 

for (i in 1 : 10) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))

