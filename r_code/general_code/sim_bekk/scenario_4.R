

# Comparing QAF, QC1 and QC2. Scenario 1 

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()


set.seed(1234)

for (i in 1 : 100) {
  simulated = simulateBEKK(3, 1500, c(1, 1), params = c(0.1, 0.1, 0.1, 0, 0, 0, 0.3, 0, 0, 0.4, 0, 0, 0.2, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0))
  series <- simulated$eps
  cluster1[[i]] <- cbind(as.vector(series[[1]]), as.vector(series[[2]]), as.vector(series[[3]]))
}

for (i in 1 : 100) {
  simulated = simulateBEKK(3, 1500, c(1, 1), params = c(0.1, 0.1, 0.1, 0, 0, 0, 0.7, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0,
                                                        0, 0, 0.2, 0, 0, 0, 0))
  series <- simulated$eps
  cluster2[[i]] <- cbind(as.vector(series[[1]]), as.vector(series[[2]]), as.vector(series[[3]]))
}

for (i in 1 : 100) {
  simulated = simulateBEKK(3, 1500, c(1, 1), params = c(0.1, 0.1, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5,
                                                        0, 0.3, 0, 0, 0, 0.1, 0))
  series <- simulated$eps
  cluster3[[i]] <- cbind(as.vector(series[[1]]), as.vector(series[[2]]), as.vector(series[[3]]))
}

for (i in 1 : 100) {
  simulated = simulateBEKK(3, 1500, c(1, 1), params = c(0.1, 0.1, 0.1, 0, -0.4, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0,
                                                        0.7, 0, 0, 0, 0, 0, 0))
  series <- simulated$eps
  cluster4[[i]] <- cbind(as.vector(series[[1]]), as.vector(series[[2]]), as.vector(series[[3]]))
}


cluster <- c(cluster1, cluster2, cluster3, cluster4)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))


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
coherence2l <- list(length = 20) 

for (i in 1 : 20) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario2)
max(listTomatrix(b))












