



# Comparing QAF, QC1 and QC2, scenario 1 

a1 <- garch.sim(alpha = c(0.01, 0.1), beta = 0.8, n = 8000)
b1 <- garch.sim(alpha = c(0.01, 0.8), beta = 0.1, n = 8000)
df1 <- cbind(a1, b1)
a2 <- garch.sim(alpha = c(0.01, 0.3), beta = 0.3, n = 8000)
b2 <- garch.sim(alpha = c(0.01, 0.4), beta = 0.5, n = 8000)
df2 <- cbind(a2, b2)
a3 <- garch.sim(alpha = c(0.01, 0.9), beta = 0.05, n = 8000)
b3 <- garch.sim(alpha = c(0.01, 0.1), beta = 0.1, n = 8000)
df3 <- cbind(a3, b3)

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()

A1 <- rbind(c(1, 1), c(0.1, 0.5))
K1 <- rbind(c(0, 0.5), c(0.5, 0.2))
spec1 <- gogarchspec(ica.fix = list(A = A1, K = K1))
fit1 <- gogarchfit(spec = spec1, data = df1)
A2 <- rbind(c(1, 1), c(0.1, 0.5))
K2 <- rbind(c(0, 0.5), c(0.5, 0.2))
spec2 <- gogarchspec(ica.fix = list(A = A2, K = K2))
fit2 <- gogarchfit(spec = spec2, data = df2)
A3 <- rbind(c(1, 1), c(0.1, 0.5))
K3 <- rbind(c(0, 0.5), c(0.5, 0.2))
spec3 <- gogarchspec(ica.fix = list(A = A3, K = K3))
fit3 <- gogarchfit(spec = spec3, data = df3)

set.seed(1234)


for (i in 1 : 100) {
  sim1 <- gogarchsim(fit1, n.sim = 1500)
  l1 <- slot(sim1, 'msim')
  cluster1[[i]] <- l1$seriesSim[[1]]
}

for (i in 1 : 100) {
  sim2 <- gogarchsim(fit2, n.sim = 1500)
  l2 <- slot(sim2, 'msim')
  cluster2[[i]] <- l2$seriesSim[[1]]
}

for (i in 1 : 100) {
  sim3 <- gogarchsim(fit3, n.sim = 1500)
  l3 <- slot(sim3, 'msim')
  cluster3[[i]] <- l3$seriesSim[[1]]
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

for (i in 1 : 20) {
  gammal[[i]] <- gamma
}


b <- numeric()
b <- parLapply(c1, gammal, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))

# Performing k-means quantile coherence (modulus)

coherence1 <- listTomatrix(lapply(cluster, quantile_coherence))
coherence1l <- list(length = 10) 

for (i in 1 : 20) {
  coherence1l[[i]] <- coherence1
}

b <- numeric()
b <- parLapply(c1, coherence1l, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))

# Performing k-means quantile coherence (real-imaginary)

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
coherence2l <- list(length = 10) 

for (i in 1 : 20) {
  coherence2l[[i]] <- coherence2
}

b <- numeric()
b <- parLapply(c1, coherence2l, kmeans_mc_av_ari_scenario1)
max(listTomatrix(b))
