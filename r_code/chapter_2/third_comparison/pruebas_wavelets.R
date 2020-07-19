



n <- 5 # Number of series per cluster
l <- 1000 # Length 
K <- 2 # Number of clusters
ground_truth <- c(rep(1, n), rep(2, n))


cluster1 <- list()


for (j in 1 : n){
  
  garch1p <- garch.sim(c(0.01, 0.01), 0.95, rnd = rt(n = l+1, df = 2),  n = l + 1)
  garch2p <- garch.sim(c(0.01, 0.01), 0.95, rnd = rt(n = l+1, df = 3),  n = l + 1)
  
  
  cluster1[[j]] <- cbind(garch1p, garch2p)
  
}





cluster2 <- list()


for (j in 1 : n){
  
  garch1p <- garch.sim(c(0.01, 0.20), 0.65, n = l + 1)
  garch2p <- garch.sim(c(0.01, 0.20), 0.65, n = l + 1)
  
  
  cluster2[[j]] <- cbind(garch1p, garch2p)
  
}



cluster <- c(cluster1, cluster2)



# QC2

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# DM

J <- 5 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
wf <- "d4"

features <- lapply(cluster, wavelet_features, wf = wf, J = J) 
dis_matrix <- proxy::dist(features, wave_dist)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# AP

features <- listTomatrix(lapply(cluster, gcc_features_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)
