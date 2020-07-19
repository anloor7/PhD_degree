
n <- 10 # Number of series per cluster
l <- 1000 # Length 
K <- 3 # Number of clusters
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n))

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
sigma <- list()
epsilon <- list()

theta_c1 <- matrix(c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
theta_c2 <- matrix(c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), nrow = 3)
sigma_c2 <- sigma_c1
theta_c3 <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), nrow = 3)
sigma_c3 <- sigma_c1

set.seed(1234)


# Cluster 1



for (i in 1 : n) {
  
  sigma[[i]] <- diag(3)
  #epsilon[[i]] <- rmvt(l,  sigma = sigma1[[i]]) # Student
  epsilon[[i]] <- rdmsn(l, 3, mean = c(0, 0, 0), cov = sigma[[i]], del = c(4, 4, 4)) # Skewed normal 
  cluster1[[i]] <- varma(l, k = 3, VMA = theta_c1, innov = epsilon[[i]] )
  
}


# Cluster 2

for (i in 1 : n) {
  
  sigma[[i]] <- diag(3)
  epsilon[[i]] <- rmvt(l,  sigma = sigma1[[i]]) # Student
  # epsilon[[i]] <- rdmsn(l, 3, mean = c(0, 0, 0), cov = sigma[[i]], del = c(4, 4, 4)) # Skewed normal 
  cluster2[[i]] <- varma(l, k = 3, VMA = theta_c2, innov = epsilon[[i]])
  
}



# Cluster 3 

for (i in 1 : n) {
  
  sigma[[i]] <- diag(3)
  epsilon[[i]] <- rmvt(l,  sigma = sigma1[[i]]) # Student
  # epsilon[[i]] <- rdmsn(l, 3, mean = c(0, 0, 0), cov = sigma[[i]], del = c(4, 4, 4)) # Skewed normal 
  cluster3[[i]] <- varma(l, k = 3, VMA = theta_c3, innov = epsilon[[i]])
  
}

cluster <- c(cluster1, cluster2, cluster3)


# QAF

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
clustering <- kmeans(gamma, K)$cluster
external_validation(ground_truth, clustering)

# QC2

coherence <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence, EuclideanDistance)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)



# Wavelets

J <- 3 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
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
