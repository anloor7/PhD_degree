

# Comparing QC2, DM y AP. Scenario 11. 

n <- 10 # Number of series per cluster
l <- 200 # Length
K <- 4
B <- 20 # Monte Carlo replicas
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n)) # Ground Truth 

set.seed(1234)

cluster <- list() 

for (k in 1 : B){
  
  cluster1 <- list()
  cluster2 <- list()
  cluster3 <- list()
  cluster4 <- list()
  
  phi_c1 <- matrix(c(0.5, 0, 0, 0, 0.3, 0.4, 0, 0, 0.3, 0, 0.5, 0, 0, 0, 0.2, 0.7), nrow = 4)
  theta_c1 <- matrix(c(0.2, 0, 0.5, 0.2, 0, 0.3, 0, 0.6, 0, 0.5, 0.4, 0, 0.2, 0.2, 0.5, 0.6), nrow = 4)
  sigma_c1 <- diag(4)
  phi_c2 <- matrix(c(0.5, 0, 0, 0, 0.3, -0.4, 0, 0, 0.3, 0, 0.5, 0, 0, 0, 0.2, 0.7), nrow = 4)
  theta_c2 <- matrix(c(0.2, 0, 0.5, 0.2, 0, 0.3, 0, 0.8, 0, 0.5, 0.4, 0, 0.2, 0.2, 0.8, 0.9), nrow = 4)
  sigma_c2 <- diag(4)
  phi_c3 <- matrix(c(0.5, 0, 0.2, 0, 0.3, -0.4, 0, 0, 0.3, 0, 0.5, 0, 0, 0, 0.2, 0.7), nrow = 4)
  theta_c3 <- matrix(c(0.2, 0, 0.5, 0.2, 0, 0.3, 0, 0.4, 0, 0.5, 0.4, 0, 0.2, 0.2, 0.3, 0.3), nrow = 4)
  sigma_c3 <- diag(4)
  phi_c4 <- matrix(c(0.5, 0, 0.2, 0, 0.3, -0.4, 0, 0, 0.3, 0, 0.5, 0, 0.3, 0, 0.5, 0.7), nrow = 4)
  theta_c4 <- matrix(c(0.2, 0, 0.5, 0.2, 0, 0.3, 0, 0.1, 0, 0.5, 0.4, 0, 0.2, 0.2, 0.1, 0.1), nrow = 4)
  sigma_c4 <- diag(4)
  sigma <- list()
  epsilon <- list()
  
  
  for (i in 1 : n) {
    
    sigma[[i]] <- diag(4)
    # epsilon[[i]] <- rmvt(l,  sigma = sigma[[i]], df = 3) # Student
    epsilon[[i]] <- rdmsn(l, mean = c(0, 0, 0, 0), cov = sigma[[i]], del = c(3, 0, 3, 0)) # Skewed normal 
    cluster1[[i]] <- varma(l, k = 4, VAR = phi_c1, VMA = theta_c1, innov = epsilon[[i]])
    
  }
  
  for (i in 1 : n) {
    
    sigma[[i]] <- diag(4)
    # epsilon[[i]] <- rmvt(l,  sigma = sigma[[i]], df = 3) # Student
    epsilon[[i]] <- rdmsn(l, 4, mean = c(0, 0, 0, 0), cov = sigma[[i]], del = c(3, 0, 3, 0)) # Skewed normal 
    cluster2[[i]] <- varma(l, k = 4, VAR = phi_c2, VMA = theta_c2, innov = epsilon[[i]])
    
  }
  
  for (i in 1 : n) {
    
    sigma[[i]] <- diag(4)
    # epsilon[[i]] <- rmvt(l,  sigma = sigma[[i]], df = 3) # Student
    epsilon[[i]] <- rdmsn(l, 4, mean = c(0, 0, 0, 0), cov = sigma[[i]], del = c(3, 0, 3, 0)) # Skewed normal 
    cluster3[[i]] <- varma(l, k = 4, VAR = phi_c3, VMA = theta_c3, innov = epsilon[[i]])
    
  }
  
  for (i in 1 : n) {
    
    
    sigma[[i]] <- diag(4)
    # epsilon[[i]] <- rmvt(l,  sigma = sigma[[i]], df = 3) # Student
    epsilon[[i]] <- rdmsn(l, 4, mean = c(0, 0, 0, 0), cov = sigma[[i]], del = c(3, 0, 3, 0)) # Skewed normal 
    cluster4[[i]] <- varma(l, k = 4, VAR = phi_c4, VMA = theta_c4, innov = epsilon[[i]])
    
  }
  
  
  cluster[[k]] <- c(cluster1, cluster2, cluster3, cluster4)
  
}




# Performing k-means quantile coherence (real-imaginary)


# coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
# dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
# clustering <- pam(dis_matrix, K)$cluster
# external_validation(ground_truth, clustering)
# coherence2l <- list(length = 20) 

ari_qc2 <- numeric()

for (i in 1 : B) {
  
  coherence2 <- listTomatrix(lapply(cluster[[i]], quantile_coherence_re_im))
  dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
  clustering <- pam(dis_matrix, K)$cluster
  ari_qc2[i] <- external_validation(ground_truth, clustering)
  
}

mean(ari_qc2)

# Wavelets 

J <- 6 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
wf <- "d4"

ari_dm <- numeric()

for (i in 1 : B) {
  
  features <- lapply(cluster[[i]], wavelet_features, wf = wf, J = J) 
  dis_matrix <- proxy::dist(features, wave_dist)  
  clustering <- pam(dis_matrix, K)$cluster
  ari_dm[i] <- external_validation(ground_truth, clustering)
  
}

mean(ari_dm)

# Alonso y Peña 

ari_ap <- numeric()

for (i in 1 : 3) {
  
  features <- listTomatrix(lapply(cluster[[i]], gcc_features_mts))
  dis_matrix <- proxy::dist(features, EuclideanDistance)
  clustering <- pam(dis_matrix, K)$cluster
  ari_ap[i] <- external_validation(ground_truth, clustering)
  
}

mean(ari_ap)
