

K <- 3
l <- 300
n <- 10
B <- 20 # Monte Carlo replicas

ground_truth <- c(rep(1, n), rep(2, n), rep(3, n))

set.seed(1234)

cluster <- list() 

for (k in 1 : B){

# Cluster 1

cluster1 <- list()

for (j in 1 : n){

Xt <- list()
Xt[[1]] <- c(0, 0)

for (i in 2 : (l + 1)) {
  
  theta <- 0*diag(2)
  u12 <- runif(1)
  u21 <- runif(1)
  theta[1, 2] <- 1.2*(u12 - 0.5)
  theta[2, 1] <- 1.2*(u21 - 0.5)
  
  theta0 <- numeric()
  theta0[1] <- qnorm(u12, 0)
  theta0[2] <- qnorm(u21, 0)
  
  Xt[[i]] <- theta %*% Xt[[i - 1]] + theta0
  
}

cluster1[[j]] <- listTomatrix(Xt)


}



# Cluster 2 

cluster2 <- list()

for (j in 1 : n){
  
  Xt <- list()
  Xt[[1]] <- c(0, 0)
  Xt[[2]] <- c(0, 0)
  
  for (i in 3 : (l + 2)) {
    
    theta1 <- 0*diag(2)
    u12 <- runif(1)
    u21 <- runif(1)
    
    theta2 <- 0*diag(2)
    theta2[1, 2] <- 1.2*(u12 - 0.5)
    theta2[2, 1] <- 1.2*(u21 - 0.5)
    
    theta0 <- numeric()
    theta0[1] <- qnorm(u12, 0)
    theta0[2] <- qnorm(u21, 0)
    
    Xt[[i]] <- theta1 %*% Xt[[i - 1]] + theta2 %*% Xt[[i - 2]] + theta0
    
  }
  
  cluster2[[j]] <- listTomatrix(Xt)[-(l+2),]
  
  
}


  

# Cluster 3


cluster3 <- list()

for (j in 1 : n){
  
  Xt <- list()
  Xt[[1]] <- c(0, 0)
  Xt[[2]] <- c(0, 0)
  Xt[[3]] <- c(0, 0)
  
  for (i in 4 : (l + 3)) {
    
    theta1 <- 0*diag(2)
    u12 <- runif(1)
    u21 <- runif(1)
    
    theta2 <- 0 * diag(2)
    
    theta3 <- 0*diag(2)
    theta3[1, 2] <- 1.2*(u12 - 0.5)
    theta3[2, 1] <- 1.2*(u21 - 0.5)
    
    theta0 <- numeric()
    theta0[1] <- qnorm(u12, 0)
    theta0[2] <- qnorm(u21, 0)
    
    Xt[[i]] <- theta1 %*% Xt[[i - 1]] + theta2 %*% Xt[[i - 2]] + theta3 %*% Xt[[i - 3]] + theta0
    
  }
  
  cluster3[[j]] <- listTomatrix(Xt)[-c((l+2), (l+3)),]
  
  
  
}


cluster[[k]] <- c(cluster1, cluster2, cluster3)

}

ari_qc2 <- numeric()

for (i in 1 : B) {
  
  coherence2 <- listTomatrix(lapply(cluster[[i]], quantile_coherence_re_im))
  dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
  clustering <- pam(dis_matrix, K)$cluster
  ari_qc2[i] <- external_validation(ground_truth, clustering)
  
}

mean(ari_qc2)

# Wavelets 

J <- 4 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
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
