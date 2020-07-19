
n <- 10 # Number of series per cluster
l <- 500 # Length 
K <- 4 # Number of clusters
B <- 20 # Monte Carlo replicas
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n))

set.seed(1234)
cluster <- list()

for (k in 1 : B){
  
  # Cluster 1 
  
  cluster1 <- list()
  sigma1 <- list()
  
  for (i in 1 : l) {
    
    sigma1[[i]] <- diag(2)
    sigma1[[i]][1, 2] <- 0.8
    sigma1[[i]][2, 1] <- sigma1[[i]][1, 2]
    
  }
  
  # Simulating GARCH processes 
  
  for (j in 1 : n){
    
    garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
    garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
    garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
    garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
    
    
    # Simulating standardized innovations 
    
    # Time variant covariance matrixes 
    
    
    epsilon <- list()
    
    for (i in 1 : l) {
      
      # epsilon[[i]] <- rmvnorm(1,  sigma = sigma1[[i]]) # Normal
       epsilon[[i]] <- rmvt(1,  sigma = sigma1[[i]], df = 3) # Student
      # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma1[[i]], del = c(4, 4)) # Skewed normal 
      # epsilon[[i]] <- mrstab(1,4,1.3,c(0.1,0.5,0.5,0.1),c(0,0))
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
    series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
    
    cluster1[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  # Cluster 2 
  
  cluster2 <- list()
  sigma2 <- list()
  
  for (i in 1 : l) {
    
    sigma2[[i]] <- diag(2)
    sigma2[[i]][1, 2] <- 0.5 + 0.4 * cos(2 * pi * i/200)
    sigma2[[i]][2, 1] <- sigma2[[i]][1, 2]
    
  }
  
  # Simulating GARCH processes 
  
  for (j in 1 : n){
    
    garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
    garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
    garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
    garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
    
    
    # Simulating standardized innovations 
    
    # Time variant covariance matrixes 
    
    
    epsilon <- list()
    
    for (i in 1 : l) {
      
      # epsilon[[i]] <- rmvnorm(1,  sigma = sigma2[[i]]) # Normal
       epsilon[[i]] <- rmvt(1,  sigma = sigma2[[i]], df = 3) # Student
      # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma2[[i]], del = c(4, 4)) # Skewed normal 
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
    series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
    
    cluster2[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  # Cluster 3 
  
  cluster3 <- list()
  sigma3 <- list()
  
  for (i in 1 : l) {
    
    sigma3[[i]] <- diag(2)
    sigma3[[i]][1, 2] <- if (i < (l/2)) {0.9} else {0.6}
    sigma3[[i]][2, 1] <- sigma3[[i]][1, 2]
    
  }
  
  # Simulating GARCH processes 
  
  for (j in 1 : n){
    
    garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
    garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
    garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
    garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
    
    
    # Simulating standardized innovations 
    
    # Time variant covariance matrixes 
    
    
    epsilon <- list()
    
    for (i in 1 : l) {
      
      # epsilon[[i]] <- rmvnorm(1,  sigma = sigma3[[i]]) # Normal
       epsilon[[i]] <- rmvt(1,  sigma = sigma3[[i]], df = 3) # Student
      # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma3[[i]], del = c(4, 4)) # Skewed normal 
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
    series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
    
    cluster3[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  
  # Cluster 4
  
  cluster4 <- list()
  sigma4 <- list()
  
  for (i in 1 : l) {
    
    sigma4[[i]] <- diag(2)
    sigma4[[i]][1, 2] <- 0.4*mod(i, 3)
    sigma4[[i]][2, 1] <- sigma4[[i]][1, 2]
    
  }
  
  # Simulating GARCH processes 
  
  for (j in 1 : n){
    
    garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
    garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
    garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
    garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
    
    
    # Simulating standardized innovations 
    
    # Time variant covariance matrixes 
    
    
    epsilon <- list()
    
    for (i in 1 : l) {
      
      # epsilon[[i]] <- rmvnorm(1,  sigma = sigma4[[i]]) # Normal
       epsilon[[i]] <- rmvt(1,  sigma = sigma4[[i]], df = 3) # Student
      # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma4[[i]], del = c(4, 4)) # Skewed normal  
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
    series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
    
    cluster4[[j]] <- cbind(series1, series2)
    
  }
  
  
  cluster[[k]] <- c(cluster1, cluster2, cluster3, cluster4)
  
}


for (k in 1 : B) {
  for (j in (1 : (4 * n))) {
    cluster[[k]][[j]] <- cluster[[k]][[j]][100 : l,]
  }
  
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

J <- 5 # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
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


