


n <- 5 # Number of series per cluster
l <- 3000 # Length 
K <- 3 # Number of clusters
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n))



cluster1 <- list()
sigma1 <- list()

for (i in 1 : l) {
  
  sigma1[[i]] <- diag(2)
  sigma1[[i]][1, 2] <- 0.8*(-l/2+i)
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
    
    epsilon[[i]] <- rmvnorm(1,  sigma = sigma1[[i]]) # Normal
    # epsilon[[i]] <- rmvt(1,  sigma = sigma1[[i]], df = 3) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma1[[i]], del = c(4, 4)) # Skewed normal 
    # epsilon[[i]] <- mrstab(1,4,1.3,c(0.1,0.5,0.5,0.1),c(0,0))
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  
  series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
  
  cluster1[[j]] <- cbind(series1, series2)
  
}




cluster2 <- list()
sigma2 <- list()

for (i in 1 : l) {
  
  sigma2[[i]] <- diag(2)
  sigma2[[i]][1, 2] <-  0.8*sin(i)
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
    
    epsilon[[i]] <- rmvnorm(1,  sigma = sigma2[[i]]) # Normal
    # epsilon[[i]] <- rmvt(1,  sigma = sigma2[[i]], df = 3) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma2[[i]], del = c(4, 4)) # Skewed normal 
    
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  
  series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
  
  cluster2[[j]] <- cbind(series1, series2)
  
}




cluster3 <- list()
sigma3 <- list()

for (i in 1 : l) {
  
  sigma3[[i]] <- diag(2)
  sigma3[[i]][1, 2] <- if (mod(i, 2) == 1) {0.99/log(i+5)} else {-0.99/log(i+5)}
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
    
    epsilon[[i]] <- rmvnorm(1,  sigma = sigma3[[i]]) # Normal
    # epsilon[[i]] <- rmvt(1,  sigma = sigma2[[i]], df = 3) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma2[[i]], del = c(4, 4)) # Skewed normal 
    
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  
  series1 <- epsilonl[,1] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2] * sqrt(garch2)[2 : (l + 1)]
  
  cluster3[[j]] <- cbind(series1, series2)
  
}


cluster <- c(cluster1, cluster2, cluster3)

# QC2

coherence2 <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# DM

features <- lapply(cluster, wavelet_features, wf = wf, J = J) 
dis_matrix <- proxy::dist(features, wave_dist)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# AP

features <- listTomatrix(lapply(cluster, gcc_features_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

# KST

dis_matrix <- matrix(0, 15, 15)
for (i in 1 : 15) {
  for (j in 1 : 15) {
    dis_matrix[i, j] <- j_divergence(cluster[[i]], cluster[[j]])
  }
}

diag(dis_matrix) <- 0 # Numerical error 
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

