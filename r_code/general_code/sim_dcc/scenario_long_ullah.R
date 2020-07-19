

n <- 5 # Number of series per cluster
l <- 1000 # Length 
K <- 3 # Number of clusters
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n))



# Cluster 1 

cluster1 <- list()
sigma1 <- list()



# Simulating GARCH processes 

for (j in 1 : n){
  
  garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
  garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
  garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
  garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
  
  
  # Simulating standardized innovations 
  
  # Time variant covariance matrixes 
  
  
  epsilon <- list()
  epsilon[[1]] <- c(0, 0)
  epsilon[[2]] <- c(0, 0)
  var_epsilon1[[1]] <- 0
  var_epsilon1[[2]] <- 0
  var_epsilon2[[1]] <- 0
  var_epsilon2[[2]] <- 0
  var_epsilon1 <- list()
  var_epsilon2 <- list()
  
  
  for (i in 3 : (l+2)) {
    
    sigma1[[i]] <- diag(2)
    var_epsilon1[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
    var_epsilon2[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
    sigma1[[i]][1, 2] <- if(epsilon[[i-1]][1]^2 <=  var_epsilon1[[i]]) {0.3} else {0.8}
    sigma1[[i]][2, 1] <- sigma1[[i]][1, 2] 
    # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma1[[i]]) # Normal
    epsilon[[i]] <- rmvt(1,  sigma = sigma1[[i]]) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma1[[i]], del = c(4, 4)) # Skewed normal 
    
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  

  series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[2 : (l + 1)]
  
  cluster1[[j]] <- cbind(series1, series2)
  
}






# Cluster 2

cluster2 <- list()
sigma2 <- list()



# Simulating GARCH processes 

for (j in 1 : n){
  
  garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
  garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
  garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
  garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
  
  
  # Simulating standardized innovations 
  
  # Time variant covariance matrixes 
  
  
  epsilon <- list()
  epsilon[[1]] <- c(0, 0)
  epsilon[[2]] <- c(0, 0)
  var_epsilon1[[1]] <- 0
  var_epsilon1[[2]] <- 0
  var_epsilon2[[1]] <- 0
  var_epsilon2[[2]] <- 0
  var_epsilon1 <- list()
  var_epsilon2 <- list()
  
  
  for (i in 3 : (l+2)) {
    
    sigma2[[i]] <- diag(2)
    var_epsilon1[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
    var_epsilon2[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
    sigma2[[i]][1, 2] <- 0.99 - 1.98/(1 + exp(0.5 * max(epsilon[[i-1]][1]^2, epsilon[[i-1]][2]^2)))
    sigma2[[i]][2, 1] <- sigma2[[i]][1, 2] 
    # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma2[[i]]) # Normal
    epsilon[[i]] <- rmvt(1,  sigma = sigma2[[i]]) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma2[[i]], del = c(4, 4)) # Skewed normal 
    
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  
  
  series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[2 : (l + 1)]
  
  cluster2[[j]] <- cbind(series1, series2)
  
}



# Cluster 3

cluster3 <- list()
sigma3 <- list()



# Simulating GARCH processes 

for (j in 1 : n){
  
  garch1p <- garch.sim(c(0.01, 0.05), 0.94, n = l + 1)
  garch2p <- garch.sim(c(0.5, 0.2), 0.5, n = l + 1)
  garch1 <- garch(garch1p, order = c(1, 1))$fitted.values[,1]
  garch2 <- garch(garch2p, order = c(1, 1))$fitted.values[,1]
  
  
  # Simulating standardized innovations 
  
  # Time variant covariance matrixes 
  
  
  epsilon <- list()
  epsilon[[1]] <- c(0, 0)
  epsilon[[2]] <- c(0, 0)
  var_epsilon1[[1]] <- 0
  var_epsilon1[[2]] <- 0
  var_epsilon2[[1]] <- 0
  var_epsilon2[[2]] <- 0
  var_epsilon1 <- list()
  var_epsilon2 <- list()
  
  
  for (i in 3 : (l+2)) {
    
    sigma3[[i]] <- diag(2)
    var_epsilon1[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
    var_epsilon2[[i]] <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
    sigma3[[i]][1, 2] <- 0.5 + 0.4*cos(pi * i/100)
    sigma3[[i]][2, 1] <- sigma3[[i]][1, 2] 
    # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma3[[i]]) # Normal
    epsilon[[i]] <- rmvt(1,  sigma = sigma3[[i]]) # Student
    # epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma3[[i]], del = c(4, 4)) # Skewed normal 
    
  }
  
  epsilonl <- listTomatrix(epsilon)
  
  
  # Creating the series 
  
  
  series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[2 : (l + 1)]
  series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[2 : (l + 1)]
  
  cluster3[[j]] <- cbind(series1, series2)
  
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
