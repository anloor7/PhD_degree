

n <- 10 # Number of series per cluster
l <- 900 # Length 
K <- 4 # Number of clusters
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n))
B <- 20 # Monte Carlo replicas

epsilon <- list()
epsilon[[1]] <- c(0, 0)
epsilon[[2]] <- c(0, 0)
a1 <- 0.1
b2 <- 0.1
sigma <- list()
garch1 <- numeric()
garch2 <- numeric()


set.seed(1234)

cluster <- list() 

for (k in 1 : B){
  
  
  # Cluster 1 
  
  cluster1 <- list()
  
  for (j in 1 : n){
    
    
    for (i in 3 : (l+2)) {
      
      b1 <- 0.3 * exp(epsilon[[i-1]][1]/(1 + exp(epsilon[[i-1]][1])))
      a2 <- 1.5 * exp(epsilon[[i-1]][2]/(1 + exp(epsilon[[i-1]][2])))
      garch1[i] <- exp(a1 + b1 * epsilon[[i-1]][1])
      garch2[i] <- exp(a2 + b2 * epsilon[[i-1]][2])
      
      sigma <- diag(2)
      var_epsilon1 <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
      # var_epsilon2 <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
      sigma[1, 2] <- 0.9
      sigma[2, 1] <- 0.9
      # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma) # Normal
      # epsilon[[i]] <- rmvt(1, df = 3,  sigma = sigma) # Student
       epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma, del = c(3, 3)) # Skewed normal 
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    
    series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[3 : (l + 2)]
    series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[3 : (l + 2)]
    
    cluster1[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  
  
  # Cluster 2
  
  cluster2 <- list()
  
  for (j in 1 : n){
    
    
    for (i in 3 : (l+2)) {
      
      b1 <- 0.3 * exp(epsilon[[i-1]][1]/(1 + exp(epsilon[[i-1]][1])))
      a2 <- 1.5 * exp(epsilon[[i-1]][2]/(1 + exp(epsilon[[i-1]][2])))
      garch1[i] <- exp(a1 + b1 * epsilon[[i-1]][1])
      garch2[i] <- exp(a2 + b2 * epsilon[[i-1]][2])
      
      sigma <- diag(2)
      var_epsilon1 <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
      # var_epsilon2 <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
      sigma[1, 2] <- if(epsilon[[i-1]][1]^2 <=  var_epsilon1) {0.3} else {0.8}
      sigma[2, 1] <- sigma[1, 2] 
      # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma) # Normal
      # epsilon[[i]] <- rmvt(1, df = 3,  sigma = sigma) # Student
       epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma, del = c(3, 3)) # Skewed normal 
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    
    series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[3 : (l + 2)]
    series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[3 : (l + 2)]
    
    cluster2[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  # Cluster 3
  
  cluster3 <- list()
  
  for (j in 1 : n){
    
    
    for (i in 3 : (l+2)) {
      
      b1 <- 0.3 * exp(epsilon[[i-1]][1]/(1 + exp(epsilon[[i-1]][1])))
      a2 <- 1.5 * exp(epsilon[[i-1]][2]/(1 + exp(epsilon[[i-1]][2])))
      garch1[i] <- exp(a1 + b1 * epsilon[[i-1]][1])
      garch2[i] <- exp(a2 + b2 * epsilon[[i-1]][2])
      
      sigma <- diag(2)
      var_epsilon1 <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
      # var_epsilon2 <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
      sigma[1, 2] <- 0.99 - 1.98/(1 + exp(0.5 * max(epsilon[[i-1]][1]^2, epsilon[[i-1]][2]^2)))
      sigma[2, 1] <- sigma[1, 2] 
      # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma) # Normal
      # epsilon[[i]] <- rmvt(1, df = 3,  sigma = sigma) # Student
      epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma, del = c(3, 3)) # Skewed normal 
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    
    series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[3 : (l + 2)]
    series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[3 : (l + 2)]
    
    cluster3[[j]] <- cbind(series1, series2)
    
  }
  
  
  
  # Cluster 4
  
  cluster4 <- list()
  
  for (j in 1 : n){
    
    
    for (i in 3 : (l+2)) {
      
      b1 <- 0.3 * exp(epsilon[[i-1]][1]/(1 + exp(epsilon[[i-1]][1])))
      a2 <- 1.5 * exp(epsilon[[i-1]][2]/(1 + exp(epsilon[[i-1]][2])))
      garch1[i] <- exp(a1 + b1 * epsilon[[i-1]][1])
      garch2[i] <- exp(a2 + b2 * epsilon[[i-1]][2])
      
      sigma <- diag(2)
      var_epsilon1 <- var(listTomatrix(epsilon[1 : (i-1)])[,1])
      # var_epsilon2 <- var(listTomatrix(epsilon[1 : (i-1)])[,2])
      sigma[1, 2] <- 0.5 + 0.4*cos(pi*(i-2)/100)
      sigma[2, 1] <- sigma[1, 2] 
      # epsilon[[i]] <- mvrnorm(1, mu = c(0, 0),  Sigma = sigma) # Normal
      # epsilon[[i]] <- rmvt(1, df = 3,  sigma = sigma) # Student
       epsilon[[i]] <- rdmsn(1, 2, mean = c(0, 0), cov = sigma, del = c(3, 3)) # Skewed normal 
      
    }
    
    epsilonl <- listTomatrix(epsilon)
    
    
    # Creating the series 
    
    
    series1 <- epsilonl[,1][3 : (l + 2)] * sqrt(garch1)[3 : (l + 2)]
    series2 <- epsilonl[,2][3 : (l + 2)] * sqrt(garch2)[3 : (l + 2)]
    
    cluster4[[j]] <- cbind(series1, series2)
    
  }
  
  cluster[[k]] <- c(cluster1, cluster2, cluster3, cluster4)
  
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

