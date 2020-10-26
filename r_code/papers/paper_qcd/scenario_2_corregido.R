
K <- 4
n <- 5
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n))
l <- 250

h1<- matrix(0, nrow = l + 1, ncol = 2)
r1 <- matrix(0, nrow = l + 1, ncol = 2)
h1[1,] <- c(0, 0)

cluster1 <- list()

for (k in (1:n)){

for (i in (1:l)) {
  
  sigma1 <- diag(2)
  sigma1[1, 2] <- 0.5
  sigma1[2, 1] <- sigma1[1, 2]
  epsilon <- mvrnorm(n = 1, mu = c(0, 0), Sigma = sigma1)
  r1[i,1] <- epsilon[1] * sqrt(h1[i, 1])
  r1[i,2] <- epsilon[2] * sqrt(h1[i, 2])
  h1[(i+1),1] <- 0.01 + 0.05*r1[i,1]^2 + 0.94*h1[i, 1]
  h1[(i+1),2] <- 0.5 + 0.2*r1[i,2]^2 + 0.5*h1[i, 2]
  cluster1[[k]] <- r1[(2 : (l+1)),]
  
}
  
  
  
}






h2<- matrix(0, nrow = l + 1, ncol = 2)
r2 <- matrix(0, nrow = l + 1, ncol = 2)
h2[1,] <- c(0, 0)

cluster2 <- list()

for (k in (1:n)){
  
  for (i in (1:l)) {
    
    sigma2 <- diag(2)
    sigma2[1, 2] <- if (i <= (l/2)) {0.7} else {-0.9}
    sigma2[2, 1] <- sigma2[1, 2]
    epsilon <- mvrnorm(n = 1, mu = c(0, 0), Sigma = sigma2)
    r2[i,1] <- epsilon[1] * sqrt(h2[i, 1])
    r2[i,2] <- epsilon[2] * sqrt(h2[i, 2])
    h2[(i+1),1] <- 0.01 + 0.05*r2[i,1]^2 + 0.94*h2[i, 1]
    h2[(i+1),2] <- 0.5 + 0.2*r2[i,2]^2 + 0.5*h2[i, 2]
    
  }
  
  cluster2[[k]] <- r2[(2 : (l+1)),]
  
}






h3<- matrix(0, nrow = l + 1, ncol = 2)
r3 <- matrix(0, nrow = l + 1, ncol = 2)
h3[1,] <- c(0, 0)

cluster3 <- list()

for (k in (1:n)){
  
  for (i in (1:l)) {
    
    sigma3 <- diag(2)
    sigma3[1, 2] <- if (i <= (l/2)) {0.9} else {-0.7}  
    sigma3[2, 1] <- sigma3[1, 2]
    epsilon <- mvrnorm(n = 1, mu = c(0, 0), Sigma = sigma3)
    r3[i,1] <- epsilon[1] * sqrt(h3[i, 1])
    r3[i,2] <- epsilon[2] * sqrt(h3[i, 2])
    h3[(i+1),1] <- 0.01 + 0.05*r3[i,1]^2 + 0.94*h3[i, 1]
    h3[(i+1),2] <- 0.5 + 0.2*r3[i,2]^2 + 0.5*h3[i, 2]
    
  }
  
  cluster3[[k]] <- r3[(2 : (l+1)),]
  
}









h4<- matrix(0, nrow = l + 1, ncol = 2)
r4 <- matrix(0, nrow = l + 1, ncol = 2)
h4[1,] <- c(0, 0)

cluster4 <- list()

for (k in (1:n)){
  
  for (i in (1:l)) {
    
    sigma4 <- diag(2)
    sigma4[1, 2] <- if (mod(i, 2) == 1) {0.99/log(i+2)} else {-0.99/log(i+2)}
    sigma4[2, 1] <- sigma4[1, 2]
    epsilon <- mvrnorm(n = 1, mu = c(0, 0), Sigma = sigma4)
    r4[i,1] <- epsilon[1] * sqrt(h4[i, 1])
    r4[i,2] <- epsilon[2] * sqrt(h4[i, 2])
    h4[(i+1),1] <- 0.01 + 0.05*r4[i,1]^2 + 0.94*h4[i, 1]
    h4[(i+1),2] <- 0.5 + 0.2*r4[i,2]^2 + 0.5*h4[i, 2]
    
  }
  
  cluster4[[k]] <- r4[(2 : (l+1)),]
  
}

cluster <- c(cluster1, cluster2, cluster3, cluster4)

coherence2 <- listTomatrix(lapply(cluster, quantile_quantities_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

features <- lapply(cluster, wavelet_features, wf = wf, J = J) 
dis_matrix <- proxy::dist(features, wave_dist)  
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

features <- listTomatrix(lapply(cluster, gcc_features_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)


dis_matrix <- matrix(0, (n*K), (n*K))
for (p in 1 : (n*K)) {
  for (j in 1 : (n*K)) {
    dis_matrix[p, j] <- j_divergence(cluster[[p]], cluster[[j]])
  }
}

diag(dis_matrix) <- 0 # Numerical error 
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)
