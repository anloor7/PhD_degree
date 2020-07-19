



# Auxiliary function to retain the minimum number os principal components such that, in each dataset, at least 95% of
# variability gets explained 

# Input parameters:
# Y: list of MTS
# lambda: desired percentage of explained variability

get_r <- function(Y, lambda = 0.95){

N <- length(Y)
pc <- list() # Objects for principal components 
exp_var <- numeric() # Number of components according to explained variability

for (i in 1 : N) {                        
  pc[[i]] <- prcomp(Y[[i]])
  s <- summary(pc[[i]])$importance
  a <- numeric()
  for (j in 1 : ncol(s)) {
    
    if (s[3, j] > lambda) {
      a <- c(a, j)
    }
  }
  exp_var[i] <- a[1]
}
r1 <- max(exp_var)
r1
}


# Auxiliary function to compute the aggregated dataset

# Input parameters:
# Y: list of MTS

laggregate <- function(Y){
  
  M <- length(Y)
  tY <- list()
  
  for (i in 1 : M) {
    
    tY[[i]] <- t(Y[[i]])
    
  }
  
  t(do.call(cbind, tY))
  
}


# Auxiliary function to compute distance between and aggregated dataset and a MTS
# Input parameters:
# A: aggregated dataset
# X MTS
# alpha: weight for SPCA

ag_mts_distance <- function(A, X, alpha = 1){
  
  
  # Similarity PCA (1)
  
  r1 <- ncol(X)
  covA <- cov(A)
  covX <- cov(X)
  sim1 <- PCAsimilarity(covX, covA, ret.dim = r1)
  
  
  # Similarity DIST (2)
  
  SVD <- svd(covA)
  eigenv <- (SVD$d^2)
  eigenvprop <- cumsum(eigenv)/sum(eigenv)
  r2 <- min(which(eigenvprop > 0.95))
  
  if (r2 == 1){
    Dr <- 1/eigenv[1:r2]
  } else {
    Dr <- diag(1/eigenv[1:r2])}
  ur <- SVD$u
  vr <- SVD$v
  Ur <- as.matrix(ur[,1:r2])
  Vr <- as.matrix(vr[,1:r2])
  gi <- Vr %*% Dr %*% t(Ur)
  
  phi <- sqrt(matrix((colMeans(X)-colMeans(A)), 
                     ncol = ncol(A)) %*% gi %*% t(matrix((colMeans(X)-colMeans(A)), ncol = ncol(A))))
  sim2 <- 2 * (1 - pnorm(phi))
  1 - (alpha * sim1 + (1-alpha) * sim2)
  
}


# This is a function to compute the K-means algorithm of Singhal and Seborg for multivariate time series 

# Input parameters:
# Y: list of MTS
# K: number of clusters
# niter: number of iterations
# tol: threshold for termination step 
# dis : function to compute distance between two MTS objects 


km_mts_ss <- function(Y, K, niter = 1000, tol = 0.01, dis = ag_mts_distance){
  
  M <- length(Y)
  
  
  if (K > M){
    stop("Number of clusters can not be greater than number of elements")
  }
  
  # Initializing aggregate dataset 
  

  
  indexes <- sample(K, M, replace = T)
  
  while (length(unique(indexes)) != K) {
    indexes <- sample(K, M, replace = T)
  }
  
  
  Z <-  vector(mode = 'list', length = K)
  
  for (i in 1 : K) {
    positions <- which(indexes == i)
    Z[[i]] <- laggregate(Y[positions])
  }
 

  clus_matrix <- matrix(nrow = K, ncol = M) # Initialization of clustering matrix
  dist_matrix <- matrix(nrow = K, ncol = M) # Initialization of distance matrix
  J <- numeric(niter)  # Initialization of objective function 
  J[1] <- Inf
  
  
  
  
  for (l in 1 : niter) {
    
    
    for (k in 1 : K) {
      for (i in 1 : M) {
        dist_matrix[k, i] <- dis(Z[[k]], Y[[i]]) # Recomputing distance matrix
      }
    }
    
    
    
    for (i in 1 : K) {
      for (j in 1 : M) {
        clus_matrix[i, j] <- dist_matrix[i, j]/min(dist_matrix[,j]) # Recomputing clustering matrix
      }
      
    }
    
    
    clus_matrix[clus_matrix != 1] <- 0
    clus_matrix[clus_matrix == 'NaN'] <- 1
    
    # Updating the cluster centers 
    
    
    indexes <- numeric(M)
    
    for (i in 1 : M) {
      
      indexes[i] <- which(clus_matrix[,i] == 1)
    }
      
    for (i in 1 : K) {
      
      if (length(which(indexes == i)) == 0) {
      } else {
      positions <- which(indexes == i)
      Z[[i]] <- laggregate(Y[positions])
      }
    }
    
    
    sum(clus_matrix*dist_matrix^2)
    J[l + 1] = sum(clus_matrix*dist_matrix^2)
    
    if ((J[l] - J[l + 1]) < tol){
      
      break
    }
    
   
  }
  
  
  
  
  
  indexes
  
  
  
  
  
}  

