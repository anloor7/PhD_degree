

# These are three functions to perform Mc2PCA

# The first function computes the common space and the second function computes the clustering solution 

# First function: cpca
# Input:
# sigma: list containing a list of covariance matrices of MTS (taking into account the variables)
# lambda: percentage of explained variability

cpca <- function(sigma, lambda = 0.70){

  n <- length(sigma)
  sigma_t <- matrix(0, nrow = nrow(sigma[[1]]), ncol = ncol(sigma[[1]]))

  for (i in 1 : n) {
    sigma_t <- sigma_t + sigma[[i]]
  }
  svd <- eigen((1/n) * sigma_t)
  u <-svd$vectors
  values <- svd$values
  prop <- cumsum(values/sum(values))
  index <-which((prop > lambda) == TRUE)[1]
  s <- u[, 1 : index]
  s
}


# The second function computes the clustering from a list of length K (numbero f clusters) indicating indexes


listToclustering <- function(l){
  k <- length(l)
  clustering <- numeric(k)
  for (i in 1 : k) {
    clustering[l[[i]]] <- i
  }
  clustering
}

# The second function computes the pre-clustering results. It performs  Mc2PCA algorithm

# Second function: Mc2PCA
# Input:
# X: list containing MTS (columns must be variables ands rows must be time points)
# K: number of clusters
# lambda: percentage of explained variability 
# niter: number of iterations
# tol: tolerance

mc2pca <- function(X, K, lambda, niter, tol){
  
  n <- length(X)
  L <- floor(n/K)
  sigma <- vector(mode = 'list', length = n)
  error_matrix <- matrix(0, nrow = n, ncol = K) # Initialization of error matrix
  
  for (i in 1 : n) {
    for (j in 1 : ncol(X[[i]])) {
      X[[i]][,j] <- X[[i]][,j] - mean((X[[i]])[,j]) 
      sigma[[i]] <- cov(X[[i]]) # Covariance matrices of each normalized MTS
    }
  }
  
  
  
  l <- vector(mode = 'list') 
  s <- vector(mode = 'list') # Initialization of projection axes
  sample <- sample(n)
  for (k in 1 : (K-1)) {
    l[[k]] <- sample[(k - 1) * L + 1 : L]
  }
  l[[K]] <- sample[(K - 1) * L + 1 : (L + (length(sample) - K * L))]
 
  for(k in 1 : K){
   s[[k]] <- cpca(sigma[l[[k]]], lambda = lambda)
  }
  
  
  for (niter in 1 : niter){
  
  for (i in 1 : n) {
    for (k in 1 : K) {
      Y <- X[[i]] %*% s[[k]] %*% t(s[[k]])
      error_matrix[i, k] <- norm(X[[i]] - Y, type = 'F') 
    }
  }
    
    discrete_matrix <- matrix(0, nrow = n, ncol = K)  
    for (i in 1 : n) {
      for (k in 1 : K) {
        discrete_matrix[i, k] <- error_matrix[i, k]/min(error_matrix[i,])
      }
    }  
    
    for (k in 1 : K) {
      l[[k]] <- which(discrete_matrix[,k] == 1)
    }
    
    
    
  }   
  
  listToclustering(l)
}


