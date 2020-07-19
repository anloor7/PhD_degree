

# This is a function to compute crisp K-means for multivariate time series 

# Input parameters:
# Y: list of MTS
# K: number of clusters
# niter: number of iterations
# tol: threshold for termination step 
# dis : function to compute distance between two MTS objects 


km_mts <- function(Y, K, niter = 1000, tol = 0.01, dis){
  
  M <- length(Y)
  
  if (K > M){
    stop("Number of clusters can not be greater than number of elements")
  }
  
  vector_centroids <- sample(M, K)
  Z <- vector(mode = 'list', length = K)
  
  for (i in 1 : length(vector_centroids)) { # K initial random centroids (MTS)
    Z[[i]] <- Y[[vector_centroids[i]]]
  }
  
  dist_matrix <- matrix(nrow = K, ncol = M) # Initialization of distance matrix
  clus_matrix <- matrix(nrow = K, ncol = M) # Initialization of clustering matrix
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
    
    
    for (k in 1 : K){
      if (is.matrix(Y[[1]])){
        Z[[k]] <- matrix(0L, nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]]))
      } else {
        Z[[k]] <- numeric(length(Y[[1]]))
      }
      for (i in 1 : M){
        Z[[k]] <- Z[[k]] + clus_matrix[k, i] * Y[[i]]  # Updating the cluster centers
      }
      Z[[k]] <- Z[[k]]/sum(clus_matrix[k,])
    }
    
    
    
    J[l + 1] = sum(clus_matrix*dist_matrix^2)
    
    if ((J[l] - J[l + 1]) < tol){
      
      break
    }
    
  }
  
  
  # c <- vector(mode = 'list', length = 2)
  # c[[1]] <- mem_matrix
  # c[[2]] <- l
  # c
  
  clus_matrix
  fuzzytocrisp(clus_matrix)
}
