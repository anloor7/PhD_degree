

# This is a function to compute fuzzy c-means for multivariate time series 

# Input parameters:
# Y: list of MTS
# K: number of clusters
# b: fuzziness coefficient
# niter: number of iterations
# tol: threshold for termination step 
# dis : function to compute distance between two MTS objects 

fcm_mts <- function(Y, K, b = 2, niter = 1000, tol = 0.01){
  
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
  mem_matrix <- matrix(nrow = K, ncol = M) # Initialization of membership matrix 
  J <- numeric(niter)  # Initialization of objective function 
  J[1] <- Inf
  
  for (l in 1 : niter) {
    
    
    for (k in 1 : K) {
      for (i in 1 : M) {
        dist_matrix[k, i] <- sw_distance(Z[[k]], Y[[i]]) # Recomputing distance matrix
      }
    }
    
    for (k in 1 : K) {
      for (i in 1 : M) {
        sum_k <- sum(dist_matrix[,i]^(-2/(b-1)))
        mem_matrix[k, i] <- dist_matrix[k, i]^(-2/(b-1))/(sum_k)  # Recomputing membership matrix
        mem_matrix[mem_matrix == 'NaN'] <- 1
      }
    }
    
    
    for (k in 1 : K){
      Z[[k]] <- matrix(0L, nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]]))
      for (i in 1 : M){
        Z[[k]] <- (Z[[k]] + mem_matrix[k, i]^b * Y[[i]]) # Updating the fuzzy cluster centers
      }
      Z[[k]] = Z[[k]]/sum(mem_matrix[k,]^b)
    }
    
    
    
   J[l + 1] = sum(mem_matrix^b*dist_matrix^2)
    
    if ((J[l] - J[l + 1]) < tol){
      
      break
    }
    
  }
  
  
# c <- vector(mode = 'list', length = 2)
# c[[1]] <- mem_matrix
# c[[2]] <- l
# c
  mem_matrix

  
}
