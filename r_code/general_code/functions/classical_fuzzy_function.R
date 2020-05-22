

# This is a function to perform classical fuzzy C-means clustering with euclidean distance and cross_sectional data. 

# Input parameters:
# X: dataset
# K: number of clusters
# b: cluster fuzziness
# niter: number of iterations
# tol: threshold for termination step 

fuzzy_cm <- function(X, K, b, niter, tol){
  
  n <- nrow(X)
  vector_centroids <- sample(n, K) 
  matrix_centroids <- X[vector_centroids,] # K initial random centroids
  dist_matrix <- matrix(nrow = K, ncol = n) # Initialization of distance matrix
  mem_matrix <- matrix(nrow = K, ncol = n) # Initialization of membership matrix 
  J <- numeric(niter)  # Initialization of objective function 
  J[1] <- 0
  
  for (t in 2 : niter) {
    
    
    for (k in 1 : K) {
      for (i in 1 : n) {
        dist_matrix[k, i] <- sqrt((sum((matrix_centroids[k,] - X[i,])^2))) # Recomputing distance matrix
      }
    }
    
    
    for (k in 1 : K) {
      for (i in 1 : n) {
        sum_k <- sum(dist_matrix[,i]^(-2/(b-1)))
        mem_matrix[k, i] <- dist_matrix[k, i]^(-2/(b-1))/(sum_k)  # Recomputing membership matrix
        mem_matrix[mem_matrix == 'NaN'] <- 1
      }
    }
    
  for (k in 1 : K){
    matrix_centroids[k,] <- numeric(ncol(X))
    for (i in 1 : n){
      matrix_centroids[k,] <- (matrix_centroids[k,] + mem_matrix[k, i]^b * X[i,]) # Updating the fuzzy cluster centers
    }
    matrix_centroids[k,] = matrix_centroids[k,]/sum(mem_matrix[k,]^b)
  }
    
    J[t] = sum(mem_matrix^b*dist_matrix^2)
    
    if (J[t] - J[t-1] < tol){
      break
    }
    
    
    
  }
  mem_matrix
}


