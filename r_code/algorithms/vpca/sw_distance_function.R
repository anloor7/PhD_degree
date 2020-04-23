

# This is a function to compute the spatial wighted distance between two reduced MTS samples

# Input parameters 

# Yi and Yj, two dimensionality reduced MTS, preferably obtained from the output of vpca_function



# Auxiliary function


spatial_distance <- function(i1, j1, i2, j2){
  sqrt((i1 - i2)^2 + (j1 - j2)^2)
}

# Auxiliary function to compute spatial matrix


spatial_matrix <- function(N, ps){
  
  len <- N * ps
  S <- matrix(nrow = len, ncol = len)
  row <- vector(mode = "list", length = len) # Future rows of matrix S
  
  for (i in 1 : len) {
    row[[i]] <- matrix(nrow = N, ncol = ps )
  }

  i1 = 1  
  j1 = 1
  for (i in 1 : len){
  
    
  for (i2 in 1 : N) {
    for (j2 in 1 : ps) {
      row[[i]][i2, j2] <- spatial_distance(i1 = i1, j1 = j1, i2 = i2, j2 = j2)
    
    }
  }
  
    if (i1 < N) {
      i1 = i1 + 1
    } else {
        i1 = i1 - N + 1
        j1 = j1 + 1
        
      }
    
  }
  
  for (i in 1 : len) {
    S[i,] <- as.vector(row[[i]])
  }
  k <- 2*(1-1/ps)^2
  S_spatial  <- 1/(pi * k) * exp(-S^2/k)
  S_spatial
}  






sw_distance <- function(Yi, Yj, spatial_matrix){
  
  yi <- as.vector(Yi)
  yj <- as.vector(Yj)
  sqrt((yi - yj) %*% spatial_matrix %*% (yi - yj))
  
}

