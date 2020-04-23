

# This is a function to perfom dimensionality reduction through VPCA  

# Input parameters

# P: list containing MTS (variables in rows and time points in columns)
# lambda: desired percentage of overall explained variability (default, 0.95)
# Warning: the numbers of MTS must be greater than number of time points

vpca <- function(P, lambda = 0.95){
  
  if (length(P) < nrow(P[[1]])){
    
    print('Number of elements must be greater than number of time points!')
  }
  
  N <- nrow(P[[1]]) # Number of variables
  M <- length(P)  # Number of MTS 
  V <- vector(mode = 'list', length = N) # Variable marices 
  pc <- vector(mode = 'list', length = N) # Objects for principal components 
  exp_var <- numeric(N) # Number of components according to explained variability
  Fscores <- vector(mode = 'list', length = N) # Matrixes of dimensionality reduction
  Y <- vector(mode = 'list', length = M) # MTS reduced sample 
  
  for (i in 1 : M) { 
    for (j in 1 : N) {
      V[[j]] <- rbind(V[[j]], P[[i]][j,]) # Filling variable matrices
      
    }
  }

 
  for (i in 1 : N) {                        
    pc[[i]] <- prcomp(V[[i]])
    s <- summary(pc[[i]])$importance
    a <- numeric()
    for (j in 1 : ncol(s)) {
      
      if (s[3, j] > lambda) {
        a <- c(a, j)
      }
    }
    exp_var[i] <- a[1]
  }
  max_exp_var <- max(exp_var)
  
  for (i in 1 : N) {
    Fscores[[i]] <- pc[[i]]$x[, 1 : max_exp_var] # Computing scores matrices
  }
  
  for (i in 1 : M) {
    for (j in 1 : N) {
      Y[[i]] <- rbind(Y[[i]], Fscores[[j]][i,])
    }
    
  }
  
  Y # Output: dimensionality reduced MTS sample 
  
  
}

