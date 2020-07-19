

# This is a function to compute J divergence and symmetric Chernoff information divergence between two MTS

# Input parameters:
# X,Y: two MTS (columns are variables and rows are time points) 

j_divergence <- function(X, Y){
  
  n <- nrow(X)
  specx <- spectral.density(X, freq = ((2 * pi * (0 : (nrow(X)/2)))/nrow(X)))$operators
  specy <- spectral.density(Y, freq = ((2 * pi * (0 : (nrow(Y)/2)))/nrow(Y)))$operators
  d <- dim(specx)[3]
  
  coeff <- (1/2)*n^(-1)
  vector <- numeric(d)
  
  for (i in 1 : d) {
    
    vector[i] <- tr(specx[,,i] %*% solve(specy[,,i])) + tr(specy[,,i] %*% solve(specx[,,i]))  - 2 * ncol(X)
  
  }
  
  Re(coeff * sum(vector))
  
}


chernoff_divergence <- function(X, Y, alpha = 0.5){
  
  n <- nrow(X)
  specx <- spectral.density(X, freq = ((2 * pi * (0 : (nrow(X)/2)))/nrow(X)))$operators
  specy <- spectral.density(Y, freq = ((2 * pi * (0 : (nrow(Y)/2)))/nrow(Y)))$operators
  d <- dim(specx)[3]
  
  coeff <- (1/2)*n^(-1)
  vector <- numeric(d)
  
  for (i in 1 : d) {
    
    a <-  log(Det(alpha*specx[,,i] + (1-alpha) * specy[,,i])/Det(specy[,,i]))
    b <-  log(Det(alpha*specy[,,i] + (1-alpha) * specx[,,i])/Det(specx[,,i]))
    vector[i] <- a + b
    
  }
  
  Re(coeff * sum(vector))
  
}


# Now we are going to develop their counterparts. The aim now is to compute j_divergence and chenorff divergence between 
# a list of MTS and a MTS. As for the list of MTSs, we are going to take into account the average spectral matrix

# Input parameters:
# L: list of MTSs
# Y: MTS

j_divergence_average <- function(L, X){
  
  l <- length(L)
  n <- nrow(X)
  specx <- spectral.density(X, freq =  ((2 * pi * (0 : (nrow(X)/2)))/nrow(X)))$operators
  specl <- array(dim = c(dim(specx), l))
  speclav <- array(0, dim = dim(specl)[1 : 3])
  
  for (i in 1 : l) {
      specl[,,,i] <- spectral.density(L[[i]], freq = ((2 * pi * (0 : (nrow(L[[i]])/2)))/nrow(L[[i]])))$operators
  }
  
  for (i in 1 : l) {
    speclav <- speclav + specl[,,,i]
  }
  
  speclav <- speclav/l



  
  d <- dim(specx)[3]
  
  coeff <- (1/2)*n^(-1)
  vector <- numeric(d)
  
  for (i in 1 : d) {
    
    vector[i] <- tr(specx[,,i] %*% solve(speclav[,,i])) + tr(speclav[,,i] %*% solve(specx[,,i]))  - 2 * ncol(X) 
    
  }
  
  Re(coeff * sum(vector))
  
}




chernoff_divergence_average <- function(L, X, alpha = 0.5){
  
  l <- length(L)
  n <- nrow(X)
  specx <- spectral.density(X, freq = ((2 * pi * (0 : (nrow(X)/2)))/nrow(X)))$operators
  specl <- array(dim = c(dim(specx), l))
  speclav <- array(0, dim = dim(specl)[1 : 3])
  
  
  for (i in 1 : l) {
    specl[,,,i] <- spectral.density(L[[i]], freq = ((2 * pi * (0 : (nrow(L[[i]])/2)))/nrow(L[[i]])))$operators
  }
  
  for (i in 1 : l) {
    speclav <- speclav + specl[,,,i]
  }
  
  speclav <- speclav/l
  
  d <- dim(specx)[3]
  
  coeff <- (1/2)*n^(-1)
  vector <- numeric(d)
  
  for (i in 1 : d) {
    
    a <-  log(Det(alpha*specx[,,i] + (1-alpha) * speclav[,,i])/Det(speclav[,,i]))
    b <-  log(Det(alpha*speclav[,,i] + (1-alpha) * specx[,,i])/Det(specx[,,i]))
    vector[i] <- a + b
    
  }
  
  Re(coeff * sum(vector))
  
}




