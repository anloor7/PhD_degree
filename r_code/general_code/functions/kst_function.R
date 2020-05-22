

# This is a function to compute J divergence and symmetric Chernoff information divergence between two MTS

# Input parameters:
# X,Y: two MTS (columns are variables and rows are time points) 

j_divergence <- function(X, Y){
  
  n <- nrow(X)
  specx <- spectral.density(X)$operators
  specy <- spectral.density(Y)$operators
  d <- dim(specx)[3]
  
  coeff <- (1/2)*n^(-1)
  vector <- numeric(d)
  
  for (i in 1 : d) {
    
    vector[i] <- tr(specx[,,i] %*% solve(specy[,,i])) + tr(specy[,,i] %*% solve(specx[,,i])) # - 2 * 
  
  }
  
  Re(coeff * sum(vector))
  
}


chernoff_divergence <- function(X, Y, alpha = 0.5){
  
  n <- nrow(X)
  specx <- spectral.density(X)$operators
  specy <- spectral.density(Y)$operators
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
