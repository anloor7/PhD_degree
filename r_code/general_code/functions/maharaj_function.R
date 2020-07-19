

# Creating a function to compute VARMA coefficients of a list of series via AIC criterion (maximum order
# VAR is taken as a reference)

# Input:
# X: a list MTS


maharaj <- function(X){
  
  l <- length(X)
  coefficients <- list()
  
  for (i in (1 : l)) {
    
    auxiliary <- numeric(4)
    p <- seq(0, 4)
    
    for (j in p) {
        
        auxiliary[j] <- VARMA(X[[i]], p = j)$aic
        
      }
    
    
    min_pos <- which(auxiliary == min(auxiliary), arr.ind = T)
    coefficients[[i]] <- as.vector(VARMA(X[[i]], p = min_pos[1])$coef)
    
  }
  
  lengths = numeric()
  lc <- length(coefficients)
  
  for (i in (1 : lc)) {
    
    lengths[i] <- length(coefficients[[i]])
    
  }
  
  max_length <- max(lengths)
  
  for (i in (1 : l)) {
    
    if (length(coefficients[[i]]) != max_length) {
      
      length(coefficients[[i]]) <- li
      add <- numeric(li)
      coefficients[[i]] <- c(coefficients[[i]], add)
      
    }
    
  }
  
  coefficients
  
  
}
