

# Given a list of MTS, we fit a VARMA(1, 1) model each one, and get the 8-dimensional vector containing the coefficients

library(MTS)

library(TSdist)

varma_coefs <- function(X){
  coefs <- list()
  n <- length(X)
  
  for (i in 1 : n){
  varmax <- VARMA(X[[i]], p = 1, q = 0)
  coefx <- c(as.vector(varmax$Phi), as.vector(varmax$Theta), as.vector(varmax$Sigma))
  coefs[[i]] <- coefx
  }
  coefs
}
