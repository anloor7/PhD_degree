

# Given a list of bivariate MTS, we fit a VARMA(1) model each one, and get the 8-dimensional vector containing
# the coefficients and AIC and BIC criteria 

library(MTS)

library(TSdist)

varma_coefs <- function(X){
  coefs <- list()
  n <- length(X)
  
  for (i in 1 : n){
  varmax <- VARMA(X[[i]], p = 1)
  coefx <- c(as.vector(varmax$coef), varmax$aic, varmax$bic)
  coefs[[i]] <- coefx
  }
  coefs
}

