
# All the MTS subjected to cluster are provided in a list named MTS
# MTS <- list()

###############################################################
# Auxiliary functions
###############################################################

# Computing weights for EROS distance according to 
# Algorithms 1 and 2 in Yang and Shahabi (2005)
###############################################################
# INPUTS:
# sigma: K*M matrix. The i-th column contains the PCA eigenvalues for i-th MTS
# method1: "mean" (default), "min", or "max"
# method2: "normalized" (normalized eigenvalues, default) or "raw" (raw eigenvalues)  
# OUTPUT: vector of weights 

w.eros <- function(sigma, method1 = "mean", method2 = "normalized")
{
  if (method2 == "normalized") sigma <- sigma %*% diag(1/colSums(sigma)) 
  statistic <- get(method1)
  pesos <- apply(sigma,1,statistic)
  w <- pesos/sum(pesos)
  return(w)
}

# Computing the EROS distance Yang and Shahabi (2005)
###############################################################
# INPUTS:
# A, B: PCA eigenvector matrices for two MTS
# w: weights from w.eros() function 
# OUTPUT: Eros distance 
Eros <- function(A,B,w) sqrt( 2- 2*sum(w*abs(diag(t(A)%*%B))) )

###############################################################
# MAIN BODY
###############################################################
# PCAs  <- lapply(MTS,prcomp)                      # PCA  
# SIGMA <- (do.call(cbind,lapply(PCAs,'[[',1)))^2  # matrix of eigenvalues
# V <- lapply(PCAs,'[[',2)                         # list with eigenvector matrices
# w <- w.eros(SIGMA)                               # weights
# d <- proxy::dist(V,V,Eros,w)                     # distances
# d <- as.dist(d)                                  # as a dist object for clustering

# Note: Do not worry about NAs. We have to find zeros in the diagonal, but due to the 
#       approximations we can potentially find values very close to zero. If these 
#       approximations are negative, sqrt() generates NAs. So, no problem. 


