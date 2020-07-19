# All the MTS subjected to cluster are provided in a list named MTS
# MTS <- list()

####################################################################
# Required libraries
####################################################################
library(waveslim)

####################################################################
# Auxiliary functions
####################################################################
# 1. wavelet.features()
# 2. wave.dist()

# 1. wavelet.features(): Replace a MTS by a wavelet feature matrix.
# The output matrix has so many columns as preselected wavelet scales
# Each colum constains the vectorized wavelet variance–correlation matrix 
# for the corresponding scale (see (3.4) in D'urso and Maharaj (2012)).
# INPUTS:
#     Xi: A MTS object
#     wf: Wavelet filter 
#     J : number of scales
# OUTPUT: Wavelet variance–correlation matrix (one column per scale)

wavelet_features <- function(Xi,wf,J)
{
N <- dim(Xi)[1]
K <- dim(Xi)[2]
Xi.modwt    <- apply(Xi,2,modwt,wf,J)
Xi.modwt.bw <- lapply(Xi.modwt,brick.wall,wf)
Xi.var <- lapply(Xi.modwt.bw,wave.variance)
Xi.var.2 <- (do.call(cbind,lapply(Xi.var,'[[',1)))
features <- t(Xi.var.2)
for (j in 1:(K-1))
  for (i in (j+1):K)
  {
    x <- wave.correlation(Xi.modwt.bw[[j]], Xi.modwt.bw[[i]], N)[,1]
    features <- rbind(features,x)
  }
return(features[,-dim(features)[2]])
}


# 2. wave.dist(): Compute the Euclidean distance between a
# pair (A and B) of wavelet variance–correlation matrices
####################################################################

wave_dist <- function(A,B)
{
  L <- ncol(A)
  d.escalas <- numeric(length=L)
  for (i in 1:ncol(A))
    d.escalas[i] <- sqrt(sum((A[,i]-B[,i])^2))
  dist <- sum(d.escalas)
  return(dist)
}


###############################################################
# MAIN BODY
###############################################################
# J <- 6      # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
# wf <- "d4"  # wavelet filter
# features <- lapply(MTS,wavelet.features,wf,J) 
# d <- proxy::dist(features,wave.dist)  





