

# Creating a function to compute the square modulus of the sequence of quantile coherences in each 
# MTS 

# Input:
# X: a MTS
# levels: vector with the levels of probability associated with the quantiles 

quantile_quantities <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  
  
  coherence <- abs(as.vector(matrix(qSPGv, ncol = 1)))^2
  coherence
  
}


quantile_quantities_re_im <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  
  
  matrix <- matrix(qSPGv, ncol = 1)
  coherence <- c(Re(matrix), Im(matrix))
  
}


quantile_coherence <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getCoherency(qSPG, frequencies = freq)
  
  
  coherence <- abs(as.vector(matrix(qSPGv, ncol = 1)))^2
  coherence
  
}


quantile_coherence_re_im <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getCoherency(qSPG, frequencies = freq)
  
  
  matrix <- matrix(qSPGv, ncol = 1)
  coherence <- c(Re(matrix), Im(matrix))
  
}


