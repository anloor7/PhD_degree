


# Extraer las features de QAF correspondientes a un par arbitrario de UTS que componen la MTS, X e Y


qaf_uts_coefs_xy_sep <- function(X, Y, levels = c(0.1, 0.5, 0.9), l = 1) {
  
  k <- length(levels)
  nx <- length(X)
  qx <- quantile(X, probs = levels)
  qy <- quantile(Y, probs = levels)
  gamma <- matrix(0, k, k) # Future vector gamma for X, regarding lag 1
  
  for (i in 1 : k) {
    for (j in 1 : k) {
      ax <- (X <= qx[i])[1 : (nx - l)]
      bx <- (Y[(1 + l) : nx] <= qy[j])
      gamma[i, j] = (1/(nx - l))*sum(ax * bx)-levels[i]*levels[j]
    }
  }
  
  as.vector(gamma)
  
}

#  Extraer las features de QAF de una MTS X 


qaf_mts_coefs_xy_sep <- function(X){
  c <- ncol(X)
  gamma <- matrix(ncol = 9) # 9 if l = 1 and there are 3 probability levels 
  for (i in 1 : c) {
    for (j in 1 : c){
      gammanew <- qaf_uts_coefs_xy_sep((X[,i]), (X[,j]))
      gamma <- rbind(gamma, gammanew)
    }
  }
  gamma <- gamma[-1,]
  as.vector(gamma)
  
}


# Extraer las features de QC1 de una MTS X


quantile_coherence <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  
  
  coherence <- abs(as.vector(matrix(qSPGv, ncol = 1)))^2
  coherence
  
}

# Extraer las features de QC2 de una MTS X


quantile_coherence_re_im <- function(X, levels = c(0.1, 0.5, 0.9)) {
  
  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  
  
  matrix <- matrix(qSPGv, ncol = 1)
  coherence <- c(Re(matrix), Im(matrix))
  
}