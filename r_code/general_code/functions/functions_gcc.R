

# Lets create a function to compute the generalized cross correlation distance (Alonso and Peña 2019)


# Generalize cross correlation based distance function between univariate time series
# Input parameters:
# X, Y: two univariate time series
# k: desired lag 

gcc_uts <- function(X, Y, k = 1){
  
  X <- as.vector(scale(X))
  Y <- as.vector(scale(Y))
  
  
  # Computing Rxx
  
  
  acfx <- acf(X, lag.max = k)$acf[,,1]
  Rxx <- matrix(0, nrow = k + 1, ncol = k + 1)
  
  for (i in 2 : (k + 1)) {
    Rxx[(i-1), i : (k + 1)] <- acfx[1 : (k + 2 - i)]
  }
  
  for (i in 1 : (k + 1)) {
    for (j in 1 : (i - 1)) {
     Rxx[i, j] <- Rxx[j, i]
    }
  }
  
  diag(Rxx) <- 1
  
  
  # Computing Ryy
  
  
  acfy <- acf(Y, lag.max = k)$acf[,,1]
  Ryy <- matrix(0, nrow = k + 1, ncol = k + 1)
  
  for (i in 2 : (k + 1)) {
    Ryy[(i-1), i : (k + 1)] <- acfy[1 : (k + 2 - i)]
  }
  
  for (i in 1 : (k + 1)) {
    for (j in 1 : (i - 1)) {
      Ryy[i, j] <- Ryy[j, i]
    }
  }
  
  diag(Ryy) <- 1
  
  
  # Computing Cxy 
  
  ccfxy <- ccf(X, Y, lag.max = k)$acf[,,1]
  Cxy <- matrix(0, nrow = k + 1, ncol = k + 1)
  
  for (i in 1 : (k + 1)) {
    
    Cxy[i,] <- rev(ccfxy[i : (k + i)])
  }
  
  
  # Computing the generalize cross correlation based distance function
  
  detn <- det(Ryy - Cxy %*% solve(Rxx)  %*% t(Cxy))
  detd <- det(Ryy)
  distance <- (detn^(1/(k + 1)))/(detd^(1/(k + 1)))
  
  if (is.na(distance)){
    0
  } else {
    distance
  }
  
}

# Generalized cross correlation based distance function between multivariate time series. Extracting features
# Input parameters:
# X: list of mts 
# k: desired lag 


gcc_features_mts <- function(X, k = 1){
  
  c <- ncol(X)
  dis <- 0
  features <- matrix(0, c, c)
  
  for (i in 1 : c) {
    for (j in 1 : c){
    features[i, j] <- gcc_uts(X[,i], X[,j], k = k)
    }
  }
  
  as.vector(features)
  
}







