

# This is a function to compute the QAC distance for a MTS (as the sum of univariate QAC distances)


# First, we construct a function to compute QAC distance between two univariate time series

# Input:
# X, Y, two univariate time series

qaf_uts <- function(X, Y, levels = c(0.1, 0.5, 0.9), l = 1) {
  
  k <- length(levels)
  nx <- length(X) 
  ny <- length(Y)
  qx <- quantile(X, probs = levels)
  qy <- quantile(Y, probs = levels)
  gammax <- matrix(0, k, k) # Vector gamma for X, regarding lag 1
  gammay <- matrix(0, k, k) # Vector gamma for Y, regarding lag 1
  
  for (i in 1 : k) {
    for (j in 1 : k) {
      ax <- (X <= qx[i])[1 : (nx - l)]
      bx <- (X[(1 + l) : nx] <= qx[j])
      ay <- (Y <= qy[i])[1 : (ny - l)]
      by <- (Y[(1 + l) : ny] <= qy[j])
      gammax[i, j] = (1/(nx - l))*sum(ax * bx)-levels[i]*levels[j]
      gammay[i, j] = (1/(ny - l))*sum(ay * by)-levels[i]*levels[j]
    }
  }
  
  EuclideanDistance(as.vector(gammax), as.vector(gammay))^2
  
}

# Then, we construct a function to compute QAC distance between two multivariate time series

# Input:
# X, Y: two MTS in matrix form

qaf_mts <- function(X, Y){
  c <- ncol(X)
  dist <- 0
  for (i in 1 : c) {
    dist <- dist + qac_uts(X[,i], Y[,i])
  }
  
  dist
  
}

# Now, we construct a function to extract the vector gammax from a univariate time series 

qaf_uts_coefs <- function(X, levels = c(0.1, 0.5, 0.9), l = 1) {
  
  k <- length(levels)
  nx <- length(X) 
  qx <- quantile(X, probs = levels)
  gammax <- matrix(0, k, k) # Vector gamma for X, regarding lag 1
  
  for (i in 1 : k) {
    for (j in 1 : k) {
      ax <- (X <= qx[i])[1 : (nx - l)]
      bx <- (X[(1 + l) : nx] <= qx[j])
      gammax[i, j] = (1/(nx - l))*sum(ax * bx)-levels[i]*levels[j]
    }
  }
  
  as.vector(gammax)
  
}

# Now, we construct a function to extract the vector of each gammax from a multivariate time series 

qaf_mts_coefs <- function(X){
  c <- ncol(X)
  gammax <- numeric()
  for (i in 1 : c) {
    gamma <- qaf_uts_coefs(X[,i])
    gammax <- c(gammax, gamma)
  }
  
  as.vector(gammax)
  
}

# We generalize computation of gamma vector 


qaf_uts_coefs_xy <- function(X, Y, levels = c(0.1, 0.5, 0.9), l = 1) {
  
  k <- length(levels)
  nx <- length(X)
  qx <- quantile(X, probs = levels)
  gamma <- matrix(0, k, k) # Vector gamma for X, regarding lag 1
  
  for (i in 1 : k) {
    for (j in 1 : k) {
      ax <- (X <= qx[i])[1 : (nx - l)]
      bx <- (Y[(1 + l) : nx] <= qx[j])
      gamma[i, j] = (1/(nx - l))*sum(ax * bx)-levels[i]*levels[j]
    }
  }
  
  as.vector(gamma)
  
}


qaf_mts_coefs_xy <- function(X){
  c <- ncol(X)
  gamma <- matrix(ncol = 9) # 9 if l = 1 and there are 3 probability levels 
  for (i in 1 : c) {
    for (j in 1 : c){
    gammanew <- qaf_uts_coefs_xy((X[,i]-mean(X[,i])), (X[,j]-mean(X[,j])))
    gamma <- rbind(gamma, gammanew)
    }
  }
  gamma <- gamma[-1,]
  as.vector(gamma)
  
}


# By considering separately the quantiles of X and Y


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
