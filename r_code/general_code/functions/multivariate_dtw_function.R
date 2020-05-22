

# This is a function to compute DTW between two MTSs

# Input:
# X, Y: two multivariate MTSs in matrix form (time points by rows and variables by columns)

dtw_mts_i <- function(X, Y){
  c <- ncol(X)
  distance <- 0
  
  for ( i in 1 : c) {
    distance <- distance + dtw(X[,i], Y[,i])$normalizedDistance
  }
  
  return(distance)
  
}


