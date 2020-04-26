

# Lets build a function to convert a fuzzy membership matrix to a clustering vector via the naive criterion

fuzzytocrisp <- function(u){
  
  n <- nrow(u)
  c <- ncol(u)
  w <- matrix(nrow = n, ncol = c)
  for (i in 1 : n) {
    for (j in 1 : c) {
      w[i, j] = u[i, j]/max(u[,j])
    }
  }
  
  w[w != 1] <- 0
  clustering <- numeric(c)
  for (j in 1 : c) {
    clustering[j] <- which.max(w[,j])
  }
  
  clustering
}

