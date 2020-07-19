

# This is a function to compute the validity index based on the one nearest neighbour classifier evaluated by
# leave one out cross validation. It ranges from 0 to 1

# Input:
# D: distance matrix
# clustering: clustering solution

loo1nn <- function(D, c){
  
  n <- nrow(D)
  index <- numeric(n)
  D <- as.matrix(D)
  
  for (i in (1 : n)) {
    
    
    j <- order(D[i,])[2]
    
   
    
    if (c[i] == c[j]){
      
      index[i] <- 1
      
    }
    
    
  }
  
  sum(index)/n
  
}
