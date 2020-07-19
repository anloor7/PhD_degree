

# Creating a function to perform k-means and obtaining the  ARI 
# Input:
# X: a dataset containing the elements to perform k-means 
# k: number of clusters  

kmeans_mc_av_ari_scenario1 <- function(X){
  k <- 3
  gt <- c(rep(1, 100), rep(2, 100), rep(3, 100))
  clustering <- kmeans(X, k)$cluster
  external_validation(gt, clustering)
}


kmeans_mc_av_ari_scenario2 <- function(X){
  k <- 4
  gt <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))
  clustering <- kmeans(X, k)$cluster
  external_validation(gt, clustering)
}


