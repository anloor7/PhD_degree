

# Preparation of dataset Libras 

libras <- read.csv('libras.txt', header = F)
colnames(libras) <- NULL


# We have to obtain 360 MTS of length 45, each one with a label from 1 to 15 (there are 15 clusters)

S <- vector(mode = "list", length = nrow(libras))

# S is going to contain each MTS with its label

for (i in 1 : length(S)) {
  # S[[i]] <- pendigits[i,] 
  S1 <- libras[i, c(seq(1, 90, by = 2), 91)] 
  S2 <- libras[i, c(seq(2, 90, by = 2), 91)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

save(S, file = 'libras.RData')


# Fuzzy C-means and clustering validation

# Removing the labels from S


M <- vector(mode = "list", length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- S[[i]][,seq(1, 45)]
}

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzy_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/sw_distance_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/vpca_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzytocrisp_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/kmeans_function.R')

# Dimensionality reduction vpca

Y <- vpca(M, lambda = 0.95)

# Implementation of fuzzy c-means algorithm

u <- fcm_mts(Y = Y, K = 15, niter = 15, b = 2, tol = 0.01, dis = dist_matrix)
clustering <- fuzzytocrisp(u)

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(TSdist)
library(ClusterR)
library(dtwclust)
cvi(ground_truth, clustering)
external_validation(ground_truth, clustering, summary_stats = T)


# Implementation of K-means algorithm

Y <- vpca(M, lambda = 0.95)
u <- km_mts(Y = Y, K = 15, niter = 40, tol = 0.01, dis = dist_matrix)
clustering <- fuzzytocrisp(u)

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(dtwclust)
cvi(ground_truth, clustering)
external_validation(ground_truth, clustering, summary_stats = T)

# Implementation of K-means with another distance


dist_matrix <- function(A, B){
  A <- as.vector(A)
  B <- as.vector(B)
  EuclideanDistance(A, B)
}

u <- km_mts(Y = Y, K = 15, niter = 60, tol = 0.01, dis = dist_matrix)
clustering <- fuzzytocrisp(u)

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(dtwclust)
cvi(ground_truth, clustering)
external_validation(ground_truth, clustering, summary_stats = T)

# Implementation of K-means with another distance

dist_pca <- function(A, B){
  n <- nrow(A)
  c <- ncol(A)
  props <- rev((1 : c)/sum(1:c))
  vector <- numeric()
  for (i in 1 : c) {
    vector <- c(vector, rep(props[i], n))
  }
  a <- as.vector(A)
  b <- as.vector(B)
  sqrt(sum(vector * (a - b)^2)/sum(vector))
}

u <- km_mts(Y = Y, K = 15, niter = 200, tol = 0.01, dis = dist_pca)
clustering <- fuzzytocrisp(u)

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(dtwclust)
cvi(ground_truth, clustering)
external_validation(ground_truth, clustering, summary_stats = T)

# Implementation of fuzzy C-means with this distance 

u <- fcm_mts(Y = Y, K = 15, b = 2, niter = 10, tol = 0.01, dis = dist_pca)
clustering <- fuzzytocrisp(u)

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(dtwclust)
cvi(ground_truth, clustering)


