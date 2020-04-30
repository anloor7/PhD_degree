

# Lets perform some analyses regarding distance between varma. The first is going to be based on simulations, and the second, on the 
# Libras dataset 

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/kmeans_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/varma_coef.R')

# Simulation of 20 VARMA(1, 1) processes stemming from two different model generators

X <- list()

theta1 <- matrix(c(0.5, 0.6, 0.6, 0.5), nrow = 2)
phi1 <- matrix(c(0.2, 0, 0.3, 0.7), nrow = 2)
epsilon1 <- diag(2)

for (i in 1 : 10){
X[[i]] <- VARMAsim(150, arlags = c(1), malags = c(1), phi = phi1, theta = theta1, sigma = epsilon1)$series}


theta2 <- matrix(c(1, 0.4, 0.4, 0.8), nrow = 2)
phi2 <- matrix(c(0.5, 0.2, 0.8, 0.7), nrow = 2)
epsilon2 <- diag(2)

for (i in 11 : 20){
  X[[i]] <- VARMAsim(150, arlags = c(1), malags = c(1), phi = phi2, theta = theta2, sigma = epsilon2)$series}

ground_truth <- c(rep(1, 10), rep(2, 10))

# To perform clustering, we can work with the corresponding vector of coefficients and the Euclidean Distance

v_c <- list()
vc <- varma_coefs(X)

clustering <- km_mts(vc, K = 2, niter = 100, tol = 0.01, dis = EuclideanDistance)
external_validation(ground_truth, clustering, summary_stats = T)


# Dataset Libras 

# Loading the data 

libras <- load('libras.RData')
M <- list(length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- t(S[[i]])
}

vc <- varma_coefs(M)
clustering <- km_mts(vc, K = 15, niter = 300, tol = 0.01, dis = EuclideanDistance)
external_validation(ground_truth, clustering, summary_stats = T)
