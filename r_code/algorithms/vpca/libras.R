

# Preparation of dataset Libras 

libras <- read.csv('libras.txt', header = F)
colnames(libras) <- NULL


# We have to obtain 360 MTS of length 45, each one with a label from 1 to 15 (there are 15 clusters)

S <- vector(mode = "list", length = nrow(libras))

# S is going to contain each MTS with its label

for (i in 1 : length(S)) {
  # S[[i]] <- pendigits[i,] 
  S1 <- pendigits[i, c(seq(1, 90, by = 2), 91)] 
  S2 <- pendigits[i, c(seq(2, 90, by = 2), 91)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

save(S, file = 'libras.RData')


# Fuzzy C-means and clustering validation

# Removing the labels from S


M <- vector(mode = "list", length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- S[[i]][,seq(1, 45)]
}

# Dimensionality reductions vpca

Y <- vpca(M, lambda = 0.95)

# Implementation of fuzzy c-means algorithm

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzy_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/sw_distance_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/vpca_function.R')

u <- fcm_mts(Y = Y, K = 15, niter = 50, b = 2, tol = 0.01 )


c  <- matrix(nrow = nrow(u), ncol = ncol(u))

for (i in 1 : nrow(c)) {
  for (j in 1 : ncol(c)) {
    c[i, j] = u[i, j]/max(u[,j])
  }
}

# Converting the matrix in a 1-0 matrix

c[c != 1] <- 0

# Converting the results to a vector, in order to compute clustering validity indexes

clustering <- numeric(ncol(c))
for (j in 1 : length(clustering)) {
  clustering[j] <- which.max(c[,j])
}

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- S[[j]][1, 46]
}

library(dtwclust)
cvi(ground_truth, clustering)

