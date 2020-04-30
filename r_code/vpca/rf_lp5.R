

# Preparation of dataset Robot Failure LP5

lp5 <- read.csv('lp5.data.txt', header = F, sep = '\t')
head(lp5)


# We have to obtain  164 MTS of length 15, each one with a given label (there are 5 labels)

S <- vector(mode = "list", length = 164)

for (i in 0 : 163) {
  S[[i + 1]] <- matrix(nrow = 16, ncol = 7)
  S[[i + 1]] <- lp5[ (16 * i + 1): (16 * (i + 1)),]
}

for (i in 1 : 164) {
  label <- as.numeric(S[[i]][1, 1])
  S[[i]] <- cbind(S[[i]], rep(label, 16))
  colnames(S[[i]])[8] <- 'label'
  S[[i]] <- S[[i]][2:16, 2:8]
  S[[i]] <- t(S[[i]])
  rownames(S[[i]]) <- NULL
  colnames(S[[i]]) <- NULL
}

# Constructing the ground truth 

ground_truth <- numeric(164)

for (j in 1: 164) {
  ground_truth[j] <- S[[j]][7, 1]
}

# Constructing the MTS dataset without labels 

M <- vector(mode = "list", length = 164)
for (i in 1 : 164) {
  M[[i]] <- matrix(nrow = 6, ncol = 15)
  M[[i]] <- S[[i]][1:6,]
}


source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzy_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/sw_distance_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/vpca_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzytocrisp_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/kmeans_function.R')


# Implementation of K-means algorithm with VPCA 

library(TSdist)
Y <- vpca(M, lambda = 0.99)
u <- km_mts(Y = Y, K = 5, niter = 1000, tol = 0.01, dis = dist_pca)
clustering <- fuzzytocrisp(u)

library(ClusterR)
library(dtwclust)
cvi(ground_truth, clustering)
external_validation(ground_truth, clustering, summary_stats = T)






