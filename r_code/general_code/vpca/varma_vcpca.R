

##### PERFORMING VPCA_SWMDFC WITH SOME SIMULATED VAR MODELS 

# Lets simulate MTS following VARMA models under the following conditions:
# N = 2: number of variables
# L = 200: number of time points
# M = 90: number of time series
# K = 3: number of clusters, 30 time series within each cluster 

# Simulating three VAR processes, each one with different parameters and different mean vector 

library(MTS)
set.seed(1234)

# First group

Y1 <- vector(mode = 'list', length = 30)
ny1 <- length(Y1)

for (i in 1 : ny1) {
  Y1[[i]] = t(VARMAsim(200, arlags = 1, phi = matrix(c(0.9, 0, 0, 0.9), nrow = 2), 
                     sigma = matrix(c(0.0001, 0, 0, 0.0001), nrow = 2))$series)
 }

# Plotting some series in the first group 

par(mfrow = c(2, 5))

for (i in 1 : 10) {
  ts.plot(t(Y1[[i]]))
}

# Second group

Y2 <- vector(mode = 'list', length = 30)
ny2 <- length(Y2)

for (i in 1 : ny2) {
  Y2[[i]] = t(VARMAsim(200, cnst = c(0.3, 0.3), arlags = 1, phi = matrix(c(0.9, 0.7, 0, 0.9), nrow = 2), 
                     sigma = matrix(c(0.3, 0, 0, 0.3), nrow = 2))$series)
}

# Plotting some series in the second group 

par(mfrow = c(2, 5))

for (i in 1 : 10) {
  ts.plot(t(Y2[[i]]))
}

# Third group

Y3 <- vector(mode = 'list', length = 30)
ny3 <- length(Y3)

for (i in 1 : ny3) {
  Y3[[i]] = t(VARMAsim(200, arlags = 1, cnst = c(0.5, 0.5), phi = matrix(c(0.9, 0, 0, -0.9), nrow = 2), 
                     sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))$series)
}

# Plotting some series in the second group 

par(mfrow = c(2, 5))

for (i in 1 : 10) {
  ts.plot(t(Y3[[i]]))
}


# Creating a whole list

Yt <- c(Y1, Y2, Y3)

# Loading the required functions 

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzy_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/sw_distance_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/vpca_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/vpca/fuzzytocrisp_function.R')

# Dimensionality reduction

Y <- vpca(Yt)

# SWMDFC

u <- fcm_mts(Y, K = 3, b = 2)

# Evaluating clustering results

clustering <- fuzzytocrisp(u)
ground_truth <- c(rep(1, ny1), rep(2, ny2), rep(3, ny3))

library(dtwclust)
cvi(ground_truth, clustering)



# Now, lets repeat the same analysis, but now with the same mean vector within the three groups. The main
# difference relies on the underlying model structure

for (i in 1 : ny1) {
  Y1[[i]] = t(VARMAsim(200, arlags = 1, phi = matrix(c(0.9, 0, 0, 0.9), nrow = 2), 
                       sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))$series)
}

for (i in 1 : ny2) {
  Y2[[i]] = t(VARMAsim(200, arlags = 1, phi = matrix(c(0.9, 0, 0, -0.9), nrow = 2), 
                       sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))$series)
}

for (i in 1 : ny3) {
  Y3[[i]] = t(VARMAsim(200, arlags = 1, phi = matrix(c(0.9, 0.9, 0, 0.9), nrow = 2), 
                       sigma = matrix(c(0.1, 0, 0, 0.1), nrow = 2))$series)
}

# Creating a whole list

Yt <- c(Y1, Y2, Y3)

# Dimensionality reduction

Y <- vpca(Yt)

# SWMDFC

u <- fcm_mts(Y, K = 3, b = 2)

# Evaluating clustering results

clustering <- fuzzytocrisp(u)
ground_truth <- c(rep(1, ny1), rep(2, ny2), rep(3, ny3))

library(dtwclust)
cvi(ground_truth, clustering)
