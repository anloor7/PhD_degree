

# Lets check the performance of our kmeans function

library(mvtnorm)

# Simulating 50 observations from a different multivariate distribution

set.seed(1234)
X1 <- rmvnorm(50, mean = c(3, 3), sigma = diag(2)/3)
X2 <- rmvnorm(50, mean = c(0, 0), sigma = diag(2)/3)
X3 <- rmvnorm(50, mean = c(-3, -3), sigma = diag(2)/3)

X <- rbind(X1, X2, X3)
X <- as.data.frame(X)
X$color <- numeric(150)
X$color[1 :50] <- 'red'
X$color[51:100] <- 'green'
X$color[101:150] <- 'blue'

plot(X[,1], X[,2], col = X$color)

Z <- X[,1:2]


# K-means algorithm

library(TSdist)

# Converting Y to a list

Y <- vector(mode = 'list', length = 150)

for (i in 1 : 150) {
  Y[[i]] <- as.numeric(Z[i,])
}

c <- km_mts(Y = Y, K = 3, dis = EuclideanDistance)
clustering <- fuzzytocrisp(c)

# Clustering evaluation

ground_truth <- c(rep(1, 50), rep(2, 50), rep(3, 50))
library(dtwclust)
cvi(ground_truth, clustering)

# The same example using fuzzy-c-means

c <- fcm_mts(Y = Y, K = 3, b = 2, dis = EuclideanDistance)
clustering <- fuzzytocrisp(c)

# Clustering evaluation

ground_truth <- c(rep(1, 50), rep(2, 50), rep(3, 50))
library(dtwclust)
cvi(ground_truth, clustering)

