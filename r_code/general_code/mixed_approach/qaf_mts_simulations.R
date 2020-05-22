

# Lets perfom simulations in the second scenario of Baragona's paper, and applying hierarchical clustering with QAF. Then,
# we are going to add simulated multivariate GARCH processes 

# Simulations

# Second scenario
# Cluster 1: VAR(1) model, 3 variables, 100 time points, 100 elements
# Cluster 2: VMA(1) model, 3 variables, 100 time points, 100 elements
# Cluster 3: VARMA(1) model, 3 variables, 100 time points, 100 elements

cluster12 <- list()
cluster22 <- list()
cluster32 <- list()

phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c3 <- matrix(c(0.5, -0.3, 0.0, -0.3, 0.5, -0.4, 0.0, 0.4, 0.5), nrow = 3)
theta_c3 <- matrix(c(-0.5, 0.4, 0.0, 0.4, -0.5, -0.3, 0.0, 0.3, -0.5), nrow = 3)
sigma_c3 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)

set.seed(1234)
for (i in 1 : 100) {
  cluster12[[i]] <- VARMAsim(45, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster22[[i]] <- VARMAsim(45, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 100) {
  cluster32[[i]] <- VARMAsim(45, arlags = 1, malags = 1, phi = phi_c3,
                             theta = theta_c3, sigma = sigma_c3)$series
}

experiment2 <- c(cluster12, cluster22, cluster32)
ground_truth2 <- c(rep(1, 100), rep(2, 100), rep(3, 100))


# Computing dissimilarity matrix 


gamma <- listTomatrix(lapply(experiment2, qaf_mts_coefs_xy_sep))
dis_matrix <- matrix(0, 300, 300)
for (i in 1 : 300) {
  for (j in 1 : 300) {
    dis_matrix[i, j] <- EuclideanDistance(gamma[i,], gamma[j,])^2
  }
}

# Performing hierarchical clustering

hierarchical <- hclust(dist(dis_matrix))
clustering <- cutree(hierarchical, 3)
plot(hierarchical)
external_validation(ground_truth2, clustering)

# Performing k-means new approach 

clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth2, clustering)


# Multivariate GARCH processes (volatilities)

library(mgarchBEKK)
cluster1 <- list()
cluster2 <- list()

set.seed(123)
for (i in 1 : 100){
mgarch1 <- simulateBEKK(3, 50, c(1,2))
cluster1[[i]] <- cbind(mgarch1$sd[[1]], mgarch1$sd[[2]], mgarch1$sd[[3]])
}

for (i in 1 : 100){
  mgarch2 <- simulateBEKK(3, 50, c(2,1))
  cluster2[[i]] <- cbind(mgarch2$sd[[1]], mgarch2$sd[[2]], mgarch2$sd[[3]])
}


ground_truth <- c(rep(1, 100), rep(2, 100))
M <- c(cluster1, cluster2)

# Performing k-means 

gamma <- listTomatrix(lapply(M, qaf_mts_coefs))
clustering <- kmeans(gamma, 2)$cluster
external_validation(ground_truth, clustering)


# Multivariate GARCH processes (processes)


library(mgarchBEKK)
library(MASS)
cluster_sigma1 <- list()
cluster_sigma2 <- list()
epsilon1 <- list()
epsilon2 <- list()
cluster1 <- list()
cluster2 <- list()
theta_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7), nrow = 3)

for (i in 1 : 100) {
  epsilon1[[i]] <- mvrnorm(n = 700, mu = c(0, 0, 0), Sigma = diag(3), tol = 1e-6, empirical = FALSE)
}

for (i in 1 : 100) {
  epsilon2[[i]] <- mvrnorm(n = 700, mu = c(0, 0, 0), Sigma = diag(3), tol = 1e-6, empirical = FALSE)
}

set.seed(123)
for (i in 1 : 100){
  mgarch1 <- simulateBEKK(3, 700, c(1,2))
  cluster_sigma1[[i]] <- cbind(mgarch1$sd[[1]], mgarch1$sd[[2]], mgarch1$sd[[3]])
}

for (i in 1 : 100){
  mgarch2 <- simulateBEKK(3, 700, c(2,1))
  cluster_sigma2[[i]] <- cbind(mgarch2$sd[[1]], mgarch2$sd[[2]], mgarch2$sd[[3]])
}

for (i in 1 : 100) {
  cluster1[[i]] <- VARMAsim(700, malags = 1, theta = theta_c1, sigma = diag(3))$series + cluster_sigma1[[i]]*epsilon1[[i]]
}

for (i in 1 : 100) {
  cluster2[[i]] <-  VARMAsim(700, malags = 1, theta = theta_c1, sigma = diag(3))$series + cluster_sigma2[[i]]*epsilon2[[i]]
}


ground_truth <- c(rep(1, 100), rep(2, 100))
M <- c(cluster1, cluster2)

# Performing k-means 

gamma <- listTomatrix(lapply(M, qaf_mts_coefs_xy_sep))
clustering <- kmeans(gamma, 2)$cluster
external_validation(ground_truth, clustering)



# Small experiment to compare naive QAF approach vs whole QAF approach

cluster12 <- list()
cluster22 <- list()

phi_c1 <- matrix(c(0.9, 0.3, 0.3, 0.1), nrow = 2)
sigma_c1 <- matrix(c(1, 0.25, 0.25, 0.25), nrow = 2)
phi_c2 <- matrix(c(0.1, 0.2, 0.2, 0.9), nrow = 2)
sigma_c2 <- matrix(c(1, 0.3, 0.3, 0.25), nrow = 2)

set.seed(1234)
for (i in 1 : 100) {
  cluster12[[i]] <- VARMAsim(30, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster22[[i]] <- VARMAsim(30, arlags  = 1, phi = phi_c2, sigma = sigma_c2)$series
}



experiment <- c(cluster12, cluster22)
ground_truth <- c(rep(1, 100), rep(2, 100))

gamma <- listTomatrix(lapply(experiment, qaf_mts_coefs))

b <- numeric(200)
for (i in 1 : 200){
  clustering <- kmeans(gamma, 2)$cluster
  b[i] <- external_validation(ground_truth, clustering)
}
mean(b)


# Testing QAF with dimensionality reduction via CPCA, with series of 6 variables

varmalist <- list()
phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7, -0.5, 0, 0.5, 0.3, -0.4, 0.3, 0.5, -0.2, 0.8, 0, 
                   -0.5, 0.4, -0.3, 0, 0, 0.4, 0.2, 0.5, -0.7, 0.2, 0, 0.5, 0.5, 0.5, 0, -0.6, 0.2), nrow = 6)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 0.10, 0.25, 0.25, 1, 0.25, 0.10, 0.25, 0.10, 0.10, 0.25, 1, 0.10, 0, 0.25, 0.25,
                     0.10, 0.10, 1, 0.10, 0.25, 0.10, 0.25, 0, 0.10, 1, 0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 1), nrow = 6)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7, 0.5, 0.6, 0.4, -0.6, 0, 0.2, 0.5, 0.4, 0.1, 0, 0.4,
                     0.4, 0.5, 0, 0.5, 0.2, 0, 0.7, 0.1, 0, 0.3, 0.4, 0.5, 0.6, 0, 0.2, 0.1), nrow = 6)
sigma_c2 <- sigma_c1
phi_c3 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.8, 0.5, 0.0, -0.3, 0.7, -0.5, 0, 0.5, 0.3, -0.4, -0.3, 0.5, -0.2, 0.8, 0, 
                   -0.5, 0.4, -0.3, 0, 0, -0.4, 0.2, 0.5, -0.7, 0.2, 0, 0.5, 0.5, -0.5, 0, -0.6, 0.2), nrow = 6)
theta_c3 <-  matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7, 0.5, 0.6, 0.4, -0.6, 0, 0.2, 0.5, 0.4, 0.1, 0, 0.4,
                      0.4, 0.5, 0, 0.5, 0.2, 0, 0.7, 0.1, 0, 0.3, 0.4, 0.5, 0.6, 0, 0.2, 0.1), nrow = 6)
sigma_c3 <- 2*sigma_c2

set.seed(1234)
for (i in 1 : 100) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 101 : 200) {
  varmalist[[i]] <- VARMAsim(100, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 201 : 300) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, malags = 1, phi = phi_c3, theta = theta_c3, sigma = sigma_c3)$series
}


sigma <- list()

for (i in 1:300) {
  sigma[[i]] <- cov(varmalist[[i]])
}

s <- cpca(sigma, lambda = 0.60)
varmalist_reduced <- list()

for (i in 1 : 300) {
  varmalist_reduced[[i]] <- varmalist[[i]] %*% s
}

# Performing k-means 

ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))
gamma <- listTomatrix(lapply(varmalist_reduced, qaf_mts_coefs_xy_sep))
clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth, clustering)


