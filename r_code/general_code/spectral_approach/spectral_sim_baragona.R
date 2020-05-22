

# Cluster analysis of scenario 2 of Baragona via Quantile Coherence 

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
  cluster12[[i]] <- VARMAsim(400, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster22[[i]] <- VARMAsim(400, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 100) {
  cluster32[[i]] <- VARMAsim(400, arlags = 1, malags = 1, phi = phi_c3,
                             theta = theta_c3, sigma = sigma_c3)$series
}

experiment2 <- c(cluster12, cluster22, cluster32)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))


# Applying K-means 

time <- system.time(coherences <- listTomatrix(lapply(experiment2, quantile_coherence))) # 6.6 seconds (100 time
# observations, 17.6 seconds, 400 time observations) 
clustering <- kmeans(coherences, 3)$cluster
external_validation(ground_truth, clustering)


# Cluster analysis of an escenario where each series has 6 variables 


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
  varmalist[[i]] <- VARMAsim(400, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 101 : 200) {
  varmalist[[i]] <- VARMAsim(400, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 201 : 300) {
  varmalist[[i]] <- VARMAsim(400, arlags = 1, malags = 1, phi = phi_c3, theta = theta_c3, sigma = sigma_c3)$series
}


# Applying K-means 

time <- system.time(coherences <- listTomatrix(lapply(varmalist, quantile_coherence))) #  37 seconds, 400 time points
clustering <- kmeans(coherences, 3)$cluster
external_validation(ground_truth, clustering)


# Applying K-means after dimensionality reduction 

sigma <- list()

for (i in 1:300) {
  sigma[[i]] <- cov(varmalist[[i]])
}

s <- cpca(sigma, lambda = 0.70)
varmalist_reduced <- list()

for (i in 1 : 300) {
  varmalist_reduced[[i]] <- varmalist[[i]] %*% s
}


gamma <- listTomatrix(lapply(varmalist_reduced, quantile_coherence))
clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth, clustering)



# Toy example: 500 MTS, 20 variables, 400 time points

M <- list()

for (i in 1 : 500) {
  M[[i]] <- mvrnorm(400, mu = rep(0, 20), Sigma = diag(20))
}

time <- system.time(coherences <- listTomatrix(lapply(M, quantile_coherence))) # 425 seconds 
time <- system.time(clustering <- kmeans(coherences, 5)$cluster) # 471 segundos 


# Parallelizing the computation 

n <- 7
c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('M', 'quantile_coherence'))


time <- system.time(v <- parLapply(c1, M, quantile_coherence))
time # 142 seconds
stopCluster(c1)

