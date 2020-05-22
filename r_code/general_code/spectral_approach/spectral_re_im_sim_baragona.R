

# Lets compare the performance of quantile coherence by using the square modulus vs by using real and imaginary parts separately 

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
  cluster12[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster22[[i]] <- VARMAsim(100, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 100) {
  cluster32[[i]] <- VARMAsim(100, arlags = 1, malags = 1, phi = phi_c3,
                             theta = theta_c3, sigma = sigma_c3)$series
}

experiment2 <- c(cluster12, cluster22, cluster32)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))

# K means square modulus

coherences <- listTomatrix(lapply(experiment2, quantile_coherence))

b <- numeric()
for (i in 1: 400){
clustering <- kmeans(coherences, 3)$cluster
b[i] <- external_validation(ground_truth, clustering)
}
mean(b)

# K means real-imaginary

coherences <- listTomatrix(lapply(experiment2, quantile_coherence_re_im))

b <- numeric()
for (i in 1: 400){
clustering <- kmeans(coherences, 3)$cluster
b[i] <- external_validation(ground_truth, clustering)
}
mean(b)

# K means quantile autocovariances 

gamma <- listTomatrix(lapply(experiment2, qaf_mts_coefs_xy_sep))

b <- numeric()
for (i in 1: 400){
  clustering <- kmeans(gamma, 3)$cluster
  b[i] <- external_validation(ground_truth, clustering)
}
mean(b)
