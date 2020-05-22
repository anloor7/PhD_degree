

# We are going to perform the simulations presented in Bandyopadhyay, Baragona and Maulik

library(MTS)




# First scenario
# Cluster 1: VAR(1) model, 2 variables, 100 time points, 334 elements
# Cluster 2: VAR(1) model, 3 variables, 100 time points, 334 elements
# Cluster 3: VAR(1) model, 4 variables, 100 time points, 334 elements

cluster11 <- list()
cluster21 <- list()
cluster31 <- list()

phi_c1 <- matrix(c(0.8, -0.4, 0.7, 0.6), nrow = 2)
sigma_c1 <- matrix(c(1, 0.25, 0.25, 1), nrow = 2)
phi_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c3 <- matrix(c(0.6, -0.4, 0.1, 0.0, 0.5, 0.5, -0.5, 0.1, 0.1, 0.3, 0.7, 0.5,
                   0.0, 0.1, 0.1, 0.6), nrow = 4)
sigma_c3 <- matrix(c(1, 0.25, 0.10, 0.0, 0.25, 1.0, 0.25, 0.1, 0.1, 0.25, 1, 0.25,
                   0.0, 0.1, 0.25, 1), nrow = 4)

for (i in 1 : 334) {
  cluster11[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 334) {
  cluster21[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

for (i in 1 : 334) {
  cluster31[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c3, sigma = sigma_c3)$series
}

experiment1 <- c(cluster11, cluster21, cluster31)
ground_truth1 <- c(rep(1, 334), rep(2, 334), rep(3, 334))

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
ground_truth2 <- c(rep(1, 100), rep(2, 100), rep(3, 100))

# Third scenario
# Cluster 1: VAR(1) model, 4 variables, 100 time points, 334 elements 
# Cluster 2: VMA(1) model, 3 variables, 100 time points, 334 elements
# Cluster 3: VAR(2) model, 2 variables, 100 time points, 334 elements 

cluster13 <- list()
cluster23 <- list()
cluster33 <- list()

phi_c1 <- matrix(c(0.6, -0.4, 0.1, 0.0, 0.5, 0.5, -0.5, 0.1, 0.1, 0.3, 0.7, 0.5,
                   0.0, 0.1, 0.1, 0.6), nrow = 4)
sigma_c1 <- matrix(c(1, 0.25, 0.10, 0.0, 0.25, 1.0, 0.25, 0.1, 0.1, 0.25, 1, 0.25,
                     0.0, 0.1, 0.25, 1), nrow = 4)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi1_c3 <- matrix(c(0.5, -0.1, 0.1, -0.5), nrow = 2)
phi2_c3 <- matrix(c(0.3, -0.1, 0.1, -0.5), nrow = 2)
sigma_c3 <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)

for (i in 1 : 334) {
  cluster13[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 334) {
  cluster23[[i]] <- VARMAsim(100, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 1 : 334) {
  cluster33[[i]] <- VARMAsim(100, arlags = 2, phi = cbind(phi1_c3, phi2_c3),
                             sigma = sigma_c3)$series
}

experiment3 <- c(cluster13, cluster23, cluster33)
ground_truth3 <- c(rep(1, 334), rep(2, 334), rep(3, 334))
