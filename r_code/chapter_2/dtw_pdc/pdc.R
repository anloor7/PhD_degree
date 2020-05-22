

# Cluster creation

cluster1 <- list()
cluster2 <- list()
phi_c1 <- matrix(c(0.3, 0, 0.2, 0.2, 0.4, 0, 0, 0.3, 0.1), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c2 <- matrix(c(0.4, 0, 0.3, 0.2, 0.4, 0, 0, 0.4, 0.2), nrow = 3)
sigma_c2 <- sigma_c1

set.seed(1234)
for (i in 1 : 100) {
  cluster1[[i]] <- VARMAsim(400, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 100) {
  cluster2[[i]] <- VARMAsim(400, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

cluster <- c(cluster1, cluster2)

# Ground truth 

ground_truth <- c(rep(1, 100), rep(2, 100))


# Distance matrix

series_matrix <- t(listTomatrix(cluster))
dis_matrix <- pdcDist(series_matrix)

# PAM 

clustering <- pam(dis_matrix, 2)$clustering
clustering
external_validation(ground_truth, clustering)

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
clustering1 <- kmeans(gamma, 2)$cluster
external_validation(ground_truth, clustering1)




# PDC with univariate series 


cluster1 <- list()
cluster2 <- list()


set.seed(1234)
for (i in 1 : 100) {
  cluster1[[i]] <- arima.sim(n = 300, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
                             sd = sqrt(0.1796))
}

for (i in 1 : 100) {
  cluster2[[i]] <- arima.sim(n = 300, list(ar = c(0.2897, -0.6858), ma = c(0.2279, 0.6488)),
                             sd = sqrt(0.1796))
}

cluster <- c(cluster1, cluster2)


series_matrix <- t(listTomatrix(cluster))
dis_matrix <- pdcDist(series_matrix)

clustering <- pam(dis_matrix, 2)$clustering
external_validation(ground_truth, clustering)
