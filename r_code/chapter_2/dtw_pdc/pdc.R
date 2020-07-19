

# Cluster creation

cluster1 <- list()
cluster2 <- list()
phi_c1 <- matrix(c(0.3, 0, 0.2, 0.2, 0.4, 0, 0, 0.3, 0.1), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c2 <- matrix(c(0.5, 0, 0.7, 0.2, -0.5, 0, 0, 0.7, 0.2), nrow = 3)
sigma_c2 <- sigma_c1

set.seed(1234)
for (i in 1 : 10) {
  cluster1[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 1 : 10) {
  cluster2[[i]] <- VARMAsim(200, arlags = 1, phi = phi_c2, sigma = sigma_c2)$series
}

cluster <- c(cluster1, cluster2)

# Ground truth 

ground_truth <- c(rep(1, 10), rep(2, 10))


# Distance matrix

series_matrix <- t(listTomatrix(cluster))
dis_matrix <- pdcDist(series_matrix)

# PAM 

clustering <- pam(dis_matrix, 2)$cluster
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



# Proofs hierarchical 

# Computing first the distance matrix

series_matrix <- t(listTomatrix(cluster))
dis_matrix <- pdcDist(series_matrix)
hierarchical <- hclust(dis_matrix, method = 'complete')
clustering <- cutree(hierarchical, 2)
external_validation(ground_truth, clustering)

# Option 2 

series_matrix1 <- series_matrix[1 : nrow(cluster[[1]]),]
series_matrix2 <- series_matrix[((nrow(cluster[[1]]) + 1) : (2 * nrow(cluster[[1]]))),]
series_matrix3 <- series_matrix[((2 * nrow(cluster[[1]]) + 1) : (3 * nrow(cluster[[1]]))),]
dis_matrix1 <- pdcDist(series_matrix1)
dis_matrix2 <- pdcDist(series_matrix2)
dis_matrix3 <- pdcDist(series_matrix3)
dis_matrix123 <- dis_matrix1 + dis_matrix2 + dis_matrix3
hierarchical <- hclust(dis_matrix123, method = 'complete')
clustering <- cutree(hierarchical, 2)
external_validation(ground_truth, clustering)

# Directly 

array_series <- array(dim = c(200, 20, 3))

for (i in 1 : 20) {
  array_series[,i,] <- cluster[[i]]
}


direct_clust <- pdclust(array_series, clustering.method = 'complete')
clustering1 <- cutree(direct_clust, 2)
external_validation(ground_truth, clustering1)

dis_matrix <- pdcDist(array_series)
hierarchical <- pam(dis_matrix, 2)
external_validation(ground_truth, clustering)


# Computing first the distance matrix. Second approach. 

h <- matrix(0, nrow = 3 * 170, ncol = 20)
for (i in 1 : length(cluster)) {
  
  h[,i] <- as.vector(t(cluster[[i]]))
  
}

dis_matrix2 <- pdcDist(h)
hierarchical2 <- hclust(dis_matrix2, method = 'complete')
clustering2 <- cutree(hierarchical2, 2)
external_validation(ground_truth, clustering2)



# Computing pdc dist between each pair of UTS within the mts


pdc_dist_xy <- function(X, Y){
  
  n <- nrow(X)
  c <- ncol(X)
  
  pdcdist <- numeric()
  
  
  for (i in 1 : c) {
    for (j in 1 : c) {
      
      pdcdist <- c(pdcdist, pdcDist(cbind(X[,i], Y[,j])))
      
    }
    
  }
  
  sum(pdcdist)
  
}


dis_matrix <- matrix(0, nrow = 20, ncol = 20)

for (i in 1 : 20) {
  for (j in 1 : 20) {
    dis_matrix[i, j] <- pdc_dist_xy(cluster[[i]], cluster[[j]])
  }
  
}


diag(dis_matrix) <- 0
hierarchical <- hclust(dist(dis_matrix), method = 'complete')
clustering <- cutree(hierarchical, 2)
external_validation(ground_truth, clustering)
