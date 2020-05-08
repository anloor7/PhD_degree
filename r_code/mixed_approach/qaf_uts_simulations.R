

# Lets perfom simulations in one scenario of Borja's Thesis (page 58), and applying hierarchical clustering with QAF

# Simulations

models <- matrix(nrow = 300, ncol = 100)

for (i in 1 : 100){ # First model 
 models[i,] <- arima.sim(n = 100, list(ar = 0.9))

 }

for (i in 101 : 200){ # Second model 
  models[i,] <- arima.sim(n = 100, list(ma = -0.7))
  
 }

for (i in 201 : 300){ # Third model 
  models[i,] <- arima.sim(n = 100, list(ar = c(0.3, -0.1)))
  
 }

ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))

# Computing dissimilarity matrix

gamma <- t(apply(models, 1, qaf_uts_coefs))
dis_matrix <- matrix(0, 300, 300)
for (i in 1 : 300) {
 for (j in 1 : 300) {
   dis_matrix[i, j] <- EuclideanDistance(gamma[i,], gamma[j,])
 }
}

# Performing hierarchical clustering 

hierarchical <- hclust(dist(dis_matrix))
clustering <- cutree(hierarchical, 3)
plot(hierarchical)
external_validation(ground_truth, clustering)

# Performing k-means 

clustering <- kmeans(gamma, 3)$cluster
external_validation(ground_truth, clustering)




