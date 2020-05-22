

# Lets check the classical fuzzy-c means function with the iris data set

# Loading the data 

data('iris')
X <- iris[, 3:4] # Only pepal length and petal width are going to be used to perform clustering analysis 

# Plotting the data 

library(ggplot2)
ggplot(X, aes(x = iris$Petal.Length, y = iris$Petal.Width, colour = iris$Species)) + geom_point() # We have three clusters


# Lets perform the clustering algorithm

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/algorithms/classical_fuzzy/classical_fuzzy_function.R')

# Parameters:
# Number of clusters, K = 3
# Fuzziness coefficient, b = 2
# Maximum number of iterations, niter = 1000
# Threshold for terminations step, tol = 0.01

mem_matrix <- fuzzy_cm(X, 3, 2, 1000, 0.0001)

# Converting the membership matrix into a 1-0 matrix 

mem_matrix_crisp <- matrix(nrow = nrow(mem_matrix), ncol = ncol(mem_matrix))

for (i in 1 : nrow(mem_matrix)) {
  for (j in 1 : ncol(mem_matrix)) {
    mem_matrix_crisp[i, j] = mem_matrix[i, j]/max(mem_matrix[,j])
  }
}

# Converting the matrix in a 1-0 matrix

mem_matrix_crisp[mem_matrix_crisp != 1] <- 0

# Converting the results to a vector, in order to compute clustering validity indexes

clustering <- numeric(ncol(mem_matrix_crisp))
for (j in 1 : length(clustering)) {
  clustering[j] <- which.max(mem_matrix_crisp[,j])
}


ground_truth <- numeric(nrow(iris))

for (j in 1 : nrow(iris)) {
  ground_truth[j] <- iris[j,5]
}

library(dtwclust)
cvi(ground_truth, clustering)

# Plotting the clustering solution of fuzzy c means

data_fuzzy <- cbind(X, clustering)

ggplot(data_fuzzy, aes(x = data_fuzzy$Petal.Length, y = data_fuzzy$Petal.Width, 
                       colour = data_fuzzy$clustering)) + geom_point() # We have three clusters

