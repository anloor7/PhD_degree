

########## PREPARATION OF DATASETS 


# We are going to load and prepare data  

# We start with the set of data called PenDigits 

pendigits <- read.csv('pendigits.txt', header = F)
colnames(pendigits) <- NULL


# We have to obtain 7494 MTS of length 8, each one with a label from 1 to 9 (there are 9 clusters)

P <- vector(mode = "list", length = nrow(pendigits))

# P is going to contain each MTS with its label

for (i in 1:length(P)) {
  # S[[i]] <- pendigits[i,] 
  P1 <- pendigits[i, seq(1, 17, by = 2)] 
  P2 <- pendigits[i, c(seq(2, 17, by = 2), 17)]
  P[[i]] <- rbind(as.numeric(P1), as.numeric(P2))
}

M <- vector(mode = "list", length = nrow(pendigits))

# M is going to contain each MTS without its label. 

for (i in 1:length(M)) {
  M[[i]] <- P[[i]][,seq(1, 8)]
}


# We can run the Fuzzy C-means algorithm by using M


# Obtaining both matrices V1 y V2, which are both 7494 x 8 matrices

V1 <- numeric()
V2 <- numeric()

for (i in 1:length(P)) {
  V1 <- rbind(V1, M[[i]][1,])
  V2 <- rbind(V2, M[[i]][2,])
}

# Dimensionality reduction via PCA

pc1 <- princomp(V1, cor = F)
summary(pc1)
pc2 <- princomp(V2, cor = F)
summary(pc2)

# In order to explain, in both instances, more than 90% of variability, we choose ps = 5 principal components 

# Computing the scores matrixes 

F1 <- pc1$scores[, 1:5]
F2 <- pc2$scores[, 1:5]

# Reconstructing the reduced sample of MTS time series

Y <- vector(mode = 'list', length = length(P))

for (i in 1 : length(P)) {
  Y[[i]] = matrix(nrow = 2, ncol = 5)
}

for (i in 1 : length(P)) {
  Y[[i]] = rbind(F1[i,], F2[i,])
}

# Computing S matrix (10 x 10 matrix)

S <- matrix(nrow = 10, ncol = 10)
constant <- 1/(2*pi*(1-(1/2))^2)
coefficients_previous <- numeric(100)

a <- 1
for (p1 in 1:5) {
  for (n1 in 1:2) {
    for (n2 in 1:2) {
      for (p2 in 1:5){
        coefficients_previous[a] <- sqrt((n1-n2)^2 + (p1 - p2)^2)
        a <- a + 1
      }
      
    }
  }
}

S_previous <- matrix(coefficients_previous, nrow = 10, ncol = 10, byrow = T) # It is almost the desired matrix, aside from the
# order of columns 

S_previous <- S_previous[,c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]
k <- 2*(1-1/5)^2
S <- 1/(pi * k) * exp(-S_previous^2/k)

# Computing the distance matrix between the 7494 initial MTS

y <- vector(mode = 'list', length = length(Y))
for (i in 1 : length(Y)) {
  y[[i]] <- as.vector(Y[[i]])
  
}

distances <- matrix(nrow = length(Y), ncol = length(Y))
for (i in 1 : length(Y)) {
  for (j in 1 : length(Y)) {
    distances[i, j] <- sqrt((y[[i]]-y[[j]]) %*% S %*% (y[[i]]-y[[j]]))
  }
}

# First proof: hierarchical clustering

# clustering <- hclust(dist(distances))
# plot(clustering)
# cutree(clustering, 10)


# Fuzzy C-means

K <- 10 # Hyperparameter, number of clusters 
b <- 2 # Hyperparameter, fuzzy Coefficient

# Selecting K centroids randomly 

vector_centroids <- sample(length(Y), K, replace = F)

# Computation of distance matrix

# First, we create a function to compute the given distance 

dist_fun <- function(x, y){
  sqrt((x-y) %*% S %*% (x-y))
}

dist_matrix <- matrix(nrow = K, ncol = length(Y)) # Initialization of distance matrix

# Initialization of centroids 

Z <- vector(mode = 'list', length = length(vector_centroids))
z <- vector(mode = 'list', length = length(vector_centroids))

for (i in 1 : length(vector_centroids)) {
  Z[[i]] <- Y[[vector_centroids[i]]]
  z[[i]] <- as.vector(Z[[i]])
}

niter <- 1 # Maximum number of iterations 
mem_matrix <- matrix(nrow = K, ncol = length(Y)) # Initialization of membership matrix 

# Fuzzy C-means algorithm

for (l in 1 : niter) {
  
  
  for (i in 1 : length(Y)) {
    for (k in 1 : length(Z)) {
      dist_matrix[k, i] <- dist_fun(y[[i]], z[[k]])
    }
  }
  for (i in 1 : length(Y)) {
    for (k in 1 : length(Z)) {
     
      sum_k <- sum(dist_matrix[,i]^(-2/(b-1)))
      mem_matrix[k, i] <- dist_matrix[k, i]^(-2/(b-1))/(sum_k)
      mem_matrix[mem_matrix == 'NaN'] <- 1
    }
  }
  for (k in 1 : length(Z)) {
    z[[k]] <- numeric(10)
    for (i in 1:length(Y)) {
      z[[k]] <- z[[k]] + mem_matrix[k, i]^b * y[[i]]
    }
    z[[k]] <- z[[k]]/sum(mem_matrix[k,]^b)
  }
  
  
}

head(mem_matrix)

# save(mem_matrix, file = 'mem_matrix_1_iter.RData')

# Lets perform hard assignment over the membership matrix

load('mem_matrix_1_iter.RData')
mem_matrix_1_iter <- mem_matrix


mem_matrix_trans2 <- matrix(nrow = nrow(mem_matrix_1_iter), ncol = ncol(mem_matrix_1_iter))

for (i in 1 : nrow(mem_matrix_1_iter)) {
  for (j in 1 : ncol(mem_matrix_1_iter)) {
    mem_matrix_trans2[i, j] = mem_matrix_1_iter[i, j]/max(mem_matrix_1_iter[,j])
  }
}

# Converting the matrix in a 1-0 matrix

mem_matrix_trans2[mem_matrix_trans2 != 1] <- 0

# Converting the results to a vector, in order to compute clustering validity indexes

clustering <- numeric(ncol(mem_matrix_trans2))
for (j in 1 : length(clustering)) {
  clustering[j] <- which.max(mem_matrix_trans2[,j])
}

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- P[[j]][1, 9]
}

library(dtwclust)
cvi(ground_truth, clustering)

# Results with  2 iterations

load('mem_matrix_2_iter.RData')
mem_matrix_2_iter <- mem_matrix

mem_matrix_trans_2_iter <- matrix(nrow = nrow(mem_matrix_2_iter), ncol = ncol(mem_matrix_2_iter))

for (i in 1 : nrow(mem_matrix_2_iter)) {
  for (j in 1 : ncol(mem_matrix)) {
    mem_matrix_trans_2_iter[i, j] = mem_matrix_2_iter[i, j]/max(mem_matrix_2_iter[,j])
  }
}

# Converting the matrix in a 1-0 matrix

mem_matrix_trans_2_iter[mem_matrix_trans_2_iter != 1] <- 0

# Converting the results to a vector, in order to compute clustering validity indexes

clustering <- numeric(ncol(mem_matrix_trans_2_iter))
for (j in 1 : length(clustering)) {
  clustering[j] <- which.max(mem_matrix_trans_2_iter[,j])
}

ground_truth <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth[j] <- P[[j]][1, 9]
}

library(dtwclust)
cvi(ground_truth, clustering) # The results get worst

# Results with hierarchical clustering 

load('hierar_clust.RData')
hierarch_clust <- cutree(clustering, 10)
cvi(ground_truth, hierarch_clust)

