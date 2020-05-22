

# Lets perform some analyses regarding distance between varma. The first is going to be based on simulations, and the second, on the 
# Libras dataset 

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/functions/kmeans_function.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/functions/varma_coef.R')
source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/functions/listTomatrix.R')

# Simulation of 20 VARMA(1, 1) processes stemming from two different model generators

X <- list()

theta1 <- matrix(c(0.5, 0.6, 0.6, 0.5), nrow = 2)
phi1 <- matrix(c(0.2, 0, 0.3, 0.7), nrow = 2)
epsilon1 <- diag(2)

for (i in 1 : 10){
X[[i]] <- VARMAsim(150, arlags = c(1), malags = c(1), phi = phi1, theta = theta1, sigma = epsilon1)$series}


theta2 <- matrix(c(1, 0.4, 0.4, 0.8), nrow = 2)
phi2 <- matrix(c(0.5, 0.2, 0.8, 0.7), nrow = 2)
epsilon2 <- diag(2)

for (i in 11 : 20){
  X[[i]] <- VARMAsim(150, arlags = c(1), malags = c(1), phi = phi2, theta = theta2, sigma = epsilon2)$series}

ground_truth <- c(rep(1, 10), rep(2, 10))

# To perform clustering, we can work with the corresponding vector of coefficients and the Euclidean Distance

v_c <- list()
vc <- varma_coefs(X)

clustering <- km_mts(vc, K = 2, niter = 100, tol = 0.01, dis = EuclideanDistance)
external_validation(ground_truth, clustering, summary_stats = T)


# Dataset Libras 

# Loading the data 

libras <- load('libras.RData')
M <- list(length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- t(S[[i]])[1:45,]
}

ground_truth_libras <- numeric(length(Y))

for (j in 1 : length(Y)) {
  ground_truth_libras[j] <- S[[j]][1, 46]
}

time <- system.time(vc <- varma_coefs(M))
clustering <- km_mts(vc, K = 15, niter = 1000, tol = 0.01, dis = EuclideanDistance)
external_validation(ground_truth_libras, clustering, method = 'jaccard_index')


# Removing outliers

index <- numeric(360)

for (i in 1 : 360) {
  if (sum(vc[[i]][9:12] > 10) > 0)
    {index <- c(index, i)}
}

index_new <- index[which(index != 0)]
vc_new <- vc[-index_new]
ground_truth_libras_new <- ground_truth_libras[-index_new]


clustering <- km_mts(vc_new, K = 15, niter = 100, tol = 0.01, dis = EuclideanDistance)
external_validation(ground_truth_libras_new, clustering, method = 'purity')

vc_matrix <- listTomatrix(vc_new)

B <- numeric(1000)
c <- list()
for (i in 1 : 1000) {
  # c[[i]] <- kmeans(vc_matrix, 15)
  # B[i] <- intCriteria(vc_matrix, c[[i]]$cluster, 'Dunn')
  # B[i] <- c[[i]]$tot.withinss
  B[i] <- external_validation(kmeans(vc_matrix, 15)$cluster, ground_truth_libras_new, method = 'purity')
}

mean(B)

B <- numeric(20)
c <- list()
for (i in 1 : 20) {
  # c[[i]] <- kmeans(vc_matrix, 15)
  # B[i] <- intCriteria(vc_matrix, c[[i]]$cluster, 'Dunn')
  # B[i] <- c[[i]]$tot.withinss
  B[i] <- external_validation(fuzzytocrisp(t(fcm(vc_matrix, 15)$u)), ground_truth_libras_new, method = 'rand_index')
}


mean(B)

index <- which.min(B)
clustering <- c[[index]]$cluster
external_validation(ground_truth_libras_new, clustering, method = 'purity')

# Simulations of the second experiment of Baragona 

source('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/git_hub/PhD_degree/r_code/varma/simulations_baragona.R')

time <- system.time(vc <- varma_coefs(experiment2))

# K-means

b <- numeric()

for (i in 1 : 100) {
clustering <- km_mts(vc, K = 3, niter = 700, tol = 0.01, dis = EuclideanDistance)
b[i] <- external_validation(ground_truth2, clustering, method = "jaccard_index")
}

mean(b)

# Fuzzy C-means

b <- numeric()

for (i in 1 : 100) {
  clustering <- fuzzytocrisp(fcm_mts(vc, K = 3, niter = 700, tol = 0.01, dis = EuclideanDistance))
  b[i] <- external_validation(ground_truth2, clustering, method = "jaccard_index")
}

mean(b)

# Internal validity indexes 

# To compute internal indexes, the input must be a kmeans (stats) object 

matrixvc <- matrix(nrow = length(vc), ncol = length(vc[[1]]))

for (i in 1 : length(vc)) {
  matrixvc[i,] <- vc[[i]] 
}

clustering <- kmeans(matrixvc, c = 15)
external_validation(ground_truth2, clustering$cluster, method = "jaccard_index")
intCriteria(matrixvc, clustering$cluster, 'Xie_Beni')

