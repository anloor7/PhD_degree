

# Lets check the performance of our approach with the Character Trajectories Dataset 

library(dtwclust)
data(uciCT)
ct <- CharTrajMV # List containing 100, multivariate time series, each one with 3 variables

ground_truth <- numeric() # Ground truth 


for (i in  1 : 20) {
  ground_truth <- c(ground_truth, rep(i, 5))
}


# Removing the rows with all-zero values from each MTS

for (i in 1 : 100) {
  ct[[i]] <- ct[[i]][as.logical(rowSums(ct[[i]] != 0)), ]
  colnames(ct[[i]]) <- NULL
}

sigma <- list()

for (i in 1 : 100) {
  sigma[[i]] <- cov(ct[[i]])
}

s <- cpca(sigma, lambda = 0.60)
ct_reduced <- list()

for (i in 1 : 100) {
  ct_reduced[[i]] <- ct[[i]] %*% s
}

coefs <- varma_coefs(ct_reduced)


d <- listTomatrix(coefs)
clustering <- kmeans(d, 20)$cluster
external_validation(ground_truth, clustering, method = 'purity')

external_validation(fuzzytocrisp(t(fcm(d, 20)$u)), ground_truth)


# Now, we are going to perform the analysis regarding the whole dataset Character Trajectories

library(R.matlab) # The data is in Matlab format 


ct_whole_read <- readMat('ct.mat')
meta <- ct_whole_read[[1]]
ct_whole <- ct_whole_read[[2]]

for (i in 1 : 2858) {
  ct_whole[[i]] <- t(ct_whole[[i]][[1]]) # Converting the MTS into matrixes
}

for (i in 1 : 2858) {
  ct_whole[[i]] <- ct_whole[[i]][as.logical(rowSums(ct_whole[[i]] != 0)), ] # Removing all-zero rows 
  colnames(ct_whole[[i]]) <- NULL
}



sigma <- list()

for (i in 1 : 2858) {
  sigma[[i]] <- cov(ct_whole[[i]])
}

s <- cpca(sigma, lambda = 0.70)
ct_whole_reduced <- list()

for (i in 1 : 2858) {
  ct_whole_reduced[[i]] <- ct_whole[[i]] %*% s
}


# Computing the clustering solution via parallel computing 

# Creating a function to fit a VARMA model to one matrix

varma_coeff <- function(a){
  a <- VARMA(a, p = 1)
 c(as.vector(a$coef), a$aic, a$bic)
}

n <- 7 # Number of cores 

c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('varma_coeff', 'ct_whole_reduced'))

time <- system.time(coefs <- parLapply(c1, ct_whole_reduced, varma_coeff))
time # 

stopCluster(c1)


d <- listTomatrix(coefs)
clustering <- kmeans(d, 20)$cluster
table(clustering)
external_validation(ground_truth, clustering, method = 'adjusted_rand_index')

