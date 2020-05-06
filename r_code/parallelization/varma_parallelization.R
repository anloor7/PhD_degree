

# We are going to make parallel compuation in order to fit a VAR(1) model to 1000 simulated MTS

# First, we simulate the multivariate time series

varmalist <- list()

phi <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)
sigma = diag(2)
for (i in 1 : 1000) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, sigma = sigma, phi = phi)$series
}

# Creating a function to fit a VARMA model to one matrix

varma_coeff <- function(a){
   VARMA(a, p = 1)$coef
}

library(parallel)
n <- 7 # Number of cores 

c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('varma_coeff', 'varmalist'))


time <- system.time(v <- parLapply(c1, varmalist, varma_coeff))
time # 38.3 seconds
stopCluster(c1)

# Time needed without parallelization

time_without <- system.time(for(i in 1 : 1000) {VARMA(varmalist[[i]], p = 1)$coefs})




