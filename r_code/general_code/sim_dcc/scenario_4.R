





# Lets perform some simulations regarding DCC models (Engle 2002)

# Situation: 4 escenarios with 100 MTS  each one, each MTS containing 4 variables. In each scenario, 
# the series are generated from a DCC(2, 2) model with different parameters. The mean vector follows, separately for
# each component, a ARMA(0, 0) model. The EGarch model
# for each individual volatillity of univariate time series is different in each one of the clusters


# Cluster 1


vectoruspec <- list(omega = 0.5, alpha1 = 0.5, beta1 = 0.5, gamma1 = 0.5) # Parameters of eGARCH model 
uspec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm", mean.model = list(armaOrder=c(0, 0), include.mean = F), 
                    fixed.pars = vectoruspec)



model <- multispec(replicate(4,uspec)) # Replication of eGarch model for the four variables 

vectormvspec <- list(dcca1 = 0.2, dccb1 = 0.2, dcca2 = 0.1, dccb2 = 0.1 ) # Parameters of DCC model 
mvspec <- dccspec(model, dccOrder = c(2, 2), model = "DCC", distribution = # Defining DCC model 
                    "mvnorm", fixed.pars = vectormvspec)

set.seed(1234)
cluster1 <- list()

for (i in 1 : 100) {
  
  sim <- dccsim(mvspec, n.sim = 100, preQ = diag(4), Qbar = diag(4), Nbar = diag(4))
  cluster1[[i]] <- slot(sim, 'msim')[[4]][[1]] # Simulated series 
  
}


# Cluster 2 


vectoruspec <- list(omega = 0.2, alpha1 = 0.6, beta1 = 0.7, gamma1 = 0.9) # Parameters of eGARCH model 
uspec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm", mean.model = list(armaOrder=c(0, 0), include.mean = F), 
                    fixed.pars = vectoruspec)



model <- multispec(replicate(4,uspec)) # Replication of eGarch model for the four variables 
vectormvspec <- list(dcca1 = 0.1, dccb1 = 0.5, dcca2 = 0.2, dccb2 = 0) # Parameters of DCC model 
mvspec <- dccspec(model, dccOrder = c(2, 2), model = "DCC", distribution = # Defining DCC model 
                    "mvnorm", fixed.pars = vectormvspec)
cluster2 <- list()

for (i in 1 : 100) {
  
  sim <- dccsim(mvspec, n.sim = 100, preQ = diag(4), Qbar = diag(4), Nbar = diag(4))
  cluster2[[i]] <- slot(sim, 'msim')[[4]][[1]] # Simulated series 
  
}


# Cluster 3

vectoruspec <- list(omega = 0.5, alpha1 = 0.1, beta1 = 0.9, gamma1 = 0.5) # Parameters of eGARCH model 
uspec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm", mean.model = list(armaOrder=c(0, 0), include.mean = F), 
                    fixed.pars = vectoruspec)



model <- multispec(replicate(4,uspec)) # Replication of eGarch model for the four variables 

vectormvspec <- list(dcca1 = 0.1, dccb1 = 0.2, dcca2 = 0.3, dccb2 = 0.3) # Parameters of DCC model 
mvspec <- dccspec(model, dccOrder = c(2, 2), model = "DCC", distribution = # Defining DCC model 
                    "mvnorm", fixed.pars = vectormvspec)
cluster3 <- list()

for (i in 1 : 100) {
  
  sim <- dccsim(mvspec, n.sim = 100, preQ = diag(4), Qbar = diag(4), Nbar = diag(4))
  cluster3[[i]] <- slot(sim, 'msim')[[4]][[1]] # Simulated series 
  
}


# Cluster 4

vectoruspec <- list(omega = 0.9, alpha1 = 0.9, beta1 = 0.45, gamma1 = 0.3) # Parameters of eGARCH model 
uspec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm", mean.model = list(armaOrder=c(0, 0), include.mean = F), 
                    fixed.pars = vectoruspec)



model <- multispec(replicate(4,uspec)) # Replication of eGarch model for the four variables 

vectormvspec <- list(dcca1 = 0.25, dccb1 = 0.1, dcca2 = 0.32, dccb2 = 0.32) # Parameters of DCC model 
mvspec <- dccspec(model, dccOrder = c(2, 2), model = "DCC", distribution = # Defining DCC model 
                    "mvnorm", fixed.pars = vectormvspec)
cluster4 <- list()

for (i in 1 : 100) {
  
  sim <- dccsim(mvspec, n.sim = 100, preQ = diag(4), Qbar = diag(4), Nbar = diag(4))
  cluster4[[i]] <- slot(sim, 'msim')[[4]][[1]] # Simulated series 
  
}



cluster <- c(cluster1, cluster2, cluster3, cluster4)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))







# Parallelization 

n <- 7 # Number of cores 

# c1 <- makeCluster(n) # Making a cluster object
clusterExport(c1, c('cluster', 'qaf_mts_coefs_xy_sep', 'quantile_coherence', 'quantile_coherence_re_im', 'kmeans_mc_av_ari_scenario2'))


#  Performing k-means QAF 

gamma <- listTomatrix(lapply(cluster, qaf_mts_coefs_xy_sep))
gammal <- list(length = 100) 

for (i in 1 : 100) {
  gammal[[i]] <- gamma
}


b <- numeric()
b <- parLapply(c1, gammal, kmeans_mc_av_ari_scenario2)
mean(listTomatrix(b))

# Performing k-means quantile coherence (modulus)

coherence <- listTomatrix(lapply(cluster, quantile_coherence))
coherencel <- list(length = 100) 

for (i in 1 : 100) {
  coherencel[[i]] <- coherence
}

b <- numeric()
b <- parLapply(c1, coherencel, kmeans_mc_av_ari_scenario2)
mean(listTomatrix(b))

# Performing k-means quantile coherence (real-imaginary)

coherence <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
coherencel <- list(length = 100) 

for (i in 1 : 100) {
  coherencel[[i]] <- coherence
}

b <- numeric()
b <- parLapply(c1, coherencel, kmeans_mc_av_ari_scenario2)
mean(listTomatrix(b))

clustering <- kmeans(coherence, 3)$cluster
external_validation(ground_truth, clustering)
