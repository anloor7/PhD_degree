

# Lets perfom simulatiomns regarding a simple scenario in order to compare QAF, QC1 and QC2 and its competitors 

# Clusters

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()

theta_c1 <- matrix(c(0.3, 0.3, 0.3, 0.0, 0.3, 0.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0, 0.3), nrow = 4)
sigma_c1 <- 0.9*diag(4) + matrix(rep(0.1, 16), ncol = 4)
theta_c2 <- matrix(c(0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5), nrow = 4)
sigma_c2 <- 0.7*diag(4) + matrix(rep(0.3, 16), ncol = 4)
theta_c3 <- matrix(c(0.2, 0, 0.2, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0.2, 0, 0.2), nrow = 4)
sigma_c3 <- 0.5*diag(4) + matrix(rep(0.5, 16), ncol = 4)
theta_c4 <- matrix(c(0.1, 0, 0.1, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0, 0.1, 0, 0.1), nrow = 4)
sigma_c4 <- 0.3*diag(4) + matrix(rep(0.7, 16), ncol = 4)

ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)) # Ground truth 


# We are going to perfom 20 trials regarding each method, so we need 20 different lists cluster. We make a big list 
# containing this small lists 


cluster_trials <- list()
set.seed(1234)

for (j in 1 : 20){

  
  for (i in 1 : 5) {
    cluster1[[i]] <- VARMAsim(2000, malags = 1, theta = theta_c1, sigma = sigma_c1)$series
  }
  
  for (i in 1 : 5) {
    cluster2[[i]] <- VARMAsim(2000, malags = 1, theta = theta_c2, sigma = sigma_c1)$series
  }
  
  for (i in 1 : 5) {
    cluster3[[i]] <- VARMAsim(2000, malags = 1, theta = theta_c3,
                              sigma = sigma_c1)$series
  }
  
  for (i in 1 : 5) {
    cluster4[[i]] <- VARMAsim(2000, malags = 1, theta = theta_c4,
                              sigma = sigma_c1)$series
  }
  
cluster_trials[[j]] <- c(cluster1, cluster2, cluster3, cluster4)

}





# QAF WITH PAM 


boxplot_qaf_pam <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  gamma <- listTomatrix(lapply(cluster_trials[[i]], qaf_mts_coefs_xy_sep))
  clustering <- pam(gamma, 4)$cluster
  boxplot_qaf_pam[[i]] <- external_validation(ground_truth, clustering)
  
}

end_time <- Sys.time()
qaf_pam_time <- start_time - end_time
mean(boxplot_qaf_pam)
save(boxplot_qaf_pam, file = 'boxplot_qaf_pam.RData')



# QAF WITH K MEANS 

boxplot_qaf_km <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  gamma <- listTomatrix(lapply(cluster_trials[[i]], qaf_mts_coefs_xy_sep))
  auxiliary <- numeric()
  
  for (j in 1 : 3){
  
  clustering <- kmeans(gamma, 4)$cluster
  auxiliary[j] <- external_validation(ground_truth, clustering)
  
  }
  
  boxplot_qaf_km[i] <- max(auxiliary)
  
}
end_time <- Sys.time()
qaf_km_time <- (start_time - end_time)/3
mean(boxplot_qaf_km)
save(boxplot_qaf_km, file = 'boxplot_qaf_km.RData')




# QC1 WITH PAM 

boxplot_qc1_pam <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  coherence1 <- listTomatrix(lapply(cluster_trials[[i]], quantile_coherence))
  clustering <- pam(coherence1, 4)$cluster
  boxplot_qc1_pam[[i]] <- external_validation(ground_truth, clustering)
  
}
end_time <- Sys.time()
qc1_pam_time <- start_time - end_time
mean(boxplot_qc1_pam)
save(boxplot_qc1_pam, file = 'boxplot_qc1_pam.RData')




# QC1 WITH K MEANS

boxplot_qc1_km <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  coherence1 <- listTomatrix(lapply(cluster_trials[[i]], quantile_coherence))
  auxiliary <- numeric()
  
  for (j in 1 : 3){
    
    clustering <- kmeans(coherence1, 4)$cluster
    auxiliary[j] <- external_validation(ground_truth, clustering)
    
  }
  
  boxplot_qc1_km[i] <- max(auxiliary)
  
}
end_time <- Sys.time()
qc1_km_time <- (start_time - end_time)/3
mean(boxplot_qc1_km)
save(boxplot_qc1_km, file = 'boxplot_qc1_km.RData')






# QC2 WITH PAM 


boxplot_qc2_pam <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  coherence2 <- listTomatrix(lapply(cluster_trials[[i]], quantile_coherence_re_im))
  clustering <- pam(coherence2, 4)$cluster
  boxplot_qc2_pam[[i]] <- external_validation(ground_truth, clustering)
  
}
end_time <- Sys.time()
qc2_pam_time <- start_time - end_time
mean(boxplot_qc2_pam)
save(boxplot_qc2_pam, file = 'boxplot_qc2_pam.RData')




# QC2 WITH K MEANS

boxplot_qc2_km <- numeric()
start_time <- Sys.time()
for (i in 1 : 20) {
  
  coherence2 <- listTomatrix(lapply(cluster_trials[[i]], quantile_coherence_re_im))
  auxiliary <- numeric()
  
  for (j in 1 : 3){
    
    clustering <- kmeans(coherence2, 4)$cluster
    auxiliary[j] <- external_validation(ground_truth, clustering)
    
  }
  
  boxplot_qc2_km[i] <- max(auxiliary)
  
}
end_time <- Sys.time()
qc2_km_time <- (start_time - end_time)/3
mean(boxplot_qc2_km)
save(boxplot_qc2_km, file = 'boxplot_qc2_km.RData')




# DTW1 WITH PAM 


# Distance matrix

boxplot_dtw1_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){

dis_matrix <- matrix(0, 20, 20)

for (i in 1 : 20) {
  for (j in 1 : 20) {
    dis_matrix[i, j] <- dtw_mts(cluster_trials[[k]][[i]], cluster_trials[[k]][[j]])
  }
}

clustering <- pam(dis_matrix, 4)$cluster
boxplot_dtw1_pam[[k]] <- external_validation(ground_truth, clustering)

}
end_time <- Sys.time()
dtw1_pam_time <- start_time - end_time
mean(boxplot_dtw1_pam)
save(boxplot_dtw1_pam, file = 'boxplot_dtw1_pam.RData')




# DTW2 WITH PAM 


# Distance matrix

boxplot_dtw2_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  dis_matrix <- matrix(0, 20, 20)
  for (i in 1 : 20) {
    for (j in 1 : 20) {
      p <- ncol(cluster1[[1]])
      d_matrix <- matrix(0, p, p)
      for (s in 1 : p){
        for (v in 1 : p){
          
            d_matrix[s, v] <- EuclideanDistance(cluster_trials[[k]][[i]][,s], cluster_trials[[k]][[j]][,v]) 
          }
        }
      dis_matrix[i, j] <- dtw(d_matrix, distance.only = TRUE)$normalizedDistance
      # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'sakoechiba')$normalizedDistance
      # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'itakura')$normalizedDistance
    }
  }
  
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_dtw2_pam[[k]] <- external_validation(ground_truth, clustering)
  
}
end_time <- Sys.time()
dtw2_pam_time <- start_time - end_time
mean(boxplot_dtw2_pam)
save(boxplot_dtw2_pam, file = 'boxplot_dtw2_pam.RData')



# PDC WITH PAM 

boxplot_pdc_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  
  array_series <- array(dim = c(15000, 20, 4))
  
  for (i in 1 : 20) {
    array_series[,i,] <- cluster_trials[[k]][[i]]
  }
  
  dis_matrix <- pdcDist(array_series)
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_pdc_pam[[k]] <- external_validation(ground_truth, clustering)
  
}

end_time <- Sys.time()
pdc_pam_time <- start_time - end_time
mean(boxplot_pdc_pam)
save(boxplot_pdc_pam, file = 'boxplot_pdc_pam.RData')





# KST WITH PAM 

boxplot_kst_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  dis_matrix <- matrix(0, 20, 20)
  
  for (i in 1 : 20) {
    for (j in 1 : 20) {
      dis_matrix[i, j] <- j_divergence(cluster_trials[[k]][[i]], cluster_trials[[k]][[j]])
    }
  }
  
  diag(dis_matrix) <- 0 # Numerical error 
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_kst_pam[[k]] <- external_validation(ground_truth, clustering)
  
}

end_time <- Sys.time()
kst_pam_time <- start_time - end_time
mean(boxplot_kst_pam)
save(boxplot_kst_pam, file = 'boxplot_kst_pam_he.RData')





# KST WITH K MEANS 

boxplot_kst_km <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  auxiliary <- numeric()
  
  for (j in 1 : 3){
    
    clustering <- km_kst(cluster_trials[[k]], 4, dis = j_divergence_average)
    auxiliary[j] <- external_validation(ground_truth, clustering)
    
  }
  
  boxplot_kst_km[k] <- max(auxiliary)
  
}

end_time <- Sys.time()
kst_km_time <- (start_time - end_time)/3
mean(boxplot_kst_km)
save(boxplot_kst_km, file = 'boxplot_kst_km.RData')





# PCA WITH PAM 

boxplot_pca_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  clustercov <- list()
  for (i in 1 : length(cluster_trials[[k]])) {
    clustercov[[i]] <- cov(cluster_trials[[k]][[i]])
  }
  
  dis_matrix <- 1 - PCAsimilarity(clustercov)
  
  for (i in 1 : nrow(dis_matrix)) {
    for (j in i : ncol(dis_matrix)) {
      dis_matrix[i, j] <- dis_matrix[j, i]
    }
  }
  
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_pca_pam[[k]] <- external_validation(ground_truth, clustering)
  
}

end_time <- Sys.time()
pca_pam_time <- start_time - end_time
mean(boxplot_pca_pam)
save(boxplot_pca_pam, file = 'boxplot_pca_pam.RData')





# PCA WITH K MEANS 

boxplot_pca_km <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  auxiliary <- numeric()
  
  for (j in 1 : 3){
    
    clustering <- km_mts_ss(cluster_trials[[k]], 4)
    auxiliary[j] <- external_validation(ground_truth, clustering)
    
  }
  
  boxplot_pca_km[k] <- max(auxiliary)
  
}

end_time <- Sys.time()
pca_km_time <- (start_time - end_time)/3
mean(boxplot_pca_km)
save(boxplot_pca_km, file = 'boxplot_pca_km.RData')





# MC2PCA with K MEANS 

boxplot_mc2pca_km <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  auxiliary <- numeric()
  
  for (j in 1 : 3){
    
    clustering <- mc2pca(cluster_trials[[k]], 4, lambda = 0.95, niter = 1000, tol = 0.01)
    auxiliary[j] <- external_validation(ground_truth, clustering)
    
  }
  
  boxplot_mc2pca_km[k] <- max(auxiliary)
  
}

end_time <- Sys.time()
mc2pca_km_time <- (start_time - end_time)/3
mean(boxplot_mc2pca_km)
save(boxplot_mc2pca_km, file = 'boxplot_mc2pca_km.RData')


# WAVELETS WITH PAM 

# Distance matrix

boxplot_wavelets_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  
  J <- 6      # number of scales (see Table 3, page 45, in D'urso and Maharaj 2012)
  wf <- "d4"
  features <- lapply(cluster_trials[[k]], wavelet.features, wf = wf, J = J) 
  dis_matrix <- proxy::dist(features, wave_dist)  
  
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_wavelets_pam[[k]] <- external_validation(ground_truth, clustering)
  
}


end_time <- Sys.time()
wavelets_pam_time <- start_time - end_time
mean(boxplot_wavelets_pam)
save(boxplot_wavelets_pam, file = 'boxplot_wavelets_pam.RData')


# GCC WITH PAM 

# Distance matrix


boxplot_gcc_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  dis_matrix <- matrix(0, 20, 20)
  
  for (i in 1 : 20) {
    for (j in 1 : 20) {
      dis_matrix[i, j] <- EuclideanDistance(gcc_features_mts(cluster_trials[[k]][[i]]), 
                                            gcc_features_mts(cluster_trials[[k]][[j]]))
    }
  }
  
  diag(dis_matrix) <- 0 # Numerical error 
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_gcc_pam[[k]] <- external_validation(ground_truth, clustering)
  
}



end_time <- Sys.time()
gcc_pam_time <- start_time - end_time
mean(boxplot_gcc_pam)
save(boxplot_gcc_pam, file = 'boxplot_gcc_pam.RData')





# EROS WITH PAM 


boxplot_eros_pam <- numeric()
start_time <- Sys.time()

for (k in 1 : 20){
  
  
  PCAs  <- lapply(cluster_trials[[k]], prcomp)
  SIGMA <- (do.call(cbind, lapply(PCAs, '[[',1)))^2 
  V <- lapply(PCAs, '[[',2)                         
  w <- w.eros(SIGMA)                               
  dis_matrix <- proxy::dist(V, V, Eros, w)                     
  dis_matrix <- as.dist(dis_matrix)     
  
  diag(dis_matrix) <- 0 # Numerical error 
  clustering <- pam(dis_matrix, 4)$cluster
  boxplot_eros_pam[[k]] <- external_validation(ground_truth, clustering)
  
}



end_time <- Sys.time()
eros_pam_time <- start_time - end_time
mean(boxplot_eros_pam)
save(boxplot_eros_pam, file = 'boxplot_eros_pam.RData')
