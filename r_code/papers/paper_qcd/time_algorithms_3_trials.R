

# Lets compute the time the 6 metrics spend in 1 trial of the algoritm (250 and 1000)

# Loading the data

load('cluster_1_250_n.RData')
load('cluster_1_500_n.RData')
load('cluster_1_1000_n.RData')

cluster_250 <- loadRData(paste0('cluster_1_', 250, '_n.RData'))[[1]] # Taking the first 25 series
cluster_1000 <- loadRData(paste0('cluster_1_', 1000, '_n.RData'))[[1]] # Taking the first 25 series



# QCD


set.seed(1234)

time_general_qcd_250 <- numeric()
time_general_qcd_1000 <- numeric()


for (i in 1 : 30) {
  
start_time <- Sys.time()
coherence2 <- listTomatrix(lapply(cluster_250, quantile_quantities_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance) 
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_qcd_250[i] <- end_time - start_time

}

for (i in 1 : 30) {

start_time <- Sys.time()
coherence2 <- listTomatrix(lapply(cluster_1000, quantile_quantities_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance) 
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_qcd_1000[i] <- end_time - start_time

}

save(time_general_qcd_250, file = 'time_general_qcd_250.RData')
save(time_general_qcd_1000, file = 'time_general_qcd_1000.RData')




# DM

set.seed(1234)

time_general_w_250 <- numeric()
time_general_w_1000 <- numeric()


J <- 5
w <- 'd4'
for (i in 1 : 100) {

start_time <- Sys.time()
features <- lapply(cluster_250, wavelet_features, wf = w, J = J)
dis_matrix <- proxy::dist(features, wave_dist)  
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_w_250[[i]] <- end_time - start_time

}



J <- 7
w <- 'd4'
for (i in 1 : 100) {
  
start_time <- Sys.time()
features <- lapply(cluster_1000, wavelet_features, wf = w, J = J)
dis_matrix <- proxy::dist(features, wave_dist)  
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_w_1000[[i]] <- end_time - start_time

}

save(time_general_w_250, file = 'time_general_w_250.RData')
save(time_general_w_1000, file = 'time_general_w_1000.RData')


# AP 



set.seed(1234)
time_general_gcc_250 <- numeric()
time_general_gcc_1000 <- numeric()

for (i in 1 : 5) {
  
start_time <- Sys.time()
features <- listTomatrix(lapply(cluster_250, gcc_features_reduced_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_gcc_250[[i]] <- end_time - start_time

}

for (i in 1 : 5) {
  
start_time <- Sys.time()
features <- listTomatrix(lapply(cluster_1000, gcc_features_reduced_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_gcc_1000[[i]] <- end_time - start_time

}

save(time_general_gcc_250, file = 'time_general_gcc_250.RData')
save(time_general_gcc_1000, file = 'time_general_gcc_1000.RData')



# KST


set.seed(1234)

time_general_kst_250 <- numeric()
time_general_kst_1000 <- numeric()
dis_matrix <- matrix(0, 25, 25)

for (i in 1 : 3) {
  
start_time <- Sys.time()
for (j1 in 1 : 25) {
  for (j2 in 1 : 25) {
    dis_matrix[j1, j2] <- j_divergence(cluster_250[[j1]], cluster_250[[j2]])
  }
}

print(i)
diag(dis_matrix) <- 0 # Numerical error 
pamc <- pam(dis_matrix, 5)
end_time <- Sys.time()
time_general_kst_250[[i]] <- end_time - start_time

}

for (i in 1 : 3) {
  
  start_time <- Sys.time()
  for (j1 in 1 : 25) {
    for (j2 in 1 : 25) {
      dis_matrix[j1, j2] <- j_divergence(cluster_1000[[j1]], cluster_1000[[j2]])
    }
  }
  
  print(i)
  diag(dis_matrix) <- 0 # Numerical error 
  pamc <- pam(dis_matrix, 5)
  end_time <- Sys.time()
  time_general_kst_1000[[i]] <- end_time - start_time
  
}

save(time_general_kst_250, file = 'time_general_kst_250.RData')
save(time_general_kst_1000, file = 'time_general_kst_1000.RData')



# MAHARAJ


set.seed(1234)

time_general_mah_250 <- numeric()
time_general_mah_1000 <- numeric()


for (i in 1 : 3) {
  
  start_time <- Sys.time()
  features <- maharaj(cluster_250)
  dis_matrix <- proxy::dist(features, EuclideanDistance) 
  pamc <- pam(dis_matrix, 5)
  end_time <- Sys.time()
  time_general_mah_250[i] <- end_time - start_time
  print(i)
  
}

for (i in 1 : 3) {
  
  start_time <- Sys.time()
  features <- maharaj(cluster_1000)
  dis_matrix <- proxy::dist(features, EuclideanDistance) 
  pamc <- pam(dis_matrix, 5)
  end_time <- Sys.time()
  time_general_mah_1000[i] <- end_time - start_time
  print(i)
  
}

save(time_general_mah_250, file = 'time_general_mah_250.RData')
save(time_general_mah_1000, file = 'time_general_mah_1000.RData')



# PCA

set.seed(1234)

time_general_pca_250 <- numeric()
time_general_pca_1000 <- numeric()
dis_matrix <- matrix(0, 25, 25)

for (i in 1 : 50) {
  
start_time <- Sys.time()
clustercov <- list()

for (m in 1 : length(cluster_250)) {
  clustercov[[m]] <- cov(cluster_250[[m]])
}

dis_matrix <- 1 - PCAsimilarity(clustercov)
dis_matrix[col(dis_matrix) >= row(dis_matrix)] <- 0
dis_matrix <- as.dist(dis_matrix)
pamc <- pam(dis_matrix, 5)$cluster
end_time <- Sys.time()
time_general_pca_250[i] <- end_time - start_time


}



for (i in 1 : 50) {
  
  start_time <- Sys.time()
  clustercov <- list()
  
  for (m in 1 : length(cluster_250)) {
    clustercov[[m]] <- cov(cluster_1000[[m]])
  }
  
  dis_matrix <- 1 - PCAsimilarity(clustercov)
  dis_matrix[col(dis_matrix) >= row(dis_matrix)] <- 0
  dis_matrix <- as.dist(dis_matrix)
  pamc <- pam(dis_matrix, 5)$cluster
  end_time <- Sys.time()
  time_general_pca_1000[i] <- end_time - start_time
  
  
}


