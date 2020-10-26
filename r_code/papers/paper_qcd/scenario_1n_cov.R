

n <- 5  # Number of series per cluster
K <- 5 # Number of clusters
B <- 100 # Number of Monte Carlos trials
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n), rep(5, n)) # Ground Truth 


cluster <-  list() 
cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()
cluster5 <- list()
epsilon <-  list()

set.seed(1234)
for (l in c(250, 500, 1000)){
  
  
  
  for (k in 1 : B){
    
    phi_c1 <- matrix(c(0.6, -0.4, 0, 0.5, 0.5, -0.5, 0, 0.3, 0.7), nrow = 3)
    sigma_c1 <- matrix(c(1, 0.5, 0.75, 0.5, 1, 0.5, 0.75, 0.5, 1), nrow = 3)
    phi_c2 <- matrix(c(0.4, -0.4, 0, 0.4, 0.5, -0.5, 0, 0.4, 0.7), nrow = 3)
    sigma_c2 <- sigma_c1
    theta_c3 <- matrix(c(0.6, -0.4, 0, 0.5, 0.5, -0.5, 0, 0.3, 0.7), nrow = 3)
    sigma_c3 <- sigma_c1
    theta_c4 <- matrix(c(0.4, -0.4, 0, 0.4, 0.5, -0.5, 0, 0.4, 0.7), nrow = 3)
    sigma_c4 <- sigma_c1
    phi_c5 <- phi_c1
    theta_c5 <- theta_c3
    sigma_c5 <- sigma_c1
    sigma <- list()
    
    for (i in 1 : n) {
      
      
      epsilon[[i]] <- rmvnorm(l,  sigma = sigma_c1) # Normal
      cluster1[[i]] <- varma(l, k = 3, VAR = phi_c1, innov = epsilon[[i]])
      
    }
    
    for (i in 1 : n) {
      
      sigma[[i]] <- diag(3)
      epsilon[[i]] <- rmvnorm(l,  sigma = sigma_c2) # Normal
      cluster2[[i]] <- varma(l, k = 3, VAR = phi_c2, innov = epsilon[[i]])
      
    }
    
    for (i in 1 : n) {
      
      sigma[[i]] <- diag(3)
      epsilon[[i]] <- rmvnorm(l,  sigma = sigma_c3) # Normal
      cluster3[[i]] <- varma(l, k = 3, VMA = theta_c3, innov = epsilon[[i]])
      
    }
    
    for (i in 1 : n) {
      
      
      sigma[[i]] <- diag(3)
      epsilon[[i]] <- rmvnorm(l,  sigma = sigma_c4) # Normal
      cluster4[[i]] <- varma(l, k = 3, VMA = theta_c4, innov = epsilon[[i]])
      
    }
    
    
    for (i in 1 : n) {
      
      sigma[[i]] <- diag(3)
      epsilon[[i]] <- rmvnorm(l,  sigma = sigma_c5) # Normal
      cluster5[[i]] <- varma(l, k = 3, VAR = phi_c5, VMA = theta_c5, innov = epsilon[[i]])
      
    }
    
    
    
    cluster[[k]] <- c(cluster1, cluster2, cluster3, cluster4, cluster5)
    
  }
  
  
  # Saving the simulations 
  
  
  if (l == 250){
    
    cluster_1_250_n_cov <- cluster 
    save(cluster_1_250_n_cov, file = 'cluster_1_250_n_cov.RData')
    
  }
  
  if (l == 500){
    
    cluster_1_500_n_cov <- cluster 
    save(cluster_1_500_n_cov, file = 'cluster_1_500_n_cov.RData')
    
  }
  
  if (l == 1000){
    
    cluster_1_1000_n_cov <- cluster 
    save(cluster_1_1000_n_cov, file = 'cluster_1_1000_n_cov.RData')
    
  }
  
}


# Loading the data and making a cluster object

n <- 5  # Number of series per cluster
K <- 5 # Number of clusters
B <- 100 # Number of Monte Carlos trials
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n), rep(5, n)) # Ground Truth
c1 <- makeCluster(7) # Making a cluster object
registerDoParallel(cores = 7)
load('cluster_1_250_n_cov.RData')
load('cluster_1_500_n_cov.RData')
load('cluster_1_1000_n_cov.RData')



# QC2

ari_qc2_1_250_n_cov <- numeric()
jaccard_qc2_1_250_n_cov <- numeric()
larsen_qc2_1_250_n_cov <- numeric()
loo_qc2_1_250_n_cov <- numeric()
ari_qc2_1_500_n_cov <- numeric()
jaccard_qc2_1_500_n_cov <- numeric()
larsen_qc2_1_500_n_cov <- numeric()
loo_qc2_1_500_n_cov <- numeric()
ari_qc2_1_1000_n_cov <- numeric()
jaccard_qc2_1_1000_n_cov <- numeric()
larsen_qc2_1_1000_n_cov <- numeric()
loo_qc2_1_1000_n_cov <- numeric()
clustering_qc2_1_250_n_cov <- list()
clustering_qc2_1_500_n_cov <- list()
clustering_qc2_1_1000_n_cov <- list()
time_qc2_1_250_n_cov <- numeric()
time_qc2_1_500_n_cov <- numeric()
time_qc2_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  clusterExport(c1, c('cluster', 'quantile_quantities_re_im'))
  start_time <- Sys.time()
  
  
  for (i in (1:B)) {
    
    coherence2 <- parLapply(c1, cluster[[i]], quantile_quantities_re_im)
    dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
    clustering <- pam(dis_matrix, K)$cluster
    
    if (l == 250) {
      
      ari_qc2_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_qc2_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_qc2_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_qc2_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_qc2_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_qc2_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_qc2_1_250_n_cov, file = 'ari_qc2_1_250_n_cov.RData')
save(ari_qc2_1_500_n_cov, file = 'ari_qc2_1_500_n_cov.RData')
save(ari_qc2_1_1000_n_cov, file = 'ari_qc2_1_1000_n_cov.RData')
save(jaccard_qc2_1_250_n_cov, file = 'jaccard_qc2_1_250_n_cov.RData')
save(jaccard_qc2_1_500_n_cov, file = 'jaccard_qc2_1_500_n_cov.RData')
save(jaccard_qc2_1_1000_n_cov, file = 'jaccard_qc2_1_1000_n_cov.RData')
save(larsen_qc2_1_250_n_cov, file = 'larsen_qc2_1_250_n_cov.RData')
save(larsen_qc2_1_500_n_cov, file = 'larsen_qc2_1_500_n_cov.RData')
save(larsen_qc2_1_1000_n_cov, file = 'larsen_qc2_1_1000_n_cov.RData')
save(loo_qc2_1_250_n_cov, file = 'loo_qc2_1_250_n_cov.RData')
save(loo_qc2_1_500_n_cov, file = 'loo_qc2_1_500_n_cov.RData')
save(loo_qc2_1_1000_n_cov, file = 'loo_qc2_1_1000_n_cov.RData')
save(clustering_qc2_1_250_n_cov, file = 'clustering_qc2_1_250_n_cov.RData')
save(clustering_qc2_1_500_n_cov, file = 'clustering_qc2_1_500_n_cov.RData')
save(clustering_qc2_1_1000_n_cov, file = 'clustering_qc2_1_1000_n_cov.RData')
save(time_qc2_1_250_n_cov, file = 'time_qc2_1_250_n_cov.RData')
save(time_qc2_1_500_n_cov, file = 'time_qc2_1_500_n_cov.RData')
save(time_qc2_1_1000_n_cov, file = 'time_qc2_1_1000_n_cov.RData')










# DM



ari_w_1_250_n_cov <- numeric()
jaccard_w_1_250_n_cov <- numeric()
larsen_w_1_250_n_cov <- numeric()
loo_w_1_250_n_cov <- numeric()
ari_w_1_500_n_cov <- numeric()
jaccard_w_1_500_n_cov <- numeric()
larsen_w_1_500_n_cov <- numeric()
loo_w_1_500_n_cov <- numeric()
ari_w_1_1000_n_cov <- numeric()
jaccard_w_1_1000_n_cov <- numeric()
larsen_w_1_1000_n_cov <- numeric()
loo_w_1_1000_n_cov <- numeric()
clustering_w_1_250_n_cov <- list()
clustering_w_1_500_n_cov <- list()
clustering_w_1_1000_n_cov <- list()
time_w_1_250_n_cov <- numeric()
time_w_1_500_n_cov <- numeric()
time_w_1_1000_n_cov <- numeric()


wf <- 'd4'


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  clusterExport(c1, c('cluster', 'wavelet_features'))
  
  if (l == 250) {
    
    J <- 5
    
  }
  
  
  if (l == 500) {
    
    J <- 6
    
  }
  
  if (l == 1000) {
    
    J <- 7
    
  }
  
  for (i in  (1 : B)) {
    
    features <- parLapply(c1, cluster[[i]], wavelet_features, J = J, wf = wf)
    dis_matrix <- proxy::dist(features, wave_dist)  
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_w_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_w_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_w_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_w_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_w_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_w_1_1000_n_cov <- start_time - end_time
    
  }
  
  
}


save(ari_w_1_250_n_cov, file = 'ari_w_1_250_n_cov.RData')
save(ari_w_1_500_n_cov, file = 'ari_w_1_500_n_cov.RData')
save(ari_w_1_1000_n_cov, file = 'ari_w_1_1000_n_cov.RData')
save(jaccard_w_1_250_n_cov, file = 'jaccard_w_1_250_n_cov.RData')
save(jaccard_w_1_500_n_cov, file = 'jaccard_w_1_500_n_cov.RData')
save(jaccard_w_1_1000_n_cov, file = 'jaccard_w_1_1000_n_cov.RData')
save(larsen_w_1_250_n_cov, file = 'larsen_w_1_250_n_cov.RData')
save(larsen_w_1_500_n_cov, file = 'larsen_w_1_500_n_cov.RData')
save(larsen_w_1_1000_n_cov, file = 'larsen_w_1_1000_n_cov.RData')
save(loo_w_1_250_n_cov, file = 'loo_w_1_250_n_cov.RData')
save(loo_w_1_500_n_cov, file = 'loo_w_1_500_n_cov.RData')
save(loo_w_1_1000_n_cov, file = 'loo_w_1_1000_n_cov.RData')
save(clustering_w_1_250_n_cov, file = 'clustering_w_1_250_n_cov.RData')
save(clustering_w_1_500_n_cov, file = 'clustering_w_1_500_n_cov.RData')
save(clustering_w_1_1000_n_cov, file = 'clustering_w_1_1000_n_cov.RData')
save(time_w_1_250_n_cov, file = 'time_w_1_250_n_cov.RData')
save(time_w_1_500_n_cov, file = 'time_w_1_500_n_cov.RData')
save(time_w_1_1000_n_cov, file = 'time_w_1_1000_n_cov.RData')










# AP

ari_gcc_1_250_n_cov <- numeric()
jaccard_gcc_1_250_n_cov <- numeric()
larsen_gcc_1_250_n_cov <- numeric()
loo_gcc_1_250_n_cov <- numeric()
ari_gcc_1_500_n_cov <- numeric()
jaccard_gcc_1_500_n_cov <- numeric()
larsen_gcc_1_500_n_cov <- numeric()
loo_gcc_1_500_n_cov <- numeric()
ari_gcc_1_1000_n_cov <- numeric()
jaccard_gcc_1_1000_n_cov <- numeric()
larsen_gcc_1_1000_n_cov <- numeric()
loo_gcc_1_1000_n_cov <- numeric()
clustering_gcc_1_250_n_cov <- list()
clustering_gcc_1_500_n_cov <- list()
clustering_gcc_1_1000_n_cov <- list()
time_gcc_1_250_n_cov <- numeric()
time_gcc_1_500_n_cov <- numeric()
time_gcc_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  clusterExport(c1, c('cluster', 'gcc_features_mts'))
  
  for (i in  (1 : B)) {
    
    features <- parLapply(c1, cluster[[i]], gcc_features_mts)
    dis_matrix <- proxy::dist(features, EuclideanDistance)
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_gcc_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_gcc_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_gcc_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_gcc_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_gcc_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_gcc_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_gcc_1_250_n_cov, file = 'ari_gcc_1_250_n_cov.RData')
save(ari_gcc_1_500_n_cov, file = 'ari_gcc_1_500_n_cov.RData')
save(ari_gcc_1_1000_n_cov, file = 'ari_gcc_1_1000_n_cov.RData')
save(jaccard_gcc_1_250_n_cov, file = 'jaccard_gcc_1_250_n_cov.RData')
save(jaccard_gcc_1_500_n_cov, file = 'jaccard_gcc_1_500_n_cov.RData')
save(jaccard_gcc_1_1000_n_cov, file = 'jaccard_gcc_1_1000_n_cov.RData')
save(larsen_gcc_1_250_n_cov, file = 'larsen_gcc_1_250_n_cov.RData')
save(larsen_gcc_1_500_n_cov, file = 'larsen_gcc_1_500_n_cov.RData')
save(larsen_gcc_1_1000_n_cov, file = 'larsen_gcc_1_1000_n_cov.RData')
save(loo_gcc_1_250_n_cov, file = 'loo_gcc_1_250_n_cov.RData')
save(loo_gcc_1_500_n_cov, file = 'loo_gcc_1_500_n_cov.RData')
save(loo_gcc_1_1000_n_cov, file = 'loo_gcc_1_1000_n_cov.RData')
save(clustering_gcc_1_250_n_cov, file = 'clustering_gcc_1_250_n_cov.RData')
save(clustering_gcc_1_500_n_cov, file = 'clustering_gcc_1_500_n_cov.RData')
save(clustering_gcc_1_1000_n_cov, file = 'clustering_gcc_1_1000_n_cov.RData')
save(time_gcc_1_250_n_cov, file = 'time_gcc_1_250_n_cov.RData')
save(time_gcc_1_500_n_cov, file = 'time_gcc_1_500_n_cov.RData')
save(time_gcc_1_1000_n_cov, file = 'time_gcc_1_1000_n_cov.RData')







# KST

ari_kst_1_250_n_cov <- numeric()
jaccard_kst_1_250_n_cov <- numeric()
larsen_kst_1_250_n_cov <- numeric()
loo_kst_1_250_n_cov <- numeric()
ari_kst_1_500_n_cov <- numeric()
jaccard_kst_1_500_n_cov <- numeric()
larsen_kst_1_500_n_cov <- numeric()
loo_kst_1_500_n_cov <- numeric()
ari_kst_1_1000_n_cov <- numeric()
jaccard_kst_1_1000_n_cov <- numeric()
larsen_kst_1_1000_n_cov <- numeric()
loo_kst_1_1000_n_cov <- numeric()
clustering_kst_1_250_n_cov <- list()
clustering_kst_1_500_n_cov <- list()
clustering_kst_1_1000_n_cov <- list()
time_kst_1_250_n_cov <- numeric()
time_kst_1_500_n_cov <- numeric()
time_kst_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  
  for (i in  (1 : B)) {
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    foreach  (p = 1 : (n*K)) %:%
      foreach (j = 1 : (n*K)) %dopar% {
        a <- j_divergence(cluster[[i]][[p]], cluster[[i]][[j]])
      }
    
    for (i in 1 : (n * K)) {
      for (j in i : (n * K)) {
        
        dis_matrix[i, j] <- a[[i]][[j]]
        dis_matrix[j, i] <- dis_matrix[i, j]
        
      }
    }
    
    
    diag(dis_matrix) <- 0 # Numerical error 
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_kst_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_kst_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_kst_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_kst_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_kst_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_kst_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_kst_1_250_n_cov, file = 'ari_kst_1_250_n_cov.RData')
save(ari_kst_1_500_n_cov, file = 'ari_kst_1_500_n_cov.RData')
save(ari_kst_1_1000_n_cov, file = 'ari_kst_1_1000_n_cov.RData')
save(jaccard_kst_1_250_n_cov, file = 'jaccard_kst_1_250_n_cov.RData')
save(jaccard_kst_1_500_n_cov, file = 'jaccard_kst_1_500_n_cov.RData')
save(jaccard_kst_1_1000_n_cov, file = 'jaccard_kst_1_1000_n_cov.RData')
save(larsen_kst_1_250_n_cov, file = 'larsen_kst_1_250_n_cov.RData')
save(larsen_kst_1_500_n_cov, file = 'larsen_kst_1_500_n_cov.RData')
save(larsen_kst_1_1000_n_cov, file = 'larsen_kst_1_1000_n_cov.RData')
save(loo_kst_1_250_n_cov, file = 'loo_kst_1_250_n_cov.RData')
save(loo_kst_1_500_n_cov, file = 'loo_kst_1_500_n_cov.RData')
save(loo_kst_1_1000_n_cov, file = 'loo_kst_1_1000_n_cov.RData')
save(clustering_kst_1_250_n_cov, file = 'clustering_kst_1_250_n_cov.RData')
save(clustering_kst_1_500_n_cov, file = 'clustering_kst_1_500_n_cov.RData')
save(clustering_kst_1_1000_n_cov, file = 'clustering_kst_1_1000_n_cov.RData')
save(time_kst_1_250_n_cov, file = 'time_kst_1_250_n_cov.RData')
save(time_kst_1_500_n_cov, file = 'time_kst_1_500_n_cov.RData')
save(time_kst_1_1000_n_cov, file = 'time_kst_1_1000_n_cov.RData')









# DTW1 WITH PAM 

ari_dtw1_1_250_n_cov <- numeric()
jaccard_dtw1_1_250_n_cov <- numeric()
larsen_dtw1_1_250_n_cov <- numeric()
loo_dtw1_1_250_n_cov <- numeric()
ari_dtw1_1_500_n_cov <- numeric()
jaccard_dtw1_1_500_n_cov <- numeric()
larsen_dtw1_1_500_n_cov <- numeric()
loo_dtw1_1_500_n_cov <- numeric()
ari_dtw1_1_1000_n_cov <- numeric()
jaccard_dtw1_1_1000_n_cov <- numeric()
larsen_dtw1_1_1000_n_cov <- numeric()
loo_dtw1_1_1000_n_cov <- numeric()
clustering_dtw1_1_250_n_cov <- list()
clustering_dtw1_1_500_n_cov <- list()
clustering_dtw1_1_1000_n_cov <- list()
time_dtw1_1_250_n_cov <- numeric()
time_dtw1_1_500_n_cov <- numeric()
time_dtw1_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  
  for (i in  (1 : B)) {
    
    
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    foreach  (p = 1 : (n*K)) %:%
      foreach (j = 1 : (n*K)) %dopar% {
        a <- dtw_mts_i(cluster[[i]][[p]], cluster[[i]][[j]])
      }
    
    for (i in 1 : (n * K)) {
      for (j in i : (n * K)) {
        
        dis_matrix[i, j] <- a[[i]][[j]]
        dis_matrix[j, i] <- dis_matrix[i, j]
        
      }
    }
    
    
    diag(dis_matrix) <- 0 # Numerical error 
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_dtw1_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_dtw1_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_dtw1_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_dtw1_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_dtw1_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_dtw1_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_dtw1_1_250_n_cov, file = 'ari_dtw1_1_250_n_cov.RData')
save(ari_dtw1_1_500_n_cov, file = 'ari_dtw1_1_500_n_cov.RData')
save(ari_dtw1_1_1000_n_cov, file = 'ari_dtw1_1_1000_n_cov.RData')
save(jaccard_dtw1_1_250_n_cov, file = 'jaccard_dtw1_1_250_n_cov.RData')
save(jaccard_dtw1_1_500_n_cov, file = 'jaccard_dtw1_1_500_n_cov.RData')
save(jaccard_dtw1_1_1000_n_cov, file = 'jaccard_dtw1_1_1000_n_cov.RData')
save(larsen_dtw1_1_250_n_cov, file = 'larsen_dtw1_1_250_n_cov.RData')
save(larsen_dtw1_1_500_n_cov, file = 'larsen_dtw1_1_500_n_cov.RData')
save(larsen_dtw1_1_1000_n_cov, file = 'larsen_dtw1_1_1000_n_cov.RData')
save(loo_dtw1_1_250_n_cov, file = 'loo_dtw1_1_250_n_cov.RData')
save(loo_dtw1_1_500_n_cov, file = 'loo_dtw1_1_500_n_cov.RData')
save(loo_dtw1_1_1000_n_cov, file = 'loo_dtw1_1_1000_n_cov.RData')
save(clustering_dtw1_1_250_n_cov, file = 'clustering_dtw1_1_250_n_cov.RData')
save(clustering_dtw1_1_500_n_cov, file = 'clustering_dtw1_1_500_n_cov.RData')
save(clustering_dtw1_1_1000_n_cov, file = 'clustering_dtw1_1_1000_n_cov.RData')
save(time_dtw1_1_250_n_cov, file = 'time_dtw1_1_250_n_cov.RData')
save(time_dtw1_1_500_n_cov, file = 'time_dtw1_1_500_n_cov.RData')
save(time_dtw1_1_1000_n_cov, file = 'time_dtw1_1_1000_n_cov.RData')






# DTW2 WITH PAM 

ari_dtw2_1_250_n_cov <- numeric()
jaccard_dtw2_1_250_n_cov <- numeric()
larsen_dtw2_1_250_n_cov <- numeric()
loo_dtw2_1_250_n_cov <- numeric()
ari_dtw2_1_500_n_cov <- numeric()
jaccard_dtw2_1_500_n_cov <- numeric()
larsen_dtw2_1_500_n_cov <- numeric()
loo_dtw2_1_500_n_cov <- numeric()
ari_dtw2_1_1000_n_cov <- numeric()
jaccard_dtw2_1_1000_n_cov <- numeric()
larsen_dtw2_1_1000_n_cov <- numeric()
loo_dtw2_1_1000_n_cov <- numeric()
clustering_dtw2_1_250_n_cov <- list()
clustering_dtw2_1_500_n_cov <- list()
clustering_dtw2_1_1000_n_cov <- list()
time_dtw2_1_250_n_cov <- numeric()
time_dtw2_1_500_n_cov <- numeric()
time_dtw2_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  
  for (i in  (1 : B)) {
    
    
    
    for (s1 in 1 : (n*K)) {
      for (v1 in 1 : (n*K)) {
        p <- ncol(cluster[[1]][[1]])
        d_matrix <- matrix(0, p, p)
        foreach (s2 = 1 : p) %:%
          foreach (v2 = 1 : p) %dopar% {
            
            d_matrix[s2, v2] <- EuclideanDistance(cluster[[i]][[s1]][,s2], cluster[[i]][[v1]][,v2]) 
          }
        
        dis_matrix[s1, v1] <- dtw(d_matrix, distance.only = TRUE)$normalizedDistance
        # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'sakoechiba')$normalizedDistance
        # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'itakura')$normalizedDistance
      }
    }
    
    
    diag(dis_matrix) <- 0 # Numerical error 
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_dtw2_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_dtw2_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_dtw2_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_dtw2_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_dtw2_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_dtw2_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_dtw2_1_250_n_cov, file = 'ari_dtw2_1_250_n_cov.RData')
save(ari_dtw2_1_500_n_cov, file = 'ari_dtw2_1_500_n_cov.RData')
save(ari_dtw2_1_1000_n_cov, file = 'ari_dtw2_1_1000_n_cov.RData')
save(jaccard_dtw2_1_250_n_cov, file = 'jaccard_dtw2_1_250_n_cov.RData')
save(jaccard_dtw2_1_500_n_cov, file = 'jaccard_dtw2_1_500_n_cov.RData')
save(jaccard_dtw2_1_1000_n_cov, file = 'jaccard_dtw2_1_1000_n_cov.RData')
save(larsen_dtw2_1_250_n_cov, file = 'larsen_dtw2_1_250_n_cov.RData')
save(larsen_dtw2_1_500_n_cov, file = 'larsen_dtw2_1_500_n_cov.RData')
save(larsen_dtw2_1_1000_n_cov, file = 'larsen_dtw2_1_1000_n_cov.RData')
save(loo_dtw2_1_250_n_cov, file = 'loo_dtw2_1_250_n_cov.RData')
save(loo_dtw2_1_500_n_cov, file = 'loo_dtw2_1_500_n_cov.RData')
save(loo_dtw2_1_1000_n_cov, file = 'loo_dtw2_1_1000_n_cov.RData')
save(clustering_dtw2_1_250_n_cov, file = 'clustering_dtw2_1_250_n_cov.RData')
save(clustering_dtw2_1_500_n_cov, file = 'clustering_dtw2_1_500_n_cov.RData')
save(clustering_dtw2_1_1000_n_cov, file = 'clustering_dtw2_1_1000_n_cov.RData')
save(time_dtw2_1_250_n_cov, file = 'time_dtw2_1_250_n_cov.RData')
save(time_dtw2_1_500_n_cov, file = 'time_dtw2_1_500_n_cov.RData')
save(time_dtw2_1_1000_n_cov, file = 'time_dtw2_1_1000_n_cov.RData')







# PCA 

ari_pca_1_250_n_cov <- numeric()
jaccard_pca_1_250_n_cov <- numeric()
larsen_pca_1_250_n_cov <- numeric()
loo_pca_1_250_n_cov <- numeric()
ari_pca_1_500_n_cov <- numeric()
jaccard_pca_1_500_n_cov <- numeric()
larsen_pca_1_500_n_cov <- numeric()
loo_pca_1_500_n_cov <- numeric()
ari_pca_1_1000_n_cov <- numeric()
jaccard_pca_1_1000_n_cov <- numeric()
larsen_pca_1_1000_n_cov <- numeric()
loo_pca_1_1000_n_cov <- numeric()
clustering_pca_1_250_n_cov <- list()
clustering_pca_1_500_n_cov <- list()
clustering_pca_1_1000_n_cov <- list()
time_pca_1_250_n_cov <- numeric()
time_pca_1_500_n_cov <- numeric()
time_pca_1_1000_n_cov <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  
  foreach (i =  (1 : B)) %dopar% {
    
    
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    clustercov <- list()
    for (m in 1 : length(cluster[[i]])) {
      clustercov[[m]] <- cov(cluster[[i]][[m]])
    }
    
    dis_matrix <- 1 - PCAsimilarity(clustercov)
    
    for (p in 1 : nrow(dis_matrix)) {
      for (v in p : ncol(dis_matrix)) {
        dis_matrix[p, v] <- dis_matrix[v, p]
      }
    }
    
    clustering <- pam(dis_matrix, K)$cluster
    
    if (l == 250) {
      
      ari_pca_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_1_250_n_cov[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_pca_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_pca_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_pca_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_pca_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_pca_1_1000_n_cov <- start_time - end_time
    
  }
  
}





save(ari_pca_1_250_n_cov, file = 'ari_pca_1_250_n_cov.RData')
save(ari_pca_1_500_n_cov, file = 'ari_pca_1_500_n_cov.RData')
save(ari_pca_1_1000_n_cov, file = 'ari_pca_1_1000_n_cov.RData')
save(jaccard_pca_1_250_n_cov, file = 'jaccard_pca_1_250_n_cov.RData')
save(jaccard_pca_1_500_n_cov, file = 'jaccard_pca_1_500_n_cov.RData')
save(jaccard_pca_1_1000_n_cov, file = 'jaccard_pca_1_1000_n_cov.RData')
save(larsen_pca_1_250_n_cov, file = 'larsen_pca_1_250_n_cov.RData')
save(larsen_pca_1_500_n_cov, file = 'larsen_pca_1_500_n_cov.RData')
save(larsen_pca_1_1000_n_cov, file = 'larsen_pca_1_1000_n_cov.RData')
save(loo_pca_1_250_n_cov, file = 'loo_pca_1_250_n_cov.RData')
save(loo_pca_1_500_n_cov, file = 'loo_pca_1_500_n_cov.RData')
save(loo_pca_1_1000_n_cov, file = 'loo_pca_1_1000_n_cov.RData')
save(clustering_pca_1_250_n_cov, file = 'clustering_pca_1_250_n_cov.RData')
save(clustering_pca_1_500_n_cov, file = 'clustering_pca_1_500_n_cov.RData')
save(clustering_pca_1_1000_n_cov, file = 'clustering_pca_1_1000_n_cov.RData')
save(time_pca_1_250_n_cov, file = 'time_pca_1_250_n_cov.RData')
save(time_pca_1_500_n_cov, file = 'time_pca_1_500_n_cov.RData')
save(time_pca_1_1000_n_cov, file = 'time_pca_1_1000_n_cov.RData')










# Maharaj

ari_mah_1_1000_n_cov <- numeric()
jaccard_mah_1_1000_n_cov <- numeric()
larsen_mah_1_1000_n_cov <- numeric()
loo_mah_1_1000_n_cov <- numeric()
clustering_mah_1_1000_n_cov <- list()
time_mah_1_1000_n_cov <- numeric()
ari_mah_1_500_n_cov <- numeric()
jaccard_mah_1_500_n_cov <- numeric()
larsen_mah_1_500_n_cov <- numeric()
loo_mah_1_500_n_cov <- numeric()
clustering_mah_1_500_n_cov <- list()
time_mah_1_500_n_cov <- numeric()
ari_mah_1_250_n_cov <- numeric()
jaccard_mah_1_250_n_cov <- numeric()
larsen_mah_1_250_n_cov <- numeric()
loo_mah_1_250_n_cov <- numeric()
clustering_mah_1_250_n_cov <- list()
time_mah_1_250_n_cov <- numeric()





for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_1_', l, '_n_cov.RData'))
  clusterExport(c1, c('cluster', 'maharaj'))
  features <- parLapply(c1, cluster, maharaj)
  
  
  for (i in (1:length(features))) {
    
    dis_matrix <- proxy::dist(features[[i]], EuclideanDistance)  
    clustering <- pam(dis_matrix, K)$cluster
  }
  
  if (l == 250) {
    
    ari_mah_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_1_250_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_1_250_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_1_250_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_1_250_n_cov[[i]] <- pam(dis_matrix, K)$cluster
    
  }
  
  if (l == 500) {
    
    ari_mah_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_1_500_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_1_500_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_1_500_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_1_500_n_cov[[i]] <- pam(dis_matrix, K)$cluster
    
  }
  
  if (l == 1000) {
    
    ari_mah_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_1_1000_n_cov[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_1_1000_n_cov[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_1_1000_n_cov[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_1_1000_n_cov[[i]] <- pam(dis_matrix, K)$cluster
    
  }
  
  
  
  
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_mah_1_250_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_mah_1_500_n_cov <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_mah_1_1000_n_cov <- start_time - end_time
    
  }
  
  
}

save(ari_mah_1_1000_n_cov, file = 'ari_mah_1_1000_n_cov.RData')
save(jaccard_mah_1_1000_n_cov, file = 'jaccard_mah_1_1000_n_cov.RData')
save(larsen_mah_1_1000_n_cov, file = 'larsen_mah_1_1000_n_cov.RData')
save(loo_mah_1_1000_n_cov, file = 'loo_mah_1_1000_n_cov.RData')
save(clustering_mah_1_1000_n_cov, file = 'clustering_mah_1_1000_n_cov.RData')
save(time_mah_1_1000_n_cov, file = 'time_mah_1_1000_n_cov.RData')
save(ari_mah_1_500_n_cov, file = 'ari_mah_1_500_n_cov.RData')
save(jaccard_mah_1_500_n_cov, file = 'jaccard_mah_1_500_n_cov.RData')
save(larsen_mah_1_500_n_cov, file = 'larsen_mah_1_500_n_cov.RData')
save(loo_mah_1_500_n_cov, file = 'loo_mah_1_500_n_cov.RData')
save(clustering_mah_1_500_n_cov, file = 'clustering_mah_1_500_n_cov.RData')
save(time_mah_1_500_n_cov, file = 'time_mah_1_500_n_cov.RData')
save(ari_mah_1_250_n_cov, file = 'ari_mah_1_250_n_cov.RData')
save(jaccard_mah_1_250_n_cov, file = 'jaccard_mah_1_250_n_cov.RData')
save(larsen_mah_1_250_n_cov, file = 'larsen_mah_1_250_n_cov.RData')
save(loo_mah_1_250_n_cov, file = 'loo_mah_1_250_n_cov.RData')
save(clustering_mah_1_250_n_cov, file = 'clustering_mah_1_250_n_cov.RData')
save(time_mah_1_250_n_cov, file = 'time_mah_1_250_n_cov.RData')
