
K <- 3
n <- 5
B <- 100
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n))
cluster <- list()

set.seed(1234)
for (l in c(250, 500, 1000)){
  
  add <- 150
  l <- l + add
  
  for (k in 1 : B){
    
    # Cluster 1
    
    cluster1 <- list()
    
    for (j in 1 : n){
      
      Xt <- list()
      Xt[[1]] <- c(0, 0)
      
      for (i in 2 : (l + 1)) {
        
        theta <- 0*diag(2)
        u12 <- runif(1)
        u21 <- runif(1)
        theta[1, 2] <- -0.5*(u12 - 0.5)
        theta[2, 1] <- -0.5*(u21 - 0.5)
        theta0 <- numeric()
        theta0[1] <- qnorm(u12, 0)
        theta0[2] <- qnorm(u21, 0)
        
        Xt[[i]] <- theta %*% Xt[[i - 1]] + theta0
        
      }
      
      cluster1[[j]] <- listTomatrix(Xt)
      cluster1[[j]] <- cluster1[[j]][((add+1):l),]
      
    }
    
    
    
    # Cluster 2 
    
    cluster2 <- list()
    
    for (j in 1 : n){
      
      Xt <- list()
      Xt[[1]] <- c(0, 0)
      
      for (i in 2 : (l + 1)) {
        
        theta <- 0*diag(2)
        u12 <- runif(1)
        u21 <- runif(1)
        theta[1, 2] <- 0.5*(u12 - 0.5)
        theta[2, 1] <- 0.5*(u21 - 0.5)
        
        theta0 <- numeric()
        theta0[1] <- qnorm(u12, 0)
        theta0[2] <- qnorm(u21, 0)
        
        Xt[[i]] <- theta %*% Xt[[i - 1]] + theta0
        
      }
      
      cluster2[[j]] <- listTomatrix(Xt)
      cluster2[[j]] <- cluster2[[j]][((add+1):l),]
      
    }
    
    
    
    # Cluster 3
    
    
    cluster3 <- list()
    
    for (j in 1 : n){
      
      Xt <- list()
      Xt[[1]] <- c(0, 0)
      
      for (i in 2 : (l + 1)) {
        
        theta <- 0*diag(2)
        u12 <- runif(1)
        u21 <- runif(1)
        theta[1, 2] <- 1.5*(u12 - 0.5)
        theta[2, 1] <- 1.5*(u21 - 0.5)
        
        theta0 <- numeric()
        theta0[1] <- qnorm(u12, 0)
        theta0[2] <- qnorm(u21, 0)
        
        Xt[[i]] <- theta %*% Xt[[i - 1]] + theta0
        
      }
      
      cluster3[[j]] <- listTomatrix(Xt)
      cluster3[[j]] <- cluster3[[j]][((add+1):l),]
      
    
      
      
    }
    
    cluster[[k]] <- c(cluster1, cluster2, cluster3)
    
  }
  
  l <- l - add
  
  # Saving the simulations
  
  
  if (l == 250){
    
    cluster_3_250_n <- cluster 
    save(cluster_3_250_n, file = 'cluster_3_250_n.RData')
    
  }
  
  if (l == 500){
    
    cluster_3_500_n<- cluster 
    save(cluster_3_500_n, file = 'cluster_3_500_n.RData')
    
  }
  
  if (l == 1000){
    
    cluster_3_1000_n <- cluster 
    save(cluster_3_1000_n, file = 'cluster_3_1000_n.RData')
    
  }
  
}


# Loading the data 

n <- 5  # Number of series per cluster
K <- 3 # Number of clusters
B <- 100 # Number of Monte Carlos trials
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n)) # Ground Truth 

c1 <- makeCluster(7) # Making a cluster object
registerDoParallel(cores = 7)
load('cluster_3_250_n.RData')
load('cluster_3_500_n.RData')
load('cluster_3_1000_n.RData')

# QC2

ari_qc2_3_250_n <- numeric()
jaccard_qc2_3_250_n <- numeric()
larsen_qc2_3_250_n <- numeric()
loo_qc2_3_250_n <- numeric()
ari_qc2_3_500_n <- numeric()
jaccard_qc2_3_500_n <- numeric()
larsen_qc2_3_500_n <- numeric()
loo_qc2_3_500_n <- numeric()
ari_qc2_3_1000_n <- numeric()
jaccard_qc2_3_1000_n <- numeric()
larsen_qc2_3_1000_n <- numeric()
loo_qc2_3_1000_n <- numeric()
clustering_qc2_3_250_n <- list()
clustering_qc2_3_500_n <- list()
clustering_qc2_3_1000_n <- list()
time_qc2_3_250_n <- numeric()
time_qc2_3_500_n <- numeric()
time_qc2_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  clusterExport(c1, c('cluster', 'quantile_quantities_re_im'))
  
  for (i in  (1 : B)) {
    
    
    coherence2 <- parLapply(c1, cluster[[i]], quantile_quantities_re_im)
    dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
    clustering <- pam(dis_matrix, K)$cluster
    external_validation(ground_truth, clustering)
    
    if (l == 250) {
      
      ari_qc2_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_qc2_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_qc2_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_qc2_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_qc2_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_qc2_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_qc2_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_qc2_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_qc2_3_500_2 <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_qc2_3_1000_n <- start_time - end_time
    
  }
  
  
  
}


end_time <- Sys.time()
time_qc2_3_250_n <- start_time - end_time
time_qc2_3_500_n <- start_time - end_time
time_qc2_3_1000_n <- start_time - end_time
save(ari_qc2_3_250_n, file = 'ari_qc2_3_250_n.RData')
save(ari_qc2_3_500_n, file = 'ari_qc2_3_500_n.RData')
save(ari_qc2_3_1000_n, file = 'ari_qc2_3_1000_n.RData')
save(jaccard_qc2_3_250_n, file = 'jaccard_qc2_3_250_n.RData')
save(jaccard_qc2_3_500_n, file = 'jaccard_qc2_3_500_n.RData')
save(jaccard_qc2_3_1000_n, file = 'jaccard_qc2_3_1000_n.RData')
save(larsen_qc2_3_250_n, file = 'larsen_qc2_3_250_n.RData')
save(larsen_qc2_3_500_n, file = 'larsen_qc2_3_500_n.RData')
save(larsen_qc2_3_1000_n, file = 'larsen_qc2_3_1000_n.RData')
save(loo_qc2_3_250_n, file = 'loo_qc2_3_250_n.RData')
save(loo_qc2_3_500_n, file = 'loo_qc2_3_500_n.RData')
save(loo_qc2_3_1000_n, file = 'loo_qc2_3_1000_n.RData')
save(clustering_qc2_3_250_n, file = 'clustering_qc2_3_250_n.RData')
save(clustering_qc2_3_500_n, file = 'clustering_qc2_3_500_n.RData')
save(clustering_qc2_3_1000_n, file = 'clustering_qc2_3_1000_n.RData')
save(time_qc2_3_250_n, file = 'time_qc2_3_250_n.RData')
save(time_qc2_3_500_n, file = 'time_qc2_3_500_n.RData')
save(time_qc2_3_1000_n, file = 'time_qc2_3_1000_n.RData')











# DM

ari_w_3_250_n <- numeric()
jaccard_w_3_250_n <- numeric()
larsen_w_3_250_n <- numeric()
loo_w_3_250_n <- numeric()
ari_w_3_500_n <- numeric()
jaccard_w_3_500_n <- numeric()
larsen_w_3_500_n <- numeric()
loo_w_3_500_n <- numeric()
ari_w_3_1000_n <- numeric()
jaccard_w_3_1000_n <- numeric()
larsen_w_3_1000_n <- numeric()
loo_w_3_1000_n <- numeric()
clustering_w_3_250_n <- list()
clustering_w_3_500_n <- list()
clustering_w_3_1000_n <- list()
time_w_3_250_n <- numeric()
time_w_3_500_n <- numeric()
time_w_3_1000_n <- numeric()



wf <- 'd4'

for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
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
    external_validation(ground_truth, clustering)
    
    
    if (l == 250) {
      
      ari_w_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_w_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_w_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_w_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_w_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_w_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_w_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  
  
  if (l == 250) {
    
    time_w_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_w_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_w_3_1000_n <- start_time - end_time
    
  }
  
  
}


save(ari_w_3_250_n, file = 'ari_w_3_250_n.RData')
save(ari_w_3_500_n, file = 'ari_w_3_500_n.RData')
save(ari_w_3_1000_n, file = 'ari_w_3_1000_n.RData')
save(jaccard_w_3_250_n, file = 'jaccard_w_3_250_n.RData')
save(jaccard_w_3_500_n, file = 'jaccard_w_3_500_n.RData')
save(jaccard_w_3_1000_n, file = 'jaccard_w_3_1000_n.RData')
save(larsen_w_3_250_n, file = 'larsen_w_3_250_n.RData')
save(larsen_w_3_500_n, file = 'larsen_w_3_500_n.RData')
save(larsen_w_3_1000_n, file = 'larsen_w_3_1000_n.RData')
save(loo_w_3_250_n, file = 'loo_w_3_250_n.RData')
save(loo_w_3_500_n, file = 'loo_w_3_500_n.RData')
save(loo_w_3_1000_n, file = 'loo_w_3_1000_n.RData')
save(clustering_w_3_250_n, file = 'clustering_w_3_250_n.RData')
save(clustering_w_3_500_n, file = 'clustering_w_3_500_n.RData')
save(clustering_w_3_1000_n, file = 'clustering_w_3_1000_n.RData')
save(time_w_3_250_n, file = 'time_w_3_250_n.RData')
save(time_w_3_500_n, file = 'time_w_3_500_n.RData')
save(time_w_3_1000_n, file = 'time_w_3_1000_n.RData')












# AP

ari_gcc_3_250_n <- numeric()
jaccard_gcc_3_250_n <- numeric()
larsen_gcc_3_250_n <- numeric()
loo_gcc_3_250_n <- numeric()
ari_gcc_3_500_n <- numeric()
jaccard_gcc_3_500_n <- numeric()
larsen_gcc_3_500_n <- numeric()
loo_gcc_3_500_n <- numeric()
ari_gcc_3_1000_n <- numeric()
jaccard_gcc_3_1000_n <- numeric()
larsen_gcc_3_1000_n <- numeric()
loo_gcc_3_1000_n <- numeric()
clustering_gcc_3_250_n <- list()
clustering_gcc_3_500_n <- list()
clustering_gcc_3_1000_n <- list()
time_gcc_3_250_n <- numeric()
time_gcc_3_500_n <- numeric()
time_gcc_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  clusterExport(c1, c('cluster', 'gcc_features_mts'))
  
  for (i in  (1 : B)) {
    
    features <- parLapply(c1, cluster[[i]], gcc_features_mts)
    dis_matrix <- proxy::dist(features, EuclideanDistance)
    clustering <- pam(dis_matrix, K)$cluster
    external_validation(ground_truth, clustering)
    
    
    if (l == 250) {
      
      ari_gcc_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_gcc_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_gcc_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_gcc_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_gcc_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_gcc_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_gcc_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_gcc_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_gcc_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_gcc_3_1000_n <- start_time - end_time
    
  }
  
}





save(ari_gcc_3_250_n, file = 'ari_gcc_3_250_n.RData')
save(ari_gcc_3_500_n, file = 'ari_gcc_3_500_n.RData')
save(ari_gcc_3_1000_n, file = 'ari_gcc_3_1000_n.RData')
save(jaccard_gcc_3_250_n, file = 'jaccard_gcc_3_250_n.RData')
save(jaccard_gcc_3_500_n, file = 'jaccard_gcc_3_500_n.RData')
save(jaccard_gcc_3_1000_n, file = 'jaccard_gcc_3_1000_n.RData')
save(larsen_gcc_3_250_n, file = 'larsen_gcc_3_250_n.RData')
save(larsen_gcc_3_500_n, file = 'larsen_gcc_3_500_n.RData')
save(larsen_gcc_3_1000_n, file = 'larsen_gcc_3_1000_n.RData')
save(loo_gcc_3_250_n, file = 'loo_gcc_3_250_n.RData')
save(loo_gcc_3_500_n, file = 'loo_gcc_3_500_n.RData')
save(loo_gcc_3_1000_n, file = 'loo_gcc_3_1000_n.RData')
save(clustering_gcc_3_250_n, file = 'clustering_gcc_3_250_n.RData')
save(clustering_gcc_3_500_n, file = 'clustering_gcc_3_500_n.RData')
save(clustering_gcc_3_1000_n, file = 'clustering_gcc_3_1000_n.RData')
save(time_gcc_3_250_n, file = 'time_gcc_3_250_n.RData')
save(time_gcc_3_500_n, file = 'time_gcc_3_500_n.RData')
save(time_gcc_3_1000_n, file = 'time_gcc_3_1000_n.RData')












# KST

ari_kst_3_250_n <- numeric()
jaccard_kst_3_250_n <- numeric()
larsen_kst_3_250_n <- numeric()
loo_kst_3_250_n <- numeric()
ari_kst_3_500_n <- numeric()
jaccard_kst_3_500_n <- numeric()
larsen_kst_3_500_n <- numeric()
loo_kst_3_500_n <- numeric()
ari_kst_3_1000_n <- numeric()
jaccard_kst_3_1000_n <- numeric()
larsen_kst_3_1000_n <- numeric()
loo_kst_3_1000_n <- numeric()
clustering_kst_3_250_n <- list()
clustering_kst_3_500_n <- list()
clustering_kst_3_1000_n <- list()
time_kst_3_250_n <- numeric()
time_kst_3_500_n <- numeric()
time_kst_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  
  for (i in  (1 : B)) {
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    a <- foreach  (p = 1 : (n*K)) %:%
      foreach (j = 1 : (n*K)) %dopar% {
        j_divergence(cluster[[i]][[p]], cluster[[i]][[j]])
      }
    
    for (v2 in (1 : (n * K))) {
      for (v1 in ((v2) : (n * K))){
        
        dis_matrix[v1, v2] <- a[[v1]][[v2]]
        
      }
      
    }
    
    
    diag(dis_matrix) <- 0 # Numerical error
    dis_matrix <- as.dist(dis_matrix)
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_kst_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_kst_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_kst_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_kst_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_kst_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_kst_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_kst_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_kst_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_kst_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_kst_3_1000_n <- start_time - end_time
    
  }
  
}





save(ari_kst_3_250_n, file = 'ari_kst_3_250_n.RData')
save(ari_kst_3_500_n, file = 'ari_kst_3_500_n.RData')
save(ari_kst_3_1000_n, file = 'ari_kst_3_1000_n.RData')
save(jaccard_kst_3_250_n, file = 'jaccard_kst_3_250_n.RData')
save(jaccard_kst_3_500_n, file = 'jaccard_kst_3_500_n.RData')
save(jaccard_kst_3_1000_n, file = 'jaccard_kst_3_1000_n.RData')
save(larsen_kst_3_250_n, file = 'larsen_kst_3_250_n.RData')
save(larsen_kst_3_500_n, file = 'larsen_kst_3_500_n.RData')
save(larsen_kst_3_1000_n, file = 'larsen_kst_3_1000_n.RData')
save(loo_kst_3_250_n, file = 'loo_kst_3_250_n.RData')
save(loo_kst_3_500_n, file = 'loo_kst_3_500_n.RData')
save(loo_kst_3_1000_n, file = 'loo_kst_3_1000_n.RData')
save(clustering_kst_3_250_n, file = 'clustering_kst_3_250_n.RData')
save(clustering_kst_3_500_n, file = 'clustering_kst_3_500_n.RData')
save(clustering_kst_3_1000_n, file = 'clustering_kst_3_1000_n.RData')
save(time_kst_3_250_n, file = 'time_kst_3_250_n.RData')
save(time_kst_3_500_n, file = 'time_kst_3_500_n.RData')
save(time_kst_3_1000_n, file = 'time_kst_3_1000_n.RData')











# DTW1 WITH PAM 

ari_dtw1_3_250_n <- numeric()
jaccard_dtw1_3_250_n <- numeric()
larsen_dtw1_3_250_n <- numeric()
loo_dtw1_3_250_n <- numeric()
ari_dtw1_3_500_n <- numeric()
jaccard_dtw1_3_500_n <- numeric()
larsen_dtw1_3_500_n <- numeric()
loo_dtw1_3_500_n <- numeric()
ari_dtw1_3_1000_n <- numeric()
jaccard_dtw1_3_1000_n <- numeric()
larsen_dtw1_3_1000_n <- numeric()
loo_dtw1_3_1000_n <- numeric()
clustering_dtw1_3_250_n <- list()
clustering_dtw1_3_500_n <- list()
clustering_dtw1_3_1000_n <- list()
time_dtw1_3_250_n <- numeric()
time_dtw1_3_500_n <- numeric()
time_dtw1_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  
  for (i in  (1 : B)) {
    
    
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    a <- foreach  (p = 1 : (n*K)) %:%
      foreach (j = 1 : (n*K)) %dopar% {
        dtw_mts_i(cluster[[i]][[p]], cluster[[i]][[j]])
      }
    
    for (v2 in (1 : (n * K))) {
      for (v1 in ((v2) : (n * K))){
        
        dis_matrix[v1, v2] <- a[[v1]][[v2]]
        
      }
      
    }
    
    diag(dis_matrix) <- 0 # Numerical error 
    dis_matrix <- as.dist(dis_matrix)
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_dtw1_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_dtw1_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_dtw1_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw1_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw1_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw1_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw1_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_dtw1_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_dtw1_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_dtw1_3_1000_n <- start_time - end_time
    
  }
  
}





save(ari_dtw1_3_250_n, file = 'ari_dtw1_3_250_n.RData')
save(ari_dtw1_3_500_n, file = 'ari_dtw1_3_500_n.RData')
save(ari_dtw1_3_1000_n, file = 'ari_dtw1_3_1000_n.RData')
save(jaccard_dtw1_3_250_n, file = 'jaccard_dtw1_3_250_n.RData')
save(jaccard_dtw1_3_500_n, file = 'jaccard_dtw1_3_500_n.RData')
save(jaccard_dtw1_3_1000_n, file = 'jaccard_dtw1_3_1000_n.RData')
save(larsen_dtw1_3_250_n, file = 'larsen_dtw1_3_250_n.RData')
save(larsen_dtw1_3_500_n, file = 'larsen_dtw1_3_500_n.RData')
save(larsen_dtw1_3_1000_n, file = 'larsen_dtw1_3_1000_n.RData')
save(loo_dtw1_3_250_n, file = 'loo_dtw1_3_250_n.RData')
save(loo_dtw1_3_500_n, file = 'loo_dtw1_3_500_n.RData')
save(loo_dtw1_3_1000_n, file = 'loo_dtw1_3_1000_n.RData')
save(clustering_dtw1_3_250_n, file = 'clustering_dtw1_3_250_n.RData')
save(clustering_dtw1_3_500_n, file = 'clustering_dtw1_3_500_n.RData')
save(clustering_dtw1_3_1000_n, file = 'clustering_dtw1_3_1000_n.RData')
save(time_dtw1_3_250_n, file = 'time_dtw1_3_250_n.RData')
save(time_dtw1_3_500_n, file = 'time_dtw1_3_500_n.RData')
save(time_dtw1_3_1000_n, file = 'time_dtw1_3_1000_n.RData')










# DTW2 WITH PAM 


ari_dtw2_3_250_n <- numeric()
jaccard_dtw2_3_250_n <- numeric()
larsen_dtw2_3_250_n <- numeric()
loo_dtw2_3_250_n <- numeric()
ari_dtw2_3_500_n <- numeric()
jaccard_dtw2_3_500_n <- numeric()
larsen_dtw2_3_500_n <- numeric()
loo_dtw2_3_500_n <- numeric()
ari_dtw2_3_1000_n <- numeric()
jaccard_dtw2_3_1000_n <- numeric()
larsen_dtw2_3_1000_n <- numeric()
loo_dtw2_3_1000_n <- numeric()
clustering_dtw2_3_250_n <- list()
clustering_dtw2_3_500_n <- list()
clustering_dtw2_3_1000_n <- list()
time_dtw2_3_250_n <- numeric()
time_dtw2_3_500_n <- numeric()
time_dtw2_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  
  for (i in  (1 : B)) {
    
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    a <- foreach  (s1 = 1 : (n*K)) %:% 
      foreach (s2 = 1 : (n*K)) %dopar% {
        
        
        proxy::dist(t(cluster[[i]][[s1]]), t(cluster[[i]][[s2]])) 
        # dis_matrix[s1, v1] <- dtw(d_matrix, distance.only = TRUE)$normalizedDistance
        # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'sakoechiba')$normalizedDistance
        # dis_matrix[i, j] <- dtw(d_matrix, window.type = 'itakura')$normalizedDistance
      }
    
    
    for (v2 in (1 : (n * K))) {
      for (v1 in ((v2) : (n * K))){
        
        dis_matrix[v1, v2] <- dtw(a[[v1]][[v2]], distance.only = T)$normalizedDistance
        
      }
      
    }
    
    
    
    
    
    diag(dis_matrix) <- 0 # Numerical error 
    dis_matrix <- as.dist(dis_matrix)
    clustering <- pam(dis_matrix, K)$cluster
    
    
    if (l == 250) {
      
      ari_dtw2_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_dtw2_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_dtw2_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_dtw2_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_dtw2_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_dtw2_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_dtw2_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_dtw2_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_dtw2_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_dtw2_3_1000_n <- start_time - end_time
    
  }
  
}





save(ari_dtw2_3_250_n, file = 'ari_dtw2_3_250_n.RData')
save(ari_dtw2_3_500_n, file = 'ari_dtw2_3_500_n.RData')
save(ari_dtw2_3_1000_n, file = 'ari_dtw2_3_1000_n.RData')
save(jaccard_dtw2_3_250_n, file = 'jaccard_dtw2_3_250_n.RData')
save(jaccard_dtw2_3_500_n, file = 'jaccard_dtw2_3_500_n.RData')
save(jaccard_dtw2_3_1000_n, file = 'jaccard_dtw2_3_1000_n.RData')
save(larsen_dtw2_3_250_n, file = 'larsen_dtw2_3_250_n.RData')
save(larsen_dtw2_3_500_n, file = 'larsen_dtw2_3_500_n.RData')
save(larsen_dtw2_3_1000_n, file = 'larsen_dtw2_3_1000_n.RData')
save(loo_dtw2_3_250_n, file = 'loo_dtw2_3_250_n.RData')
save(loo_dtw2_3_500_n, file = 'loo_dtw2_3_500_n.RData')
save(loo_dtw2_3_1000_n, file = 'loo_dtw2_3_1000_n.RData')
save(clustering_dtw2_3_250_n, file = 'clustering_dtw2_3_250_n.RData')
save(clustering_dtw2_3_500_n, file = 'clustering_dtw2_3_500_n.RData')
save(clustering_dtw2_3_1000_n, file = 'clustering_dtw2_3_1000_n.RData')
save(time_dtw2_3_250_n, file = 'time_dtw2_3_250_n.RData')
save(time_dtw2_3_500_n, file = 'time_dtw2_3_500_n.RData')
save(time_dtw2_3_1000_n, file = 'time_dtw2_3_1000_n.RData')










# PCA 

ari_pca_3_250_n <- numeric()
jaccard_pca_3_250_n <- numeric()
larsen_pca_3_250_n <- numeric()
loo_pca_3_250_n <- numeric()
ari_pca_3_500_n <- numeric()
jaccard_pca_3_500_n <- numeric()
larsen_pca_3_500_n <- numeric()
loo_pca_3_500_n <- numeric()
ari_pca_3_1000_n <- numeric()
jaccard_pca_3_1000_n <- numeric()
larsen_pca_3_1000_n <- numeric()
loo_pca_3_1000_n <- numeric()
clustering_pca_3_250_n <- list()
clustering_pca_3_500_n <- list()
clustering_pca_3_1000_n <- list()
time_pca_3_250_n <- numeric()
time_pca_3_500_n <- numeric()
time_pca_3_1000_n <- numeric()


for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  
  for (i in  (1 : B)) {
    
    
    
    dis_matrix <- matrix(0, (n*K), (n*K))
    clustercov <- foreach (m = 1 : length(cluster[[i]])) %dopar% {
       cov(cluster[[i]][[m]])
    }
    
    dis_matrix <- 1 - PCAsimilarity(clustercov)
    dis_matrix[col(dis_matrix) >= row(dis_matrix)] <- 0
    dis_matrix <- as.dist(dis_matrix)
    clustering <- pam(dis_matrix, K)$cluster
    
    if (l == 250) {
      
      ari_pca_3_250_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_3_250_n[[i]] <- clustering
      
    }
    
    
    if (l == 500) {
      
      ari_pca_3_500_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
    
    if (l == 1000) {
      
      ari_pca_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
      jaccard_pca_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
      larsen_pca_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
      loo_pca_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
      clustering_pca_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
      
    }
    
  }
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_pca_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_pca_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_pca_3_1000_n <- start_time - end_time
    
  }
  
}





save(ari_pca_3_250_n, file = 'ari_pca_3_250_n.RData')
save(ari_pca_3_500_n, file = 'ari_pca_3_500_n.RData')
save(ari_pca_3_1000_n, file = 'ari_pca_3_1000_n.RData')
save(jaccard_pca_3_250_n, file = 'jaccard_pca_3_250_n.RData')
save(jaccard_pca_3_500_n, file = 'jaccard_pca_3_500_n.RData')
save(jaccard_pca_3_1000_n, file = 'jaccard_pca_3_1000_n.RData')
save(larsen_pca_3_250_n, file = 'larsen_pca_3_250_n.RData')
save(larsen_pca_3_500_n, file = 'larsen_pca_3_500_n.RData')
save(larsen_pca_3_1000_n, file = 'larsen_pca_3_1000_n.RData')
save(loo_pca_3_250_n, file = 'loo_pca_3_250_n.RData')
save(loo_pca_3_500_n, file = 'loo_pca_3_500_n.RData')
save(loo_pca_3_1000_n, file = 'loo_pca_3_1000_n.RData')
save(clustering_pca_3_250_n, file = 'clustering_pca_3_250_n.RData')
save(clustering_pca_3_500_n, file = 'clustering_pca_3_500_n.RData')
save(clustering_pca_3_1000_n, file = 'clustering_pca_3_1000_n.RData')
save(time_pca_3_250_n, file = 'time_pca_3_250_n.RData')
save(time_pca_3_500_n, file = 'time_pca_3_500_n.RData')
save(time_pca_3_1000_n, file = 'time_pca_3_1000_n.RData')










# Maharaj

ari_mah_3_1000_n <- numeric()
jaccard_mah_3_1000_n <- numeric()
larsen_mah_3_1000_n <- numeric()
loo_mah_3_1000_n <- numeric()
clustering_mah_3_1000_n <- list()
time_mah_3_1000_n <- numeric()
ari_mah_3_500_n <- numeric()
jaccard_mah_3_500_n <- numeric()
larsen_mah_3_500_n <- numeric()
loo_mah_3_500_n <- numeric()
clustering_mah_3_500_n <- list()
time_mah_3_500_n <- numeric()
ari_mah_3_250_n <- numeric()
jaccard_mah_3_250_n <- numeric()
larsen_mah_3_250_n <- numeric()
loo_mah_3_250_n <- numeric()
clustering_mah_3_250_n <- list()
time_mah_3_250_n <- numeric()





for (l in c(250, 500, 1000)){
  
  start_time <- Sys.time()
  cluster <- loadRData(paste0('cluster_3_', l, '_n.RData'))
  clusterExport(c1, c('cluster', 'maharaj'))
  features <- parLapply(c1, cluster, maharaj)
  
  
  for (i in (1:length(features))) {
    
    dis_matrix <- proxy::dist(features[[i]], EuclideanDistance)  
    clustering <- pam(dis_matrix, K)$cluster
  
  
  if (l == 250) {
    
    ari_mah_3_250_n[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_3_250_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_3_250_n[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_3_250_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_3_250_n[[i]] <- pam(dis_matrix, K)$cluster
    
  }
  
  if (l == 500) {
    
    ari_mah_3_500_n[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_3_500_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_3_500_n[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_3_500_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_3_500_n[[i]] <- pam(dis_matrix, K)$cluster
    
  }
  
  if (l == 1000) {
    
    ari_mah_3_1000_n[[i]] <- external_validation(ground_truth, clustering)
    jaccard_mah_3_1000_n[[i]] <- external_validation(ground_truth, clustering, 'jaccard_index')
    larsen_mah_3_1000_n[[i]] <- cluster.evaluation(ground_truth, clustering)
    loo_mah_3_1000_n[[i]] <- loo1nn.cv(dis_matrix, ground_truth)
    clustering_mah_3_1000_n[[i]] <- pam(dis_matrix, K)$cluster
    
  }
    
  }
  
  
  
  
  
  end_time <- Sys.time()
  
  if (l == 250) {
    
    time_mah_3_250_n <- start_time - end_time
    
  }
  
  
  if (l == 500) {
    
    time_mah_3_500_n <- start_time - end_time
    
  }
  
  
  if (l == 1000) {
    
    time_mah_3_1000_n <- start_time - end_time
    
  }
  
  
}

save(ari_mah_3_1000_n, file = 'ari_mah_3_1000_n.RData')
save(jaccard_mah_3_1000_n, file = 'jaccard_mah_3_1000_n.RData')
save(larsen_mah_3_1000_n, file = 'larsen_mah_3_1000_n.RData')
save(loo_mah_3_1000_n, file = 'loo_mah_3_1000_n.RData')
save(clustering_mah_3_1000_n, file = 'clustering_mah_3_1000_n.RData')
save(time_mah_3_1000_n, file = 'time_mah_3_1000_n.RData')
save(ari_mah_3_500_n, file = 'ari_mah_3_500_n.RData')
save(jaccard_mah_3_500_n, file = 'jaccard_mah_3_500_n.RData')
save(larsen_mah_3_500_n, file = 'larsen_mah_3_500_n.RData')
save(loo_mah_3_500_n, file = 'loo_mah_3_500_n.RData')
save(clustering_mah_3_500_n, file = 'clustering_mah_3_500_n.RData')
save(time_mah_3_500_n, file = 'time_mah_3_500_n.RData')
save(ari_mah_3_250_n, file = 'ari_mah_3_250_n.RData')
save(jaccard_mah_3_250_n, file = 'jaccard_mah_3_250_n.RData')
save(larsen_mah_3_250_n, file = 'larsen_mah_3_250_n.RData')
save(loo_mah_3_250_n, file = 'loo_mah_3_250_n.RData')
save(clustering_mah_3_250_n, file = 'clustering_mah_3_250_n.RData')
save(time_mah_3_250_n, file = 'time_mah_3_250_n.RData')

# Constructing the tables of results

files_rdata <- list.files(pattern = "*.RData")
l_rdata <- length(files_rdata)

for (i in 1:l_rdata) {
  
  load(files_rdata[i])
  
}

table_250 <- matrix(nrow = 8, ncol = 6)
table_250[1, 1] <- round(mean(ari_qc2_3_250_n), 3)
table_250[1, 2] <- round(mean(ari_w_3_250_n), 3)
table_250[1, 3] <- round(mean(ari_gcc_3_250_n), 3)
table_250[1, 4] <- round(mean(ari_kst_3_250_n), 3)
table_250[1, 5] <- round(mean(ari_mah_3_250_n), 3)
table_250[1, 6] <- round(mean(ari_pca_3_250_n), 3)
table_250[2, 1] <- round(sd(ari_qc2_3_250_n), 3)
table_250[2, 2] <- round(sd(ari_w_3_250_n), 3)
table_250[2, 3] <- round(sd(ari_gcc_3_250_n), 3)
table_250[2, 4] <- round(sd(ari_kst_3_250_n), 3)
table_250[2, 5] <- round(sd(ari_mah_3_250_n), 3)
table_250[2, 6] <- round(sd(ari_pca_3_250_n), 3)
table_250[3, 1] <- round(mean(larsen_qc2_3_250_n), 3)
table_250[3, 2] <- round(mean(larsen_w_3_250_n), 3)
table_250[3, 3] <- round(mean(larsen_gcc_3_250_n), 3)
table_250[3, 4] <- round(mean(larsen_kst_3_250_n), 3)
table_250[3, 5] <- round(mean(larsen_mah_3_250_n), 3)
table_250[3, 6] <- round(mean(larsen_pca_3_250_n), 3)
table_250[4, 1] <- round(sd(larsen_qc2_3_250_n), 3)
table_250[4, 2] <- round(sd(larsen_w_3_250_n), 3)
table_250[4, 3] <- round(sd(larsen_gcc_3_250_n), 3)
table_250[4, 4] <- round(sd(larsen_kst_3_250_n), 3)
table_250[4, 5] <- round(sd(larsen_mah_3_250_n), 3)
table_250[4, 6] <- round(sd(larsen_pca_3_250_n), 3)
table_250[5, 1] <- round(mean(loo_qc2_3_250_n), 3)
table_250[5, 2] <- round(mean(loo_w_3_250_n), 3)
table_250[5, 3] <- round(mean(loo_gcc_3_250_n), 3)
table_250[5, 4] <- round(mean(loo_kst_3_250_n), 3)
table_250[5, 5] <- round(mean(loo_mah_3_250_n), 3)
table_250[5, 6] <- round(mean(loo_pca_3_250_n), 3)
table_250[6, 1] <- round(sd(loo_qc2_3_250_n), 3)
table_250[6, 2] <- round(sd(loo_w_3_250_n), 3)
table_250[6, 3] <- round(sd(loo_gcc_3_250_n), 3)
table_250[6, 4] <- round(sd(loo_kst_3_250_n), 3)
table_250[6, 5] <- round(sd(loo_mah_3_250_n), 3)
table_250[6, 6] <- round(sd(loo_pca_3_250_n), 3)
table_250[7, 1] <- round(mean(jaccard_qc2_3_250_n), 3)
table_250[7, 2] <- round(mean(jaccard_w_3_250_n), 3)
table_250[7, 3] <- round(mean(jaccard_gcc_3_250_n), 3)
table_250[7, 4] <- round(mean(jaccard_kst_3_250_n), 3)
table_250[7, 5] <- round(mean(jaccard_mah_3_250_n), 3)
table_250[7, 6] <- round(mean(jaccard_pca_3_250_n), 3)
table_250[8, 1] <- round(sd(jaccard_qc2_3_250_n), 3)
table_250[8, 2] <- round(sd(jaccard_w_3_250_n), 3)
table_250[8, 3] <- round(sd(jaccard_gcc_3_250_n), 3)
table_250[8, 4] <- round(sd(jaccard_kst_3_250_n), 3)
table_250[8, 5] <- round(sd(jaccard_mah_3_250_n), 3)
table_250[8, 6] <- round(sd(jaccard_pca_3_250_n), 3)



table_1000 <- matrix(nrow = 8, ncol = 6)
table_1000[1, 1] <- round(mean(ari_qc2_3_1000_n), 3)
table_1000[1, 2] <- round(mean(ari_w_3_1000_n), 3)
table_1000[1, 3] <- round(mean(ari_gcc_3_1000_n), 3)
table_1000[1, 4] <- round(mean(ari_kst_3_1000_n), 3)
table_1000[1, 5] <- round(mean(ari_mah_3_1000_n), 3)
table_1000[1, 6] <- round(mean(ari_pca_3_1000_n), 3)
table_1000[2, 1] <- round(sd(ari_qc2_3_1000_n), 3)
table_1000[2, 2] <- round(sd(ari_w_3_1000_n), 3)
table_1000[2, 3] <- round(sd(ari_gcc_3_1000_n), 3)
table_1000[2, 4] <- round(sd(ari_kst_3_1000_n), 3)
table_1000[2, 5] <- round(sd(ari_mah_3_1000_n), 3)
table_1000[2, 6] <- round(sd(ari_pca_3_1000_n), 3)
table_1000[3, 1] <- round(mean(larsen_qc2_3_1000_n), 3)
table_1000[3, 2] <- round(mean(larsen_w_3_1000_n), 3)
table_1000[3, 3] <- round(mean(larsen_gcc_3_1000_n), 3)
table_1000[3, 4] <- round(mean(larsen_kst_3_1000_n), 3)
table_1000[3, 5] <- round(mean(larsen_mah_3_1000_n), 3)
table_1000[3, 6] <- round(mean(larsen_pca_3_1000_n), 3)
table_1000[4, 1] <- round(sd(larsen_qc2_3_1000_n), 3)
table_1000[4, 2] <- round(sd(larsen_w_3_1000_n), 3)
table_1000[4, 3] <- round(sd(larsen_gcc_3_1000_n), 3)
table_1000[4, 4] <- round(sd(larsen_kst_3_1000_n), 3)
table_1000[4, 5] <- round(sd(larsen_mah_3_1000_n), 3)
table_1000[4, 6] <- round(sd(larsen_pca_3_1000_n), 3)
table_1000[5, 1] <- round(mean(loo_qc2_3_1000_n), 3)
table_1000[5, 2] <- round(mean(loo_w_3_1000_n), 3)
table_1000[5, 3] <- round(mean(loo_gcc_3_1000_n), 3)
table_1000[5, 4] <- round(mean(loo_kst_3_1000_n), 3)
table_1000[5, 5] <- round(mean(loo_mah_3_1000_n), 3)
table_1000[5, 6] <- round(mean(loo_pca_3_1000_n), 3)
table_1000[6, 1] <- round(sd(loo_qc2_3_1000_n), 3)
table_1000[6, 2] <- round(sd(loo_w_3_1000_n), 3)
table_1000[6, 3] <- round(sd(loo_gcc_3_1000_n), 3)
table_1000[6, 4] <- round(sd(loo_kst_3_1000_n), 3)
table_1000[6, 5] <- round(sd(loo_mah_3_1000_n), 3)
table_1000[6, 6] <- round(sd(loo_pca_3_1000_n), 3)
table_1000[7, 1] <- round(mean(jaccard_qc2_3_1000_n), 3)
table_1000[7, 2] <- round(mean(jaccard_w_3_1000_n), 3)
table_1000[7, 3] <- round(mean(jaccard_gcc_3_1000_n), 3)
table_1000[7, 4] <- round(mean(jaccard_kst_3_1000_n), 3)
table_1000[7, 5] <- round(mean(jaccard_mah_3_1000_n), 3)
table_1000[7, 6] <- round(mean(jaccard_pca_3_1000_n), 3)
table_1000[8, 1] <- round(sd(jaccard_qc2_3_1000_n), 3)
table_1000[8, 2] <- round(sd(jaccard_w_3_1000_n), 3)
table_1000[8, 3] <- round(sd(jaccard_gcc_3_1000_n), 3)
table_1000[8, 4] <- round(sd(jaccard_kst_3_1000_n), 3)
table_1000[8, 5] <- round(sd(jaccard_mah_3_1000_n), 3)
table_1000[8, 6] <- round(sd(jaccard_pca_3_1000_n), 3)

