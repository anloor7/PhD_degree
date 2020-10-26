
n <- 5  # Number of series per cluster
K <- 5 # Number of clusters
B <- 100 # Number of Monte Carlos trials
ground_truth <- c(rep(1, n), rep(2, n), rep(3, n), rep(4, n), rep(5, n)) # Ground Truth

load('cluster_1_250_n.RData')
load('cluster_1_500_n.RData')
load('cluster_1_1000_n.RData')


cluster <- loadRData(paste0('cluster_1_', 250, '_n.RData'))
cluster <- cluster[[15]]



coherence2_quantities <- listTomatrix(lapply(cluster, quantile_quantities_re_im))
dis_matrix_quantities <- proxy::dist(coherence2_quantities, EuclideanDistance) 
pamc <- pam(dis_matrix_quantities, K)
clustering <- pamc$clustering
external_validation(ground_truth, clustering)

coherence2_coherence <- listTomatrix(lapply(cluster, quantile_coherence_re_im))
dis_matrix_coherence <- proxy::dist(coherence2_coherence, EuclideanDistance) 
pamc <- pam(dis_matrix_coherence, K)
clustering <- pamc$clustering
external_validation(ground_truth, clustering)


dis_matrix <- matrix(0, 25, 25)
clustercov <- foreach (m = 1 : length(cluster)) %dopar% {
  cov(cluster[[m]])
}

dis_matrix <- 1 - PCAsimilarity(clustercov)
dis_matrix[col(dis_matrix) >= row(dis_matrix)] <- 0
dis_matrix_pca <- as.dist(dis_matrix)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering, summary_stats = T)






pamc <- pam(0.6*dis_matrix_quantities/max(dis_matrix_quantities) + 0.4*dis_matrix_pca/max(dis_matrix_pca), K)
clustering <- pamc$clustering
external_validation(ground_truth, clustering)



features <- listTomatrix(lapply(cluster, gcc_features_mts))
dis_matrix_gcc <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, K)$cluster
external_validation(ground_truth, clustering)

pamc <- pam(0.90*dis_matrix_quantities/max(dis_matrix_quantities) + 0.1*dis_matrix_gcc/max(dis_matrix_gcc), K)
clustering <- pamc$clustering
external_validation(ground_truth, clustering)


coherence_t <- cbind(coherence2_quantities, coherence2_coherence)
dis_matrix_coherence_t <- proxy::dist(coherence_t, EuclideanDistance)
dis_matrix_coherence_t_normalized <- dis_matrix_coherence_t/max(dis_matrix_coherence_t)
pamc <- pam(dis_matrix_coherence_t_normalized, K)
clustering <- pamc$clustering
external_validation(ground_truth, clustering)

