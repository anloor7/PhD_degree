

# Simulating an Scenario of NDARMA processes 

K <- 4

# Data for Cluster 1

marg_prob_1 <- c(0.8, 0.2, 0.2)
probs_multin_1 <- c(0.7, 0.2, 0.1)
p_1 <- 2
q_1 <- 0

# Data for Cluster 2

marg_prob_2 <- c(0.1, 0.1, 0.8)
probs_multin_2 <- c(0.1, 0.45, 0.45)
p_2 <- 2
q_2 <- 0


# Data for Cluster 3

marg_prob_3 <- c(0.25, 0.5, 0.25)
probs_multin_3 <- c(0.5, 0.25, 0.25)
p_3 <- 2
q_3 <- 0


# Data for Cluster 4

marg_prob_4 <- c(0.4, 0.4, 0.2)
probs_multin_4 <- c(0.2, 0.8)
p_4 <- 1
q_4 <- 0

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
cluster4 <- list()

# Simulating 5 series from each cluster

vector_l <- c(600)
m <- 2 # Fuzziness coefficient


set.seed(321)





ari_cramer <- numeric()
ari_indicators <- numeric()
ari_cramer_cohen <- numeric()
ari_total <- numeric()
ari_markov <- numeric()
ari_additional <- numeric()
ari_cadez <- numeric()
ari_magariños_vilar <- numeric()
ari_mv <- numeric()
ari_weib <- numeric()
jaccard_cramer <- numeric()
jaccard_indicators <- numeric()
jaccard_cramer_cohen <- numeric()
jaccard_total <- numeric()
jaccard_markov <- numeric()
jaccard_additional <- numeric()
jaccard_cadez <- numeric()
jaccard_magariños_vilar <- numeric()
jaccard_mv <- numeric()
jaccard_weib <- numeric()
larsen_cramer <- numeric()
larsen_indicators <- numeric()
larsen_cramer_cohen <- numeric()
larsen_total <- numeric()
larsen_markov <- numeric()
larsen_additional <- numeric()
larsen_cadez <- numeric()
larsen_magariños_vilar <- numeric()
larsen_mv <- numeric()
larsen_weib <- numeric()

time_cramer <- list()
time_indicators <- list()
time_weib <- list()
time_cadez <- list()
time_mv <- list()

for (l in vector_l) {
  
  for (i in 1 : 500) {
    
    
    for (p in 1 : 5) {
      
      cluster1[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_1, p_1, q_1, series_length = l)
      cluster2[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_2, p_2, q_2, series_length = l)
      cluster3[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_3, p_3, q_3, series_length = l)
      cluster4[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_4, p_4, q_4, series_length = l)
      
      
    }
    
    # Converting the series in factors
    
    cluster1 <- lapply(cluster1, factor, levels = c('0', '1', '2'))
    cluster2 <- lapply(cluster2, factor, levels = c('0', '1', '2'))
    cluster3 <- lapply(cluster3, factor, levels = c('0', '1', '2'))
    cluster4 <- lapply(cluster4, factor, levels = c('0', '1', '2'))
    cluster <- c(cluster1, cluster2, cluster3, cluster4)
    n_cluster1 <- lapply(cluster1, as.numeric)
    n_cluster2 <- lapply(cluster2, as.numeric)
    n_cluster3 <- lapply(cluster3, as.numeric)
    n_cluster4 <- lapply(cluster4, as.numeric)
    n_cluster <- c(n_cluster1, n_cluster2, n_cluster3, n_cluster4)
    ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5))
    
    
    
    
    
    # Algorithms
    
    # cramer
    
    categories <- factor(c('0', '1', '2'))
    clusterExport(c1, c('cluster', 'p_i_j_k_function'))
    
    
    
    # clustering_cramer <- pam(features_cramer, k = K)
    # fuzzy_matrix_cramer <- clustering_cramer$clustering
    # ari_cramer[i] <- 
    # external_validation(ground_truth, fuzzy_matrix_cramer)
    # jaccard_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer,
    #                                        method = 'jaccard_index')
    # larsen_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer,
    #                                      method = 'fowlkes_mallows_index')
    
    
    
    # Indicators
    
    
    initial_time <- Sys.time()
    features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
    features_indicators <- listTomatrix(parLapply(c1, cluster, features_categorical_series, 
                                                  levels_dataset = categories, l = 2))
    
    clustering_indicators <- pam((cbind(features_indicators, features_marginal)), k = K)
    fuzzy_matrix_indicators <- clustering_indicators$clustering
    final_time <- Sys.time()
    time_indicators[[i]] <- final_time - initial_time
    ari_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators)
    jaccard_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                                 method = 'jaccard_index')
    larsen_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                                method = 'fowlkes_mallows_index')
    
    # Association measures
    
    initial_time <- Sys.time()
    features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
    features_cramer <- listTomatrix(parLapply(c1, cluster, features_cramer_v_function, categories = categories,
                                               max_lag = 2))
    features_cohen <- listTomatrix(parLapply(c1, cluster, features_cohen_function, categories = categories,
                                             max_lag = 2))
    
    # clustering_additional <- pam(scale(cbind(features_cramer, features_cohen, features_marginal, features_indicators)), k = K)
    # fuzzy_matrix_additional <- clustering_additional$clustering
    # ari_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional)
    # jaccard_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional,
    #                                            method = 'jaccard_index')
    
    
    
    
    # cramer and Cohens
    
    # categories <- factor(c('a', 'b', 'c'))
    # clusterExport(c1, c('cluster', 'p_i_j_k_function'))
    
    # clustering_cramer_cohen <- pam(scale(cbind(features_cramer, features_cohen)), k = K)
    # fuzzy_matrix_cramer_cohen <- clustering_cramer_cohen$clustering
    # ari_cramer_cohen[i] <- external_validation(ground_truth, fuzzy_matrix_cramer_cohen)
    
    clustering_total <- pam((cbind(features_cramer, features_cohen, features_marginal)), k = K)
    final_time <- Sys.time()
    time_cramer[[i]] <- final_time - initial_time
    fuzzy_matrix_total <- clustering_total$clustering
    
    ari_total[i] <- external_validation(ground_truth, fuzzy_matrix_total)
    jaccard_total[i] <- external_validation(ground_truth, fuzzy_matrix_total, 
                                            method = 'jaccard_index')
    larsen_total[i] <- external_validation(ground_truth, fuzzy_matrix_total,
                                           method = 'fowlkes_mallows_index')
    
    # MLE 
    
    # Binarizing all the series into matrixes
    
    initial_time <- Sys.time()
    list_matrices <- list()
    
    for (h in 1 : 20) {
      
      list_matrices[[h]] <- bincodes[cluster[[h]],][, (1 : 3)]
      
    }
    
    matrix_features <- matrix(0, nrow = 20, ncol = 6)
    
    for (h in 1 : 20) {
      
      auxiliary_computation <- suppressWarnings(constrOptim(c(1/4, 1/4, 1/4, 1/4), ll_darp, NULL, ui=Amat(3,2), ci=bvec(3,2), databin = list_matrices[[h]]))
      matrix_features[h,] <- c(auxiliary_computation$par, 1-sum(auxiliary_computation$par[1 : 2]), 
                               1-sum(auxiliary_computation$par[3 : 4]))
      
    }
    
    clustering_weib <- pam(matrix_features[,c(1 : 4)], k = K)
    fuzzy_matrix_weib <- clustering_weib$clustering
    final_time <- Sys.time()
    time_weib[[i]] <- final_time - initial_time
    ari_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib)
    jaccard_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib, 
                                           method = 'jaccard_index')
    larsen_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib,
                                          method = 'fowlkes_mallows_index')
    
    # Cadez
    
    initial_time <- Sys.time()
    data_cadez <- click.read(n_cluster)
    cadez <- click.EM(X = data_cadez$X, K = K, r = 1500)
    final_time <- Sys.time()
    time_cadez[[i]] <- final_time - initial_time
    ari_cadez[i] <- external_validation(ground_truth, cadez$id)
    jaccard_cadez[i] <- external_validation(ground_truth, cadez$id,
                                            method = 'jaccard_index')
    larsen_cadez[i] <- external_validation(ground_truth, cadez$id,
                                           method = 'fowlkes_mallows_index')
    
    # Magariños-Vilar
    
    # initial_time <- Sys.time()
    # cluster_matrix <- matrix(NA, nrow = 20, ncol = l)
    
    # for (p in 1 : 20) {
      
      # cluster_matrix[p,] <- as.numeric(cluster[[p]])
      
    # }
    
    # magariños_vilar <- Diss.Mat(cluster_matrix, opt = "AD-Corr", kk = 1)
    # pam_magariños_vilar <- pam(magariños_vilar$Dist, k = K)
    # fuzzy_matrix_mv <- pam_magariños_vilar$clustering
    # final_time <- Sys.time()
    # time_mv[[i]] <- final_time - initial_time
    # ari_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv)
    # jaccard_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv,
    #                                     method = 'jaccard_index')
    # larsen_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv,
    #                                    method = 'fowlkes_mallows_index')
    
    
    
    print(i)
    
  }
  
  if (l == 200) {
    
    save(ari_markov, file = 'ari_markov_200.RData')
    save(ari_total, file = 'ari_total_200.RData')
    save(ari_indicators, file = 'ari_indicators_200.RData')
    save(ari_cadez, file = 'ari_cadez_200.RData')
    save(ari_mv, file = 'ari_mv_200.RData')
    save(ari_weib, file = 'ari_weib_200.RData')
    save(jaccard_markov, file = 'jaccard_markov_200.RData')
    save(jaccard_total, file = 'jaccard_total_200.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_200.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_200.RData')
    save(jaccard_mv, file = 'jaccard_mv_200.RData')
    save(jaccard_weib, file = 'jaccard_weib_200.RData')
    save(larsen_markov, file = 'larsen_markov_200.RData')
    save(larsen_total, file = 'larsen_total_200.RData')
    save(larsen_indicators, file = 'larsen_indicators_200.RData')
    save(larsen_cadez, file = 'larsen_cadez_200.RData')
    save(larsen_mv, file = 'larsen_mv_200.RData')
    save(larsen_weib, file = 'larsen_weib_200.RData')
    
  }
  
  
  if (l == 600) {
    
    save(ari_markov, file = 'ari_markov_600.RData')
    save(ari_total, file = 'ari_total_600.RData')
    save(ari_indicators, file = 'ari_indicators_600.RData')
    save(ari_cadez, file = 'ari_cadez_600.RData')
    save(ari_mv, file = 'ari_mv_600.RData')
    save(ari_weib, file = 'ari_weib_600.RData')
    save(jaccard_markov, file = 'jaccard_markov_600.RData')
    save(jaccard_total, file = 'jaccard_total_600.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_600.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_600.RData')
    save(jaccard_mv, file = 'jaccard_mv_600.RData')
    save(jaccard_weib, file = 'jaccard_weib_600.RData')
    save(larsen_markov, file = 'larsen_markov_600.RData')
    save(larsen_total, file = 'larsen_total_600.RData')
    save(larsen_indicators, file = 'larsen_indicators_600.RData')
    save(larsen_cadez, file = 'larsen_cadez_600.RData')
    save(larsen_mv, file = 'larsen_mv_600.RData')
    save(larsen_weib, file = 'larsen_weib_600.RData')
    
  }
  
}



# Creating the table

setwd('/Users/angel/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/crisp_algorithms/scenario_3')





table_previous <- matrix(NA, nrow = 5, ncol = 6)
load('ari_markov_200.RData')
load('ari_total_200.RData')
load('ari_indicators_200.RData')
load('ari_cadez_200.RData')
load('ari_mv_200.RData')
load('ari_weib_200.RData')
load('jaccard_markov_200.RData')
load('jaccard_total_200.RData')
load('jaccard_indicators_200.RData')
load('jaccard_cadez_200.RData')
load('jaccard_mv_200.RData')
load('jaccard_weib_200.RData')
load('larsen_mv_200.RData')
load('larsen_markov_200.RData')
load('larsen_total_200.RData')
load('larsen_indicators_200.RData')
load('larsen_cadez_200.RData')
load('larsen_mv_200.RData')
load('larsen_weib_200.RData')
table_previous[1, 1] <- mean(ari_total)
table_previous[1, 2] <- mean(jaccard_total)
table_previous[1, 3] <- mean(larsen_total)
table_previous[2, 1] <- mean(ari_indicators)
table_previous[2, 2] <- mean(jaccard_indicators)
table_previous[2, 3] <- mean(larsen_indicators)
table_previous[3, 1] <- mean(ari_weib)
table_previous[3, 2] <- mean(jaccard_weib)
table_previous[3, 3] <- mean(larsen_weib)
table_previous[4, 1] <- mean(ari_cadez)
table_previous[4, 2] <- mean(jaccard_cadez)
table_previous[4, 3] <- mean(larsen_cadez)
table_previous[5, 1] <- mean(ari_mv)
table_previous[5, 2] <- mean(jaccard_mv)
table_previous[5, 3] <- mean(larsen_mv)




load('ari_markov_600.RData')
load('ari_total_600.RData')
load('ari_indicators_600.RData')
load('ari_cadez_600.RData')
load('ari_mv_600.RData')
load('ari_weib_600.RData')
load('jaccard_markov_600.RData')
load('jaccard_total_600.RData')
load('jaccard_indicators_600.RData')
load('jaccard_cadez_600.RData')
load('jaccard_mv_600.RData')
load('jaccard_weib_600.RData')
load('larsen_markov_600.RData')
load('larsen_total_600.RData')
load('larsen_indicators_600.RData')
load('larsen_cadez_600.RData')
load('larsen_mv_600.RData')
load('larsen_weib_600.RData')
table_previous[1, 4] <- mean(ari_total)
table_previous[1, 5] <- mean(jaccard_total)
table_previous[1, 6] <- mean(larsen_total)
table_previous[2, 4] <- mean(ari_indicators)
table_previous[2, 5] <- mean(jaccard_indicators)
table_previous[2, 6] <- mean(larsen_indicators)
table_previous[3, 4] <- mean(ari_weib)
table_previous[3, 5] <- mean(jaccard_weib)
table_previous[3, 6] <- mean(larsen_weib)
table_previous[4, 4] <- mean(ari_cadez)
table_previous[4, 5] <- mean(jaccard_cadez)
table_previous[4, 6] <- mean(larsen_cadez)
table_previous[5, 4] <- mean(ari_mv)
table_previous[5, 5] <- mean(jaccard_mv)
table_previous[5, 6] <- mean(larsen_mv)


xtable(table_previous)
