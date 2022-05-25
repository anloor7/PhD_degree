

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_6')

K <- 2

# Data for Cluster 1

marg_prob_1 <- c(1/3, 1/3, 1/3)
probs_multin_1 <- c(0.7, 0.15, 0.15)
p_1 <- 2
q_1 <- 0

# Data for Cluster 2

marg_prob_2 <- c(1/3, 1/3, 1/3)
probs_multin_2 <- c(0.1, 0.45, 0.45)
p_2 <- 2
q_2 <- 0


# Data for Cluster 3

marg_prob_3 <- c(1/3, 1/3, 1/3)
probs_multin_3 <- (probs_multin_1 + probs_multin_2)/2
p_3 <- 2
q_3 <- 0

cutoff <- 0.7
l <- c(200)
vector_m <- c(1.5, 1.8, 2, 2.2)
w_last <- 500
m <- 2 # Fuzziness coefficient

set.seed(123467)

rate_am <- numeric()
rate_am_m <- numeric()
rate_ip <- numeric()
rate_ip_m <- numeric()
rate_mle <- numeric()
rate_mle_m <- numeric()
rate_cadez <- numeric()
rate_cadez_m <- numeric()
cluster <- list()
n_cluster <- list()
cluster_df <- list()

cluster1 <- list()
cluster2 <- list()
cluster3 <- list()
  
  for (w in 1 : w_last) {
    
    
    for (p in 1 : 5) {
      
      cluster1[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_1, p_1, q_1, series_length = l)
      cluster2[[p]] <- ndarma_sim_function(marg_prob_1, probs_multin_2, p_2, q_2, series_length = l)
      cluster3[[1]] <- ndarma_sim_function(marg_prob_1, probs_multin_3, p_3, q_3, series_length = l)
      
      
    }
    
    # Converting the series in factors
    
    cluster1 <- lapply(cluster1, factor, levels = c('0', '1', '2'))
    cluster2 <- lapply(cluster2, factor, levels = c('0', '1', '2'))
    cluster3 <- lapply(cluster3, factor, levels = c('0', '1', '2'))
    cluster[[w]] <- c(cluster1, cluster2, cluster3)
    n_cluster1 <- lapply(cluster1, as.numeric)
    n_cluster2 <- lapply(cluster2, as.numeric)
    n_cluster3 <- lapply(cluster3, as.numeric)
    n_cluster[[w]] <- c(n_cluster1, n_cluster2, n_cluster3)
    ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5))
    
    print(w)
    
    
  }


# Algorithms





# Indicators

categories <- factor(c('0', '1', '2'))
clusterExport(c1, c('cluster', 'p_i_j_k_function'))
set.seed(321)

k = 1
for (m in c(vector_m)) {
  
  
  
  for (w in (1 : w_last)) {
    
    features_marginal <- listTomatrix(parLapply(c1, cluster[[w]], marginal_probabilities, categories = categories))
    features_indicators <- listTomatrix(parLapply(c1, cluster[[w]], features_categorical_series, l = 2,
                                                  levels_dataset = categories))
    
    clustering_indicators <- FKM(X = cbind(features_indicators, features_marginal), k = K, m = m)
    fuzzy_matrix_indicators <- clustering_indicators$U
    fuzzy_matrix_indicators
    
    clustering_solution <- clustering_indicators$clus
    l_1 <- length(unique(clustering_solution[1 : 5]))
    l_2 <- length(unique(clustering_solution[6 : 10]))
    l_3 <- length(unique(clustering_solution[1 : 10]))
    
    if (l_1 == 1 & l_2 == 1 & l_3 == 2 & 
        sum(clustering_indicators$clus[,2][1 : 10] > cutoff) == 10 &
        clustering_indicators$clus[,2][11] < cutoff) {
      
      rate_ip[w] <- 1
      
    } else {
      
      rate_ip[w] <- 0
      
    }
    
    
    print(w)
  }
  
  rate_ip_m[k] = mean(rate_ip)
  k = k + 1
  print(m)
  
}



# Total

categories <- factor(c('0', '1', '2'))
clusterExport(c1, c('cluster', 'p_i_j_k_function'))

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : w_last)) {
    
    features_kruskal <- listTomatrix(parLapply(c1, cluster[[w]], features_kruskal_tau_function, categories = categories,
                                               max_lag = 2))
    features_marginal <- listTomatrix(parLapply(c1, cluster[[w]], marginal_probabilities, categories = categories))
    features_cohen <- listTomatrix(parLapply(c1, cluster[[w]], features_cohen_function, categories = categories,
                                             max_lag = 2))
    clustering_total <- FKM((cbind(features_kruskal, features_cohen, features_marginal)), k = K, m = m)
    
    fuzzy_matrix_total <- clustering_total$U
    
    clustering_solution <- clustering_total$clus
    l_1 <- length(unique(clustering_solution[1 : 5]))
    l_2 <- length(unique(clustering_solution[6 : 10]))
    l_3 <- length(unique(clustering_solution[1 : 10]))
    
    if (l_1 == 1 & l_2 == 1 & l_3 == 2 & 
        sum(clustering_total$clus[,2][1 : 10] > cutoff) == 10 &
        clustering_total$clus[,2][11] < cutoff) {
      
      rate_am[w] <- 1
      
    } else {
      
      rate_am[w] <- 0
      
    }
    
    
    print(w)
  }
  
  rate_am_m[k] = mean(rate_am)
  k = k + 1
  print(m)
  
}


# MLE

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : w_last)) {
    
    # MLE 
    
    # Binarizing all the series into matrixes
    
    list_matrices <- list()
    
    for (h in 1 : 11) {
      
      list_matrices[[h]] <- bincodes[cluster[[w]][[h]],][, (1 : 3)]
      
    }
    
    matrix_features <- matrix(0, nrow = 11, ncol = 6)
    
    for (h in 1 : 11) {
      
      auxiliary_computation <- suppressWarnings(constrOptim(c(1/4, 1/4, 1/4, 1/4), ll_darp, NULL, ui=Amat(3,2), ci=bvec(3,2), databin = list_matrices[[h]]))
      matrix_features[h,] <- c(auxiliary_computation$par, 1-sum(auxiliary_computation$par[1 : 2]), 
                               1-sum(auxiliary_computation$par[3 : 4]))
      
    }
    
    clustering_mle <- FKM(X = matrix_features[,c(1 : 4)], k = K, m = m)
  
    fuzzy_matrix_mle <- clustering_mle$U
    fuzzy_matrix_mle
    clustering_solution <- clustering_mle$clus
    l_1 <- length(unique(clustering_solution[1 : 5]))
    l_2 <- length(unique(clustering_solution[6 : 10]))
    l_3 <- length(unique(clustering_solution[1 : 10]))
    
    if (l_1 == 1 & l_2 == 1 & l_3 == 2 & 
        sum(clustering_mle$clus[,2][1 : 10] > cutoff) == 10 &
        clustering_mle$clus[,2][11] < cutoff) {
      
      rate_mle[w] <- 1
      
    } else {
      
      rate_mle[w] <- 0
      
    }
    
    
    print(w)
  }
  
  rate_mle_m[k] = mean(rate_mle)
  k = k + 1
  print(m)
  
}



# Cadez

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : w_last)) {
    
    data_cadez <- click.read(n_cluster[[w]])
    cadez <- click.EM(X = data_cadez$X, K = K, r = 10)
    
    clustering_cadez <- cadez$id
    fuzzy_matrix_cadez <- cadez$z
    fuzzy_matrix_indicators
    fuzzy_matrix_cadez_max <- apply(fuzzy_matrix_cadez, 1, max)
    
    clustering_solution <- clustering_cadez
    l_1 <- length(unique(clustering_solution[1 : 5]))
    l_2 <- length(unique(clustering_solution[6 : 10]))
    l_3 <- length(unique(clustering_solution[1 : 10]))
    
    if (l_1 == 1 & l_2 == 1 & l_3 == 2 & 
        sum(fuzzy_matrix_cadez_max[1 : 10] > cutoff) == 10 &
        fuzzy_matrix_cadez_max[11] < cutoff) {
      
      rate_cadez[w] <- 1
      
    } else {
      
      rate_cadez[w] <- 0
      
    }
    
    
    
  }
  
  rate_cadez_m[k] = mean(rate_cadez)
  k = k + 1
  print(m)
  
}




# Creating the table

table_l_200 <- matrix(0, nrow = 3, ncol = 4)
table_l_200[1, (1 : 4)] <- rate_am_m   
table_l_200[2, (1 : 4)] <- rate_ip_m   
table_l_200[3, (1 : 4)] <- rate_mle_m

save(table_l_200, file = 'table_l_200.RData')   

# Creating the whole table

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_6')
load('table_l_200.RData')
load('table_l_600.RData')

xtable(rbind(table_l_200, table_l_600))     
