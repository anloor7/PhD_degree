

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_5')


# Simulating an Scenario of Markov processes of order 1

E <- c('a', 'b', 'c') # State space
s_e <- length(E)
K <- 2


# Data for Cluster 1

p1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
e1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
d1 <- round(matrix.power(p1, 2000)[1,], 3)


# Data for Cluster 2

p2 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
e2 <- matrix(c(0.1, 0.4, 0.5, 0.25, 0.25, 0.5, 0.2, 0.5, 0.3), nrow = 3, byrow = T)
d2 <- round(matrix.power(p2, 2000)[1,], 3)


# Data for Cluster 3

p3 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
e3 <- (e1 + e2)/2
d3 <- round(matrix.power(p3, 2000)[1,], 3)


cutoff <- 0.7
l <- c(600)
vector_m <- c(1.5, 1.8, 2, 2.2)
w_last <- 500
m <- 2 # Fuzziness coefficient

TP <- array(rep(NA, 3 * 3 * 4), c(3, 3, 3)) 
TP[,,1] <- p1
TP[,,2] <- p2
TP[,,3] <- p3

# Simulating 5 series from each cluster



set.seed(1234)

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

for (w in 1 : w_last) {
  
  cluster1 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = p1, emission_probs = e1, sequence_length = l)
  cluster2 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = p2, emission_probs = e2, sequence_length = l)
  cluster3 <- simulate_hmm(n_sequences = 1, initial_probs = d1, transition_probs = p3, emission_probs = e3, sequence_length = l)
  
  
  # Converting the series in factors
  
  
  cluster1 <- apply(cluster1$observations, 2, factor)
  cluster2 <- apply(cluster2$observations, 2, factor)
  cluster3 <- apply(cluster3$observations, 2, factor)
  cluster[[w]] <- rbind(cluster1, cluster2, cluster3)
  n_cluster1 <- apply(cluster1, 2, as.numeric)
  n_cluster2 <- apply(cluster2, 2, as.numeric)
  n_cluster3 <- cluster3
  n_cluster[[w]] <- rbind(n_cluster1, n_cluster2, n_cluster3)
  
  cluster_df[[w]] <- matrix(0, nrow = nrow(cluster[[w]]), ncol = ncol(cluster[[w]]))
  cluster_df[[w]] <- as.data.frame(cluster_df[[w]])
  
  for (q in 1 : (11)) {
    
    cluster_df[[w]][q,] <- factor(cluster[[w]][q,])
    
  }
  
  
  print(w)
  
}

save(cluster, file = 'cluster_l_200.RData')
save(cluster_df, file = 'cluster_df_l_200.RData')

# Algorithms

# Total

categories <- factor(c('1', '2', '3'))
clusterExport(c1, c('cluster', 'p_i_j_k_function'))

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : w_last)) {
    
    features_cramer <- t(parApply(c1, cluster_df[[w]], 1, features_cramer_v_function, categories = categories,
                                   max_lag = 1))
    features_marginal <- t(parApply(c1, cluster_df[[w]], 1, marginal_probabilities, categories = categories))
    features_cohen <- t(parApply(c1, cluster_df[[w]], 1, features_cohen_function, categories = categories,
                                 max_lag = 1))
    clustering_total <- FKM((cbind(features_cramer, features_cohen, features_marginal)), k = K, m = m)
    
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
    
    
    
  }
  
  rate_am_m[k] = mean(rate_am)
  k = k + 1
  print(m)
  
}

save(cluster, file = '')


# Indicators

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : w_last)) {
    
    features_marginal <- features_marginal <- t(parApply(c1, cluster_df[[w]], 1, marginal_probabilities, categories = categories))
    features_indicators <- t(parApply(c1, cluster_df[[w]], 1, features_categorical_series, l = 1, levels_dataset = categories))
    
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
    
    
    
  }
  
  rate_ip_m[k] = mean(rate_ip)
  k = k + 1
  print(m)
  
}


# MLE

set.seed(321)

k = 1
for (m in vector_m) {
  
  
  
  for (w in (1 : 500)) {
    
    features_mle <- matrix(0, nrow = 11, ncol = 18)
    
    for (v in 1 :11){
      
      features_mle[v,] <- hmm_fit(cluster_df[[w]][v,])
      
    }
    
    clustering_mle <- pam(features_mle, k = K)
    
    clustering_mle <- FKM(X = features_mle, k = K, m = m)
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

table_l_600 <- matrix(0, nrow = 3, ncol = 4)
table_l_600[1, (1 : 4)] <- rate_am_m   
table_l_600[2, (1 : 4)] <- rate_ip_m   
table_l_600[3, (1 : 4)] <- rate_mle_m

save(table_l_600, file = 'table_l_600.RData')   

# Creating the whole table

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_5')
load('table_l_200.RData')
load('table_l_600.RData')

xtable(rbind(table_l_200, table_l_600))    
