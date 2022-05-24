
library(SMM)
library(markovchain)
library(cluster)
library(ClusterR)
library(ClickClust)
library(seqHMM)
library(matrixcalc)
library(parallel)
c1 <- makeCluster(7)
setwd('/Users/angel/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/crisp_algorithms/scenario_2')


# set.seed(1234)


E <- c('a', 'b', 'c') # State space
s_e <- length(E)
K <- 4
n <- 5

# Data for Cluster 1

p1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
e1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
d1 <- round(matrix.power(p1, 2000)[1,], 3)


# Data for Cluster 2

p2 <- matrix(c(0.1, 0.5, 0.4, 0.5, 0.4, 0.1, 0.3, 0.3, 0.4), nrow = 3, byrow = T)
e2 <- matrix(c(0.1, 0.7, 0.2, 0.4, 0.4, 0.2, 0.4, 0.3, 0.3), nrow = 3, byrow = T)
d2 <- round(matrix.power(p2, 2000)[1,], 3)


# Data for Cluster 3

p3 <- matrix(c(0.025, 0.95, 0.025, 0.025, 0.025, 0.95, 0.4, 0.4, 0.2), nrow = 3, byrow = T)
e3 <- matrix(c(0.05, 0.90, 0.05, 0.05, 0.05, 0.90, 0.90, 0.05, 0.05), nrow = 3, byrow = T)
d3 <- round(matrix.power(p3, 2000)[1,], 3)

# Data for Cluster 4

p4 <- matrix(c(0.9, 0.05, 0.05, 0.1, 0.8, 0.1, 0.4, 0.4, 0.2), nrow = 3, byrow = T)
e4 <- matrix(rep(1/3,9), nrow = 3, byrow = T)
d4 <- round(matrix.power(p4, 2000)[1,], 3)

# Simulating 5 series from each cluster

vector_l <- c(200)

m <- 2 # Fuzziness coefficient

TP <- array(rep(NA, 3 * 3 * 4), c(3, 3, 4)) 
TP[,,1] <- p1
TP[,,2] <- p2
TP[,,3] <- p3
TP[,,4] <- p4

# set.seed(321)





ari_cramer <- numeric()
ari_indicators <- numeric()
ari_cramer_cohen <- numeric()
ari_total <- numeric()
ari_hmm <- numeric()
ari_additional <- numeric()
ari_cadez <- numeric()
ari_magariños_vilar <- numeric()
ari_mv <- numeric()
jaccard_cramer <- numeric()
jaccard_indicators <- numeric()
jaccard_cramer_cohen <- numeric()
jaccard_total <- numeric()
jaccard_hmm <- numeric()
jaccard_additional <- numeric()
jaccard_cadez <- numeric()
jaccard_magariños_vilar <- numeric()
jaccard_mv <- numeric()
larsen_cramer <- numeric()
larsen_indicators <- numeric()
larsen_cramer_cohen <- numeric()
larsen_total <- numeric()
larsen_hmm <- numeric()
larsen_additional <- numeric()
larsen_cadez <- numeric()
larsen_magariños_vilar <- numeric()
larsen_mv <- numeric()

time_mv <- list()


for (l in vector_l) {
  
  for (i in 1 : 500) {

  
  
  cluster1 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = e3, emission_probs = e3, sequence_length = l)
  cluster2 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = e3, emission_probs = p1, sequence_length = l)
  cluster3 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = e2, emission_probs = p1, sequence_length = l)
  cluster4 <- simulate_hmm(n_sequences = 5, initial_probs = d1, transition_probs = e4, emission_probs = e4, sequence_length = l)
  
  # Converting the series in factors
  

  cluster1 <- apply(cluster1$observations, 2, factor)
  cluster2 <- apply(cluster2$observations, 2, factor)
  cluster3 <- apply(cluster3$observations, 2, factor)
  cluster4 <- apply(cluster4$observations, 2, factor)
  cluster <- rbind(cluster1, cluster2, cluster3, cluster4)
  n_cluster1 <- apply(cluster1, 2, as.numeric)
  n_cluster2 <- apply(cluster2, 2, as.numeric)
  n_cluster3 <- apply(cluster3, 2, as.numeric)
  n_cluster4 <- apply(cluster4, 2, as.numeric)
  n_cluster <- rbind(n_cluster1, n_cluster2, n_cluster3, n_cluster4)
  ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5))
  
  cluster_df <- matrix(0, nrow = nrow(cluster), ncol = ncol(cluster))
  cluster_df <- as.data.frame(cluster_df)
  
  for (q in 1 : (K*n)) {
    
    cluster_df[q,] <- factor(cluster[q,])
    
  }


    
    # Algorithms
    
    # cramer
    
    categories <- factor(c('1', '2', '3'))
    clusterExport(c1, c('cluster', 'p_i_j_k_function'))
    
    features_cramer <- t(parApply(c1, cluster_df, 1, features_cramer_v_function, categories = categories,
                                   max_lag = 1))
    
    clustering_cramer <- pam(features_cramer, k = K)
    fuzzy_matrix_cramer <- clustering_cramer$clustering
    ari_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer)
    jaccard_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer,
                                              method = 'jaccard_index')
    larsen_cramer[[i]] <- external_validation(ground_truth, fuzzy_matrix_cramer,
                                               method = 'fowlkes_mallows_index')
    
    # HMM fit
    
  features_hmm <- matrix(0, nrow = 20, ncol = 18)
    
    for (v in 1 :20){
      
     features_hmm[v,] <- hmm_fit(cluster_df[v,])
      
    }
    
    
    
   clustering_hmm <-pam(features_hmm, k = K)
   fuzzy_matrix_hmm <- clustering_hmm$clustering
   ari_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm)
   jaccard_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm,
                                             method = 'jaccard_index')
   larsen_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm,
                                             method = 'fowlkes_mallows_index')
    
    # Indicators
    
    features_marginal <-  t(parApply(c1, cluster_df, 1, marginal_probabilities, categories = categories))
    features_indicators <- t(parApply(c1, cluster_df, 1, features_categorical_series, levels_dataset = categories,
                                      l = 3))
    
    clustering_indicators <- pam((cbind(features_indicators, features_marginal)), k = K)
    fuzzy_matrix_indicators <- clustering_indicators$clustering
    ari_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators)
    jaccard_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                                 method = 'jaccard_index')
    larsen_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                                  method = 'fowlkes_mallows_index')
    
    # Indicators + total
    
    features_cohen <-  t(parApply(c1, cluster_df, 1, features_cohen_function, categories = categories,
                                  max_lag = 2))
    clustering_additional <- pam(scale(cbind(features_cramer, features_cohen, features_marginal, features_indicators)), k = K)
    fuzzy_matrix_additional <- clustering_additional$clustering
    ari_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional)
    jaccard_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional,
                                                 method = 'jaccard_index')
    
    
    
    
    # cramer and Cohens
    
    
    
    
    clustering_total <- pam((cbind(features_cramer, features_cohen, features_marginal)), k = K)
    fuzzy_matrix_total <- clustering_total$clustering
    ari_total[i] <- external_validation(ground_truth, fuzzy_matrix_total)
    jaccard_total[i] <- external_validation(ground_truth, fuzzy_matrix_total, 
                                            method = 'jaccard_index')
    larsen_total[[i]] <- external_validation(ground_truth, fuzzy_matrix_total,
                                             method = 'fowlkes_mallows_index')
    
    
    # Cadez
    
    n_cluster_list <- list()
    
    for (r in 1 : 20) {
      
      n_cluster_list[[r]] <- n_cluster[r,]
      
    }
    
    data_cadez <- click.read(n_cluster_list)
    cadez <- click.EM(X = data_cadez$X, K = K, r = 10)
    ari_cadez[i] <- external_validation(ground_truth, cadez$id)
    jaccard_cadez[i] <- external_validation(ground_truth, cadez$id,
                                            method = 'jaccard_index')
    larsen_cadez[[i]] <- external_validation(ground_truth, cadez$id,
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
     #                                    method = 'jaccard_index')
    # larsen_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv,
    #                                    method = 'fowlkes_mallows_index')
    
    print(i)
    
    
    
  }
  
  if (l == 200) {
    
    save(ari_hmm, file = 'ari_hmm_200.RData')
    save(ari_total, file = 'ari_total_200.RData')
    save(ari_indicators, file = 'ari_indicators_200.RData')
    save(ari_cadez, file = 'ari_cadez_200.RData')
    save(ari_mv, file = 'ari_mv_200.RData')
    save(jaccard_hmm, file = 'jaccard_hmm_200.RData')
    save(jaccard_total, file = 'jaccard_total_200.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_200.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_200.RData')
    save(jaccard_mv, file = 'jaccard_mv_200.RData')
    save(larsen_hmm, file = 'larsen_hmm_200.RData')
    save(larsen_total, file = 'larsen_total_200.RData')
    save(larsen_indicators, file = 'larsen_indicators_200.RData')
    save(larsen_cadez, file = 'larsen_cadez_200.RData')
    save(larsen_mv, file = 'larsen_mv_200.RData')
    
  }
  
  
  if (l == 600) {
    
    save(ari_hmm, file = 'ari_hmm_600.RData')
    save(ari_total, file = 'ari_total_600.RData')
    save(ari_indicators, file = 'ari_indicators_600.RData')
    save(ari_cadez, file = 'ari_cadez_600.RData')
    save(ari_mv, file = 'ari_mv_600.RData')
    save(jaccard_hmm, file = 'jaccard_hmm_600.RData')
    save(jaccard_total, file = 'jaccard_total_600.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_600.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_600.RData')
    save(jaccard_mv, file = 'jaccard_mv_600.RData')
    save(larsen_hmm, file = 'larsen_hmm_600.RData')
    save(larsen_total, file = 'larsen_total_600.RData')
    save(larsen_indicators, file = 'larsen_indicators_600.RData')
    save(larsen_cadez, file = 'larsen_cadez_600.RData')
    save(larsen_mv, file = 'larsen_mv_600.RData')
    
  }
  
}



# Creating the table

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/crisp_algorithms/scenario_2')





table_previous <- matrix(NA, nrow = 5, ncol = 6)
load('ari_hmm_200.RData')
load('ari_total_200.RData')
load('ari_indicators_200.RData')
load('ari_cadez_200.RData')
load('ari_mv_200.RData')
load('jaccard_hmm_200.RData')
load('jaccard_total_200.RData')
load('jaccard_indicators_200.RData')
load('jaccard_cadez_200.RData')
loead('jaccard_mv_200.Rdata')
load('larsen_mv_200.RData')
load('larsen_hmm_200.RData')
load('larsen_total_200.RData')
load('larsen_indicators_200.RData')
load('larsen_cadez_200.RData')
load('larsen_mv_200.RData')
table_previous[1, 1] <- mean(ari_total)
table_previous[1, 2] <- mean(jaccard_total)
table_previous[1, 3] <- mean(larsen_total)
table_previous[2, 1] <- mean(ari_indicators)
table_previous[2, 2] <- mean(jaccard_indicators)
table_previous[2, 3] <- mean(larsen_indicators)
table_previous[3, 1] <- mean(ari_hmm)
table_previous[3, 2] <- mean(jaccard_hmm)
table_previous[3, 3] <- mean(larsen_hmm)
table_previous[4, 1] <- mean(ari_cadez)
table_previous[4, 2] <- mean(jaccard_cadez)
table_previous[4, 3] <- mean(larsen_cadez)
table_previous[5, 1] <- mean(ari_mv)
table_previous[5, 2] <- mean(jaccard_mv)
table_previous[5, 3] <- mean(larsen_mv)



load('ari_hmm_600.RData')
load('ari_total_600.RData')
load('ari_indicators_600.RData')
load('ari_cadez_600.RData')
load('ari_mv_600.RData')
load('jaccard_hmm_600.RData')
load('jaccard_total_600.RData')
load('jaccard_indicators_600.RData')
load('jaccard_cadez_600.RData')
load('jaccard_mv_600.RData')
load('larsen_hmm_600.RData')
load('larsen_total_600.RData')
load('larsen_indicators_600.RData')
load('larsen_cadez_600.RData')
load('larsen_mv_600.RData')
table_previous[1, 4] <- mean(ari_total)
table_previous[1, 5] <- mean(jaccard_total)
table_previous[1, 6] <- mean(larsen_total)
table_previous[2, 4] <- mean(ari_indicators)
table_previous[2, 5] <- mean(jaccard_indicators)
table_previous[2, 6] <- mean(larsen_indicators)
table_previous[3, 4] <- mean(ari_hmm)
table_previous[3, 5] <- mean(jaccard_hmm)
table_previous[3, 6] <- mean(larsen_hmm)
table_previous[4, 4] <- mean(ari_cadez)
table_previous[4, 5] <- mean(jaccard_cadez)
table_previous[4, 6] <- mean(larsen_cadez)
table_previous[5, 4] <- mean(ari_mv)
table_previous[5, 5] <- mean(jaccard_mv)
table_previous[5, 6] <- mean(larsen_mv)


xtable(table_previous)
