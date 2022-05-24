
setwd('/Users/angel/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/crisp_algorithms/scenario_1')

library(SMM)
library(markovchain)
library(cluster)
library(ClusterR)
library(ClickClust)

# Simulating an Scenario of Markov processes of order 1

E <- c('a', 'b', 'c') # State space
s_e <- length(E)
K <- 4


# Data for Cluster 1

p1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
d1 <- matrix.power(p1, 5000)[1,]


# Data for Cluster 2

p2 <- matrix(c(0.1, 0.8, 0.1, 0.6, 0.3, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
d2 <- matrix.power(p2, 1000)[1,]


# Data for Cluster 3

p3 <- matrix(c(0.05, 0.90, 0.05, 0.05, 0.05, 0.90, 0.90, 0.05, 0.05), nrow = 3, byrow = T)
d3 <- matrix.power(p3, 1000)[1,]

# Data for Cluster 4

p4 <- matrix(rep(1/3, 9), nrow = 3, byrow = T)
d4 <- matrix.power(p4, 1000)[1,]


vector_l <- c(200, 600)
m <- 2 # Fuzziness coefficient

TP <- array(rep(NA, 3 * 3 * 4), c(3, 3, 4)) 
TP[,,1] <- p1
TP[,,2] <- p2
TP[,,3] <- p3
TP[,,4] <- p4

# Simulating 5 series from each cluster



set.seed(1234)



ari_cramer <- numeric()
ari_indicators <- numeric()
ari_cramer_cohen <- numeric()
ari_total <- numeric()
ari_markov <- numeric()
ari_additional <- numeric()
ari_cadez <- numeric()
ari_magariños_vilar <- numeric()
ari_mv <- numeric()
jaccard_cramer <- numeric()
jaccard_indicators <- numeric()
jaccard_cramer_cohen <- numeric()
jaccard_total <- numeric()
jaccard_markov <- numeric()
jaccard_additional <- numeric()
jaccard_cadez <- numeric()
jaccard_magariños_vilar <- numeric()
jaccard_mv <- numeric()
larsen_cramer <- numeric()
larsen_indicators <- numeric()
larsen_cramer_cohen <- numeric()
larsen_total <- numeric()
larsen_markov <- numeric()
larsen_additional <- numeric()
larsen_cadez <- numeric()
larsen_magariños_vilar <- numeric()
larsen_mv <- numeric()

for (l in vector_l) {

for (i in 1 : 500) {
  
  cluster1 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p1, init = c(0.365, 0.523, 1-0.365-0.523), k = 1)
  cluster2 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p2, init = c(0.25, 0.25, 0.5), k = 1)
  cluster3 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p3, init = c(1/3, 1/3, 1/3), k = 1)
  cluster4 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p4, init = c(1/3, 1/3, 1/3), k = 1)
  
  
  # Converting the series in factors
  
  cluster1 <- lapply(cluster1, factor)
  cluster2 <- lapply(cluster2, factor)
  cluster3 <- lapply(cluster3, factor)
  cluster4 <- lapply(cluster4, factor)
  n_cluster1 <- lapply(cluster1, as.numeric)
  n_cluster2 <- lapply(cluster2, as.numeric)
  n_cluster3 <- lapply(cluster3, as.numeric)
  n_cluster4 <- lapply(cluster4, as.numeric)
  cluster <- c(cluster1, cluster2, cluster3, cluster4)
  n_cluster <- c(n_cluster1, n_cluster2, n_cluster3, n_cluster4)
  ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5))
  
  # Algorithms
  
  # cramer
  
  categories <- factor(c('a', 'b', 'c'))
  clusterExport(c1, c('cluster', 'p_i_j_k_function'))
  
  features_cramer <- listTomatrix(parLapply(c1, cluster, features_cramer_v_function, categories = categories,
                                             max_lag = 1))
  
  clustering_cramer <- pam(features_cramer, k = K)
  fuzzy_matrix_cramer <- clustering_cramer$clustering
  ari_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer)
  jaccard_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer,
                                            method = 'jaccard_index')
  larsen_cramer[i] <- external_validation(ground_truth, fuzzy_matrix_cramer,
                                             method = 'fowlkes_mallows_index')
  
  # Markov fit
  
  features_markov <- listTomatrix(parLapply(c1, cluster, markov_vector))
  clustering_markov <- pam(features_markov, k = K)
  fuzzy_matrix_markov <- clustering_markov$clustering
  ari_markov[i] <- external_validation(ground_truth, fuzzy_matrix_markov)
  jaccard_markov[i] <- external_validation(ground_truth, fuzzy_matrix_markov,
                                           method = 'jaccard_index')
  larsen_markov[i] <- external_validation(ground_truth, fuzzy_matrix_markov,
                                            method = 'fowlkes_mallows_index')
  
  # Indicators
  
  features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
  features_indicators <- listTomatrix(parLapply(c1, cluster, features_categorical_series, 
                                                levels_dataset = categories))
  
  clustering_indicators <- pam((cbind(features_indicators, features_marginal)), k = K)
  fuzzy_matrix_indicators <- clustering_indicators$clustering
  ari_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators)
  jaccard_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                               method = 'jaccard_index')
  larsen_indicators[i] <- external_validation(ground_truth, fuzzy_matrix_indicators,
                                                method = 'fowlkes_mallows_index')
  
  # Indicators + total
  
  features_cohen <- listTomatrix(parLapply(c1, cluster, features_cohen_function, categories = categories,
                                           max_lag = 1))
  features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
  clustering_additional <- pam(scale(cbind(features_cramer, features_cohen, features_marginal, features_indicators)), k = K)
  fuzzy_matrix_additional <- clustering_additional$clustering
  ari_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional)
  jaccard_additional[i] <- external_validation(ground_truth, fuzzy_matrix_additional,
                                               method = 'jaccard_index')
  
  
  
  
  # cramer and Cohens
  
  categories <- factor(c('a', 'b', 'c'))
  clusterExport(c1, c('cluster', 'p_i_j_k_function'))
  
  clustering_cramer_cohen <- pam(scale(cbind(features_cramer, features_cohen)), k = K)
  fuzzy_matrix_cramer_cohen <- clustering_cramer_cohen$clustering
  ari_cramer_cohen[i] <- external_validation(ground_truth, fuzzy_matrix_cramer_cohen)
  
  clustering_total <- pam((cbind(features_cramer, features_cohen, features_marginal)), k = K)
  fuzzy_matrix_total <- clustering_total$clustering
  ari_total[i] <- external_validation(ground_truth, fuzzy_matrix_total)
  jaccard_total[i] <- external_validation(ground_truth, fuzzy_matrix_total, 
                                          method = 'jaccard_index')
  larsen_total[i] <- external_validation(ground_truth, fuzzy_matrix_total,
                                            method = 'fowlkes_mallows_index')
  
  
  # Cadez
  
  data_cadez <- click.read(n_cluster)
  cadez <- click.EM(X = data_cadez$X, K = K, r = 10)
  ari_cadez[i] <- external_validation(ground_truth, cadez$id)
  jaccard_cadez[i] <- external_validation(ground_truth, cadez$id,
                                          method = 'jaccard_index')
  larsen_cadez[i] <- external_validation(ground_truth, cadez$id,
                                          method = 'fowlkes_mallows_index')
  
  # Magariños-Vilar
  
  if (i <= 20) {
  
  cluster_matrix <- matrix(NA, nrow = 20, ncol = l)
  
  for (p in 1 : 20) {
  
  cluster_matrix[p,] <- as.numeric(cluster[[p]])
  
  }
  
  magariños_vilar <- Diss.Mat(cluster_matrix, opt = "AD-Corr", kk = 1)
  pam_magariños_vilar <- pam(magariños_vilar$Dist, k = K)
  fuzzy_matrix_mv <- pam_magariños_vilar$clustering
  ari_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv)
  jaccard_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv,
                                          method = 'jaccard_index')
  larsen_mv[i] <- external_validation(ground_truth, fuzzy_matrix_mv,
                                         method = 'fowlkes_mallows_index')
  
  }
  
  print(i)

  
}
  
  if (l == 200) {
    
    save(ari_markov, file = 'ari_markov_200.RData')
    save(ari_total, file = 'ari_total_200.RData')
    save(ari_indicators, file = 'ari_indicators_200.RData')
    save(ari_cadez, file = 'ari_cadez_200.RData')
    save(ari_mv, file = 'ari_mv_200.RData')
    save(jaccard_markov, file = 'jaccard_markov_200.RData')
    save(jaccard_total, file = 'jaccard_total_200.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_200.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_200.RData')
    save(jaccard_mv, file = 'jaccard_mv_200.RData')
    save(larsen_markov, file = 'larsen_markov_200.RData')
    save(larsen_total, file = 'larsen_total_200.RData')
    save(larsen_indicators, file = 'larsen_indicators_200.RData')
    save(larsen_cadez, file = 'larsen_cadez_200.RData')
    save(larsen_mv, file = 'larsen_mv_200.RData')
    
  }
  
  
  if (l == 600) {
    
    save(ari_markov, file = 'ari_markov_600.RData')
    save(ari_total, file = 'ari_total_600.RData')
    save(ari_indicators, file = 'ari_indicators_600.RData')
    save(ari_cadez, file = 'ari_cadez_600.RData')
    save(ari_mv, file = 'ari_mv_600.RData')
    save(jaccard_markov, file = 'jaccard_markov_600.RData')
    save(jaccard_total, file = 'jaccard_total_600.RData')
    save(jaccard_indicators, file = 'jaccard_indicators_600.RData')
    save(jaccard_cadez, file = 'jaccard_cadez_600.RData')
    save(jaccard_mv, file = 'jaccard_mv_600.RData')
    save(larsen_markov, file = 'larsen_markov_600.RData')
    save(larsen_total, file = 'larsen_total_600.RData')
    save(larsen_indicators, file = 'larsen_indicators_600.RData')
    save(larsen_cadez, file = 'larsen_cadez_600.RData')
    save(larsen_mv, file = 'larsen_mv_600.RData')
    
  }
  
}



# Creating the table

setwd('/Users/angel/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/crisp_algorithms/scenario_1')





table_previous <- matrix(NA, nrow = 5, ncol = 6)
load('ari_markov_200.RData')
load('ari_total_200.RData')
load('ari_indicators_200.RData')
load('ari_cadez_200.RData')
load('ari_mv_200.RData')
load('jaccard_markov_200.RData')
load('jaccard_total_200.RData')
load('jaccard_indicators_200.RData')
load('jaccard_cadez_200.RData')
load('jaccard_mv_200.RData')
load('larsen_mv_200.RData')
load('larsen_markov_200.RData')
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
table_previous[3, 1] <- mean(ari_markov)
table_previous[3, 2] <- mean(jaccard_markov)
table_previous[3, 3] <- mean(larsen_markov)
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
load('jaccard_markov_600.RData')
load('jaccard_total_600.RData')
load('jaccard_indicators_600.RData')
load('jaccard_cadez_600.RData')
load('jaccard_mv_600.RData')
load('larsen_markov_600.RData')
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
table_previous[3, 4] <- mean(ari_markov)
table_previous[3, 5] <- mean(jaccard_markov)
table_previous[3, 6] <- mean(larsen_markov)
table_previous[4, 4] <- mean(ari_cadez)
table_previous[4, 5] <- mean(jaccard_cadez)
table_previous[4, 6] <- mean(larsen_cadez)
table_previous[5, 4] <- mean(ari_mv)
table_previous[5, 5] <- mean(jaccard_mv)
table_previous[5, 6] <- mean(larsen_mv)


xtable(table_previous)
