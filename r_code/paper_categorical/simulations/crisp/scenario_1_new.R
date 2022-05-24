
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



ari_total <- numeric()
jaccard_total <- numeric()
larsen_total <- numeric()

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

  # Cramer

  categories <- factor(c('a', 'b', 'c'))
  clusterExport(c1, c('cluster', 'p_i_j_k_function'))

  features_cramer <- listTomatrix(parLapply(c1, cluster, features_cramer_v_function, categories = categories,
                                             max_lag = 1))



  features_cohen <- listTomatrix(parLapply(c1, cluster, features_cohen_function, categories = categories,
                                           max_lag = 1))
  features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))

  clustering_total <- pam((cbind(features_cramer, features_cohen, features_marginal)), k = K)
  fuzzy_matrix_total <- clustering_total$clustering
  ari_total[i] <- external_validation(ground_truth, fuzzy_matrix_total)
  jaccard_total[i] <- external_validation(ground_truth, fuzzy_matrix_total,
                                          method = 'jaccard_index')
  larsen_total[i] <- external_validation(ground_truth, fuzzy_matrix_total,
                                            method = 'fowlkes_mallows_index')

  print(i)


  if (l == 200) {

    save(ari_total, file = 'ari_total_200.RData')
    save(jaccard_total, file = 'jaccard_total_200.RData')
    save(larsen_total, file = 'larsen_total_200.RData')


  }


  if (l == 600) {

    save(ari_total, file = 'ari_total_600.RData')
    save(jaccard_total, file = 'jaccard_total_600.RData')
    save(larsen_total, file = 'larsen_total_600.RData')

  }

}

}

