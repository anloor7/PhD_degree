

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_4')

# Simulating an Scenario of Markov processes of order 1

E <- c('a', 'b', 'c') # State space
s_e <- length(E)
K <- 2


# Data for Cluster 1

p1 <- matrix(c(0.1, 0.8, 0.1, 0.5, 0.4, 0.1, 0.6, 0.2, 0.2), nrow = 3, byrow = T)
d1 <- matrix.power(p1, 5000)[1,]


# Data for Cluster 2

p2 <- matrix(c(0.2, 0.7, 0.1, 0.4, 0.3, 0.3, 0.2, 0.4, 0.4), nrow = 3, byrow = T)
d2 <- matrix.power(p2, 1000)[1,]


# Data for Cluster 3

p3 <- (p1 + p2)/2
d3 <- matrix.power(p3, 1000)[1,]


cutoff <- 0.7
l <- c(200)
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

for (w in 1 : w_last) {
    
    cluster1 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p1, init = c(0.365, 0.523, 1-0.365-0.523), k = 1)
    cluster2 <- simulMk(E = E, nbSeq = 5, lengthSeq = rep(l, 5), Ptrans = p2, init = c(0.25, 0.25, 0.5), k = 1)
    cluster3 <- simulMk(E = E, nbSeq = 1, lengthSeq = rep(l, 1), Ptrans = p3, init = c(0.332, 0.406, 0.262), k = 1)
    
    
    # Converting the series in factors
    
    cluster1 <- lapply(cluster1, factor)
    cluster2 <- lapply(cluster2, factor)
    cluster3 <- lapply(cluster3, factor)
    n_cluster1 <- lapply(cluster1, as.numeric)
    n_cluster2 <- lapply(cluster2, as.numeric)
    n_cluster3 <- lapply(cluster3, as.numeric)
    cluster[[w]] <- c(cluster1, cluster2, cluster3)
    n_cluster[[w]] <- c(n_cluster1, n_cluster2, n_cluster3)
    ground_truth <- c(rep(1, 5), rep(2, 5), rep(3, 5))

 
}

    # Algorithms
    
    # Total
    
    categories <- factor(c('a', 'b', 'c'))
    clusterExport(c1, c('cluster', 'p_i_j_k_function'))
    
    set.seed(321)
    
    k = 1
    for (m in vector_m) {
      
      
      
      for (w in (1 : w_last)) {
        
        features_cramer <- listTomatrix(parLapply(c1, cluster[[w]], features_cramer_v_function, categories = categories,
                                                   max_lag = 1))
        features_marginal <- listTomatrix(parLapply(c1, cluster[[w]], marginal_probabilities, categories = categories))
        features_cohen <- listTomatrix(parLapply(c1, cluster[[w]], features_cohen_function, categories = categories,
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
    
    
    
    # Indicators
    
    set.seed(321)
    
    k = 1
    for (m in vector_m) {
      
      
      
      for (w in (1 : w_last)) {
    
    features_marginal <- listTomatrix(parLapply(c1, cluster[[w]], marginal_probabilities, categories = categories))
    features_indicators <- listTomatrix(parLapply(c1, cluster[[w]], features_categorical_series, 
                                                  levels_dataset = categories))
    
    clustering_indicators <- FKM(X = features_indicators, k = K, m = m)
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
      
      
      
      for (w in (1 : w_last)) {
        
        features_mle <- listTomatrix(parLapply(c1, cluster[[w]], markov_vector))
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

setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/simulation_study/fuzzy_algorithms/scenario_4')
load('table_l_200.RData')
load('table_l_600.RData')

xtable(rbind(table_l_200, table_l_600))    
      