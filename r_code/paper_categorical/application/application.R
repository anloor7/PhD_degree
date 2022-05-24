

setwd('/Users/angel/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2021/paper_categorical/application')

# Loading the data

rodent_1 <- factor(as.vector(read.fasta('rodent_1.fasta')[[1]]))
rodent_2 <- factor(as.vector(read.fasta('rodent_2.fasta')[[1]]))
rodent_3 <- factor(as.vector(read.fasta('rodent_3.fasta')[[1]]))
rodent_4 <- factor(as.vector(read.fasta('rodent_4.fasta')[[1]]))
rodent_5 <- factor(as.vector(read.fasta('rodent_5.fasta')[[1]]))
# rodent_6 <- factor(as.vector(read.fasta('rodent_6.fasta')[[1]]))
# rodent_7 <- factor(as.vector(read.fasta('rodent_7.fasta')[[1]]))
# rodent_8 <- factor(as.vector(read.fasta('rodent_8.fasta')[[1]]))


circo_10 <- factor(as.vector(read.fasta('circo_10.fasta')[[1]]))
circo_10[circo_10 == 'n'] <-  'a'
circo_10[circo_10 == 'y'] <-  'a'
circo_10 <- factor(circo_10)
circo_11 <- factor(as.vector(read.fasta('circo_11.fasta')[[1]]))
circo_11[circo_11 == 'k'] <-  'a'
circo_11[circo_11 == 'm'] <-  'a'
circo_11[circo_11 == 'n'] <-  'a'
circo_11[circo_11 == 'r'] <-  'a'
circo_11[circo_11 == 'w'] <-  'a'
circo_11 <- factor(circo_11)
circo_13 <- factor(as.vector(read.fasta('circo_13.fasta')[[1]]))
circo_14 <- factor(as.vector(read.fasta('circo_14.fasta')[[1]]))
circo_14[circo_14 == 'n'] <-  'a'
circo_14 <- factor(circo_14)
circo_15 <- factor(as.vector(read.fasta('circo_15.fasta')[[1]]))
circo_15[circo_15 == 'k'] <-  'a'
circo_15[circo_15 == 'm'] <-  'a'
circo_15[circo_15 == 'n'] <-  'a'
circo_15[circo_15 == 'r'] <-  'a'
circo_15[circo_15 == 's'] <-  'a'
circo_15[circo_15 == 'w'] <-  'a'
circo_15[circo_15 == 'y'] <-  'a'
circo_15 <- factor(circo_15)
circo_16 <- factor(as.vector(read.fasta('circo_16.fasta')[[1]]))
circo_16[circo_16 == 'k'] <-  'a'
circo_16[circo_16 == 'n'] <-  'a'
circo_16[circo_16 == 's'] <-  'a'
circo_16 <- factor(circo_16)
circo_18 <- factor(as.vector(read.fasta('circo_18.fasta')[[1]]))
circo_18[circo_18 == 'm'] <-  'a'
circo_18[circo_18 == 'n'] <-  'a'
circo_18 <- factor(circo_18)
circo_19 <- factor(as.vector(read.fasta('circo_19.fasta')[[1]]))
circo_19[circo_19 == 'n'] <-  'a'
circo_19 <- factor(circo_19)
circo_21 <- factor(as.vector(read.fasta('circo_21.fasta')[[1]]))
circo_21[circo_21 == 'n'] <-  'a'
circo_21[circo_21 == 'r'] <-  'a'
circo_21[circo_21 == 'w'] <-  'a'
circo_21[circo_21 == 'y'] <-  'a'
circo_21 <- factor(circo_21)
circo_2 <- factor(as.vector(read.fasta('circo_2.fasta')[[1]]))
circo_2[circo_2 == 'n'] <-  'a'
circo_2[circo_2 == 's'] <-  'a'
circo_2 <- factor(circo_2)


cosavirus_10 <- factor(as.vector(read.fasta('cosavirus_10.fasta')[[1]]))
cosavirus_11 <- factor(as.vector(read.fasta('cosavirus_11.fasta')[[1]]))
cosavirus_12 <- factor(as.vector(read.fasta('cosavirus_12.fasta')[[1]]))
cosavirus_13 <- factor(as.vector(read.fasta('cosavirus_13.fasta')[[1]]))
cosavirus_14 <- factor(as.vector(read.fasta('cosavirus_14.fasta')[[1]]))
cosavirus_15 <- factor(as.vector(read.fasta('cosavirus_15.fasta')[[1]]))
cosavirus_15[cosavirus_15 == 'r'] <-  'a'
cosavirus_15 <- factor(cosavirus_15)
cosavirus_16 <- factor(as.vector(read.fasta('cosavirus_16.fasta')[[1]]))
cosavirus_17 <- factor(as.vector(read.fasta('cosavirus_17.fasta')[[1]]))
cosavirus_18 <- factor(as.vector(read.fasta('cosavirus_18.fasta')[[1]]))
cosavirus_19 <- factor(as.vector(read.fasta('cosavirus_19.fasta')[[1]]))


parecho_11 <- factor(as.vector(read.fasta('parecho_11.fasta')[[1]]))
parecho_13 <- factor(as.vector(read.fasta('parecho_13.fasta')[[1]]))
parecho_16 <- factor(as.vector(read.fasta('parecho_16.fasta')[[1]]))
parecho_10 <- factor(as.vector(read.fasta('parecho_10.fasta')[[1]]))
parecho_12 <- factor(as.vector(read.fasta('parecho_12.fasta')[[1]]))
parecho_9 <- factor(as.vector(read.fasta('parecho_9.fasta')[[1]]))
parecho_15 <- factor(as.vector(read.fasta('parecho_15.fasta')[[1]]))


# Applying the clustering procedures

cluster <- list()
cluster[[1]] <- rodent_1
cluster[[2]] <- rodent_2
cluster[[3]] <- rodent_3
cluster[[4]] <- rodent_4
cluster[[5]] <- rodent_5
cluster[[6]] <- circo_10
cluster[[7]] <- circo_11
cluster[[8]] <- circo_13
cluster[[9]] <- circo_14
cluster[[10]] <- circo_15
cluster[[11]] <- circo_16
cluster[[12]] <- circo_18
cluster[[13]] <- circo_19
cluster[[14]] <- circo_21
cluster[[15]] <- circo_2
cluster[[16]] <- cosavirus_10
cluster[[17]] <- cosavirus_11
cluster[[18]] <- cosavirus_12
cluster[[19]] <- cosavirus_13
cluster[[20]] <- cosavirus_14
cluster[[21]] <- cosavirus_15
cluster[[22]] <- cosavirus_16
cluster[[23]] <- cosavirus_17
cluster[[24]] <- cosavirus_18
cluster[[25]] <- cosavirus_19
cluster[[26]] <- parecho_11
cluster[[27]] <- parecho_13
cluster[[28]] <- parecho_16
cluster[[29]] <- parecho_10
cluster[[30]] <- parecho_12
cluster[[31]] <- parecho_9
cluster[[32]] <- parecho_15


n_cluster <- list()
n_cluster[[1]] <- as.numeric(rodent_1)
n_cluster[[2]] <- as.numeric(rodent_2)
n_cluster[[3]] <- as.numeric(rodent_3)
n_cluster[[4]] <- as.numeric(rodent_4)
n_cluster[[5]] <- as.numeric(rodent_5)
n_cluster[[6]] <- as.numeric(circo_10)
n_cluster[[7]] <- as.numeric(circo_11)
n_cluster[[8]] <- as.numeric(circo_13)
n_cluster[[9]] <- as.numeric(circo_14)
n_cluster[[10]] <- as.numeric(circo_15)
n_cluster[[11]] <- as.numeric(circo_16)
n_cluster[[12]] <- as.numeric(circo_18)
n_cluster[[13]] <- as.numeric(circo_19)
n_cluster[[14]] <- as.numeric(circo_21)
n_cluster[[15]] <- as.numeric(circo_2)
n_cluster[[16]] <- as.numeric(cosavirus_10)
n_cluster[[17]] <- as.numeric(cosavirus_11)
n_cluster[[18]] <- as.numeric(cosavirus_12)
n_cluster[[19]] <- as.numeric(cosavirus_13)
n_cluster[[20]] <- as.numeric(cosavirus_14)
n_cluster[[21]] <- as.numeric(cosavirus_15)
n_cluster[[22]] <- as.numeric(cosavirus_16)
n_cluster[[23]] <- as.numeric(cosavirus_17)
n_cluster[[24]] <- as.numeric(cosavirus_18)
n_cluster[[25]] <- as.numeric(cosavirus_19)
n_cluster[[26]] <- as.numeric(parecho_11)
n_cluster[[27]] <- as.numeric(parecho_13)
n_cluster[[28]] <- as.numeric(parecho_16)
n_cluster[[29]] <- as.numeric(parecho_10)
n_cluster[[30]] <- as.numeric(parecho_12)
n_cluster[[31]] <- as.numeric(parecho_9)
n_cluster[[32]] <- as.numeric(parecho_15)


ground_truth <- c(rep(1, 5), rep(2, 10), rep(3, 10), rep(4, 7))


# Total

categories <- factor(c('a', 'g', 'c', 't'))
clusterExport(c1, c('cluster', 'p_i_j_k_function'))

max_lag <- 1
features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
features_cramer <- listTomatrix(parLapply(c1, cluster, features_cramer_v_function, categories = categories,
                                           max_lag = max_lag))
features_cohen <- listTomatrix(parLapply(c1, cluster, features_cohen_function, categories = categories,
                                         max_lag = max_lag))

clustering_total <- pam((cbind(features_cramer, features_cohen, features_marginal)), k = 4)
fuzzy_matrix_total <- clustering_total$clustering
external_validation(ground_truth, fuzzy_matrix_total)
external_validation(ground_truth, fuzzy_matrix_total, method = 'jaccard_index')
external_validation(ground_truth, fuzzy_matrix_total, method = 'fowlkes_mallows_index')
intCriteria(cbind(features_cramer, features_cohen, features_marginal), part = fuzzy_matrix_total, crit = 'Xie_Beni')$xie_beni


# Indicators

categories <- factor(c('a', 'g', 'c', 't'))

max_lag <- 6
clusterExport(c1, c('cluster', 'features_categorical_series'))
features_marginal <- listTomatrix(parLapply(c1, cluster, marginal_probabilities, categories = categories))
features_indicators <- listTomatrix(lapply(cluster, 
                                              function(x) {features_categorical_series(x, 
                                                                                             levels_dataset = categories, l = max_lag)$whole_features}))

clustering_indicators <- pam((cbind(features_indicators, features_marginal)), k = 4)
fuzzy_matrix_indicators <- clustering_indicators$clustering
external_validation(ground_truth, fuzzy_matrix_indicators)
external_validation(ground_truth, fuzzy_matrix_indicators, method = 'jaccard_index')
external_validation(ground_truth, fuzzy_matrix_indicators, method = 'fowlkes_mallows_index')
intCriteria(cbind(features_indicators, features_marginal), part = fuzzy_matrix_indicators, crit = 'Xie_Beni')$xie_beni


# Cadez

data_cadez <- click.read(n_cluster)
cadez <- click.EM(X = data_cadez$X, K = 4, r = 2000, y = ground_truth)
ari_cadez[i] <- external_validation(ground_truth, cadez$id)
jaccard_cadez[i] <- external_validation(ground_truth, cadez$id,
                                        method = 'jaccard_index')
larsen_cadez[i] <- external_validation(ground_truth, cadez$id,
                                       method = 'fowlkes_mallows_index')

# MLE-Markov chain

features_markov <- listTomatrix(parLapply(c1, cluster, markov_vector))
clustering_markov <- pam((features_markov), k = 4)
fuzzy_matrix_markov <- clustering_markov$clustering
external_validation(ground_truth, fuzzy_matrix_markov)
jaccard_markov[i] <- external_validation(ground_truth, fuzzy_matrix_markov,
                                         method = 'jaccard_index')
larsen_markov[i] <- external_validation(ground_truth, fuzzy_matrix_markov,
                                        method = 'fowlkes_mallows_index')


# HMM

# Converting the series in factors


cluster1 <- apply(cluster$observations, 2, factor)
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

# HMM fit

cluster_df <- list()

for (i in 1 : 32) {
  
  cluster_df[[i]] <- t(as.data.frame(cluster[[i]]))
  
}



features_hmm <- matrix(0, nrow = 32, ncol = 21)

for (v in 1 :32){
  
  features_hmm[v,] <- hmm_fit(cluster_df[[v]])
  print(v)
  
}



clustering_hmm <- pam(features_hmm, k = 4)
fuzzy_matrix_hmm <- clustering_hmm$clustering
ari_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm)
jaccard_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm,
                                      method = 'jaccard_index')
larsen_hmm[i] <- external_validation(ground_truth, fuzzy_matrix_hmm,
                                     method = 'fowlkes_mallows_index')

# MLE-DAR(1)

list_matrices <- list()

for (h in 1 : 32) {
  
  list_matrices[[h]] <- bincodes[cluster[[h]],][, (1 : 4)]
  
}

matrix_features <- matrix(0, nrow = 32, ncol = 5)

for (h in 1 : 32) {
  
  auxiliary_computation <- suppressWarnings(constrOptim(c(1/5, 1/5, 1/5, 1/5, 1/5), ll_darp, NULL, ui=Amat(4,2), ci=bvec(4,2), databin = list_matrices[[h]]))
  matrix_features[h,] <- c(auxiliary_computation$par)
  print(h)
  
}

clustering_weib <- pam(matrix_features[,c(1 : 5)], k = 4)
fuzzy_matrix_weib <- clustering_weib$clustering
ari_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib)
jaccard_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib, 
                                       method = 'jaccard_index')
larsen_weib[i] <- external_validation(ground_truth, fuzzy_matrix_weib,
                                      method = 'fowlkes_mallows_index')



# 2 dimensional scaling plots

#d_{AM}

# Multidimensional scaling

dis_matrix <- proxy::dist((cbind(features_cramer, features_cohen, features_marginal)))^2
mds <- cmdscale(dis_matrix, 2, list. = F)
# mds$GOF
factor_color <- c(rep(1, 5), rep(2, 10), rep(3, 10), rep(4, 7))
df12 <- data.frame(cbind(mds), factor(factor_color))
plot(df12$X1, df12$X2, col = df12$factor.factor_color)
legend(-1.5 ,0.2 , unique(df12$factor.factor_color), col = 1 : 4, pch = 1)
colnames(df12)[3] <- 'series'



# Distance matrix for 2d reduced space

dis_matrix_reduced <- proxy::dist(mds)^2

sqrt(sum((dis_matrix_reduced-dis_matrix)^2)/sum((dis_matrix)^2))


plot_am <- ggplot(df12, aes(x = X1, y = X2, col = series )) + geom_point(size = 2) + 
  xlab('Coordinate 1') + ylab('Coordinate 2') + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,  size = 18),
        legend.position = '') + ggtitle(TeX('$d_{CC}$')) + xlim(c(-0.017, 0.017))


plot_auxiliary <- ggplot(df12, aes(x = X1, y = X2, col = series)) + geom_point(size = 2) + 
  xlab('Coordinate 1') + ylab('Coordinate 2') + 
  scale_color_discrete(labels = c("RAC", "CLDMD", "HC", "HP")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,  size = 18),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) + ggtitle(TeX('$d_{CC}$'))



extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(plot_auxiliary)


#d_{IP}

# Multidimensional scaling

dis_matrix <- proxy::dist((cbind(features_indicators, features_marginal)))^2
mds <- cmdscale(dis_matrix, 2, list. = F)
# mds$GOF
factor_color <- c(rep(1, 5), rep(2, 10), rep(3, 10), rep(4, 7))
df12 <- data.frame(cbind(mds), factor(factor_color))
plot(df12$X1, df12$X2, col = df12$factor.factor_color)
legend(-1.5 ,0.2 , unique(df12$factor.factor_color), col = 1 : 4, pch = 1)
colnames(df12)[3] <- 'series'



# Distance matrix for 2d reduced space

dis_matrix_reduced <- proxy::dist(mds)^2

sqrt(sum((dis_matrix_reduced-dis_matrix)^2)/sum((dis_matrix)^2))


plot_ip <- ggplot(df12, aes(x = X1, y = X2, col = series )) + geom_point(size = 2) + 
  xlab('Coordinate 1') + ylab('Coordinate 2') + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,  size = 18),
        legend.position = '') + ggtitle(TeX('$d_{B}$'))


plot_auxiliary <- ggplot(df12, aes(x = X1, y = X2, col = series)) + geom_point(size = 2) + 
  xlab('Coordinate 1') + ylab('Coordinate 2') + 
  scale_color_discrete(labels = c("RAC", "CLDMD", "HC", "HP")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,  size = 18),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) + ggtitle(TeX('$d_{B}$'))



extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(plot_auxiliary)


plot_am_ip <- grid.arrange(plot_am, plot_ip, ncol = 2)

plot_final <- grid.arrange(
  arrangeGrob(plot_am_ip), nrow = 2, ncol = 1, shared_legend,
  heights = c(40, 8))


# Fuzzy clustering methods

# Selection of m

# d_{AM}

set.seed(123456)
vector_m <- seq(1.1, 4, by = 0.1)
list_m_c <- list()
k <- 1
xie_beni_d_am <- numeric()

for (i in vector_m) {
    
    clustering_d_am <- FKM.med(X = cbind(features_cramer, features_cohen, features_marginal), k = 4, m = i, RS = 50)
    xie_beni_d_am[k] <- XB(cbind(features_cramer, features_cohen, features_marginal), 
                           clustering_d_am$U, clustering_d_am$H, m = 2)
    
    k <- k + 1
    print(i)
    
  }


# d_{IP}

set.seed(123456)
vector_m <- seq(1.1, 1.8, by = 0.1)
list_m_c <- list()
k <- 1
xie_beni_d_ip <- numeric()

for (i in vector_m) {
  
  clustering_d_ip <- FKM.med(X = cbind(features_indicators, features_marginal), k = 4, m = i, RS = 50, maxit = 100)
  xie_beni_d_ip[k] <- XB(cbind(features_indicators, features_marginal), 
                         clustering_d_ip$U, clustering_d_ip$H, m = 2)
  
  k <- k + 1
  print(i)
  
}


xtable(clustering_d_ip$U, digits = 3)
