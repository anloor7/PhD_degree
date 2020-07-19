

# Loading the datasets

sp500 <- read.csv('all_stocks_5yr.csv')
head(sp500)
sector <- read.csv('sp5_sector.csv')
head(sector)


# Performing an inner join 

library(sqldf)
sp500_labels <- sqldf('SELECT close, volume, sp500.Name as symbol, sector.Name as name, sector.Sector as sector
                      from sp500 inner join sector on sector.Symbol = sp500.Name ')
head(sp500_labels)

# Ordering the datasets by symbols 

sp500 <- sqldf('SELECT * from sp500 order by Name')
sector <- sqldf('SELECT * from sector order by Symbol')
sp500_labels <- sqldf('SELECT * from sp500_labels order by symbol')
head(sp500)
head(sector)
head(sp500_labels)
unique_symbols <- sqldf('SELECT distinct(symbol) from sp500_labels')

# There are only 443 companies appearing in both databases 

# Creating the 443 bivariate MTS 

mts <- list()

for (i in unique(sp500_labels$symbol)) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Removing the the series whose length is less than 1258 and the corresponding labels 

g <- numeric()

for (i in 1:443) {
  
  g[i] <- nrow(mts[[i]])
  
}

indexnull <- which(g != 1258)
indexnotnull <- which(g == 1258)
length(indexnotnull)

# Creating the 414 bivariate MTS 

mts <- list()

for (i in unique(sp500_labels$symbol)[indexnotnull]) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Lets perform the clustering 

coherence2 <- listTomatrix(lapply(mts, quantile_coherence_re_im))
clustering <- kmeans(coherence2, 3)$cluster

# Ground truth 

unique_symbols_filtered <- sqldf(' SELECT unique_symbols.symbol, sector.Sector from unique_symbols 
                                 inner join sector on sector.Symbol = unique_symbols.symbol ')


ground_truth <- as.numeric(factor((unique_symbols_filtered$Sector)[indexnotnull]))
external_validation(ground_truth, clustering)


# ARI in 3 sectors 

index_3_sectors_df <- unique_symbols_filtered[indexnotnull,]
index_3_sectors <- index_3_sectors_df$Sector[index_3_sectors_df$Sector == 'Utilities' 
                                                  | index_3_sectors_df$Sector == 'Consumer Discretionary' 
                                 | index_3_sectors_df$Sector == 'Financials']
ground_truth_3_sectors <- as.numeric(factor((index_3_sectors)))

where_3_sectors <- index_3_sectors_df[index_3_sectors_df$Sector == 'Utilities' 
                                             | index_3_sectors_df$Sector == 'Consumer Discretionary' 
                                             | index_3_sectors_df$Sector == 'Financials',]
where_3_sectors_symbol <- where_3_sectors$symbol

# Creating the 136 bivariate MTS 

mts <- list()

for (i in unique(where_3_sectors_symbol)) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Removing the the series whose length is less than 1258 and the corresponding labels 

g <- numeric()

for (i in 1:136) {
  
  g[i] <- nrow(mts[[i]])
  
}

indexnull <- which(g != 1258)
indexnotnull <- which(g == 1258)
length(indexnotnull)

# Creating the 136 bivariate MTS 

mts <- list()

for (i in unique(where_3_sectors_symbol)[indexnotnull]) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Lets perform the clustering 

coherence2 <- listTomatrix(lapply(mts, quantile_coherence_re_im))
clustering <- kmeans(coherence2, 3)$cluster
external_validation(ground_truth_3_sectors, clustering, summary_stats = T)



# ARI in 2 sectors 

sp500 <- read.csv('all_stocks_5yr.csv')
head(sp500)
sector <- read.csv('sp5_sector.csv')
head(sector)


index_2_sectors_df <- index_3_sectors_df
index_2_sectors <- index_2_sectors_df$Sector[index_2_sectors_df$Sector == 'Utilities' 
                                             | index_2_sectors_df$Sector == 'Financials']
ground_truth_2_sectors <- as.numeric(factor((index_2_sectors)))

where_2_sectors <- index_2_sectors_df[index_2_sectors_df$Sector == 'Utilities' 
                                      | index_2_sectors_df$Sector == 'Financials',]
where_2_sectors_symbol <- where_2_sectors$symbol

# Creating the 80 bivariate MTS 

mts <- list()

for (i in unique(where_2_sectors_symbol)) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Removing the the series whose length is less than 1258 and the corresponding labels 

g <- numeric()

for (i in 1:80) {
  
  g[i] <- nrow(mts[[i]])
  
}

indexnull <- which(g != 1258)
indexnotnull <- which(g == 1258)
length(indexnotnull)

# Creating the 80 bivariate MTS 

mts <- list()

for (i in unique(where_2_sectors_symbol)[indexnotnull]) {
  
  indexes <- which(sp500_labels$symbol == i)
  returns <- scale(Delt(sp500_labels$close[indexes], type = 'log'))
  volumes <- scale(Delt(sp500_labels$volume[indexes], type = 'log'))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Lets perform the clustering 

coherence2 <- listTomatrix(lapply(mts, quantile_coherence_re_im))
clustering <- pam(coherence2, 2)$cluster
external_validation(ground_truth_2_sectors, clustering, summary_stats = T)

# Wavelets

features <- lapply(mts, wavelet_features, wf = wf, J = J) 
dis_matrix <- proxy::dist(features, wave_dist)  
clustering <- pam(dis_matrix, 2)$cluster
external_validation(ground_truth_2_sectors, clustering, summary_stats = T)

# Alonso y Peña 

features <- listTomatrix(lapply(mts, gcc_features_mts))
dis_matrix <- proxy::dist(features, EuclideanDistance)
clustering <- pam(dis_matrix, 2)$cluster
external_validation(ground_truth_2_sectors, clustering,summary_stats = T)

# KST

dis_matrix <- matrix(0, 80, 80)
for (i in 1 : 80) {
  for (j in 1 : 80) {
    dis_matrix[i, j] <- j_divergence(mts[[i]], mts[[j]])
  }
}

diag(dis_matrix) <- 0 # Numerical error 
clustering <- pam(dis_matrix, 2)$cluster
external_validation(ground_truth_2_sectors, clustering)

# Graphical representations of 9 first series in each cluster

financials_9 <- where_2_sectors$symbol[where_2_sectors$Sector == 'Financials'][1:9]
utilities_9 <- where_2_sectors$symbol[where_2_sectors$Sector == 'Utilities'][1:9]

indexes_financials_9 <- which(where_2_sectors$symbol %in% financials_9)
indexes_utilities_9 <- which(where_2_sectors$symbol %in% utilities_9)

mts_financials_9 <- mts[indexes_financials_9]
mts_utilities_9 <- mts[indexes_utilities_9]

plot_financials_9 <- list()
plot_utilities_9 <- list()

for (i in 1:9) {
  
  plot_financials_9[[i]] <- mts_plot(mts_financials_9[[i]], title = financials_9[i])
  plot_utilities_9[[i]] <- mts_plot(mts_utilities_9[[i]], title = utilities_9[i])
  
}

do.call("grid.arrange", c(plot_financials_9[1 : 9], ncol = 3))
do.call("grid.arrange", c(plot_utilities_9[1 : 9], ncol = 3))

