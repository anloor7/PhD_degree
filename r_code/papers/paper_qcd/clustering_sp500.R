

# Lets cluster the 500 companies of S&P 500

sp500 <- read.csv('all_stocks_5yr.csv')
head(sp500)

# Creating the 505 bivariate MTS 

mts <- list()

for (i in unique(sp500$Name)) {
  
  indexes <- which(sp500$Name == i)
  returns <- scale(Delt(sp500$close[indexes]))
  volumes <- scale(Delt(sp500$volume[indexes]))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Removing the the series whose length is less than 1258 and the corresponding labels 

g <- numeric()

for (i in 1:505) {
  
  g[i] <- nrow(mts[[i]])
  
}

indexnull <- which(g != 1258)
indexnotnull <- which(g == 1258)


# Creating the 470 bivariate MTS 

mts <- list()

for (i in unique(sp500$Name)[indexnotnull]) {
  
  indexes <- which(sp500$Name == i)
  returns <- scale(Delt(sp500$close[indexes]))
  volumes <- scale(Delt(sp500$volume[indexes]))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

# Lets perform the clustering 

coherence2 <- listTomatrix(lapply(mts, quantile_coherence_re_im))
clustering <- kmeans(coherence2, 5)$cluster


# Series in cluster 1

index1 <- which(clustering == 3)
companies1 <- unique(sp500$Name)[indexnotnull][index1]

# Series in cluster 4 (apple cluster)

index4 <- which(clustering == 86)
companies4 <- unique(sp500$Name)[indexnotnull][index4]

index5 <- which(clustering == 7)
companies5 <- unique(sp500$Name)[indexnotnull][index5]
