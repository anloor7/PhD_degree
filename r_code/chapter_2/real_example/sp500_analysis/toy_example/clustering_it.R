

sp500 <- read.csv('all_stocks_5yr.csv')
head(sp500)
sector <- read.csv('sp5_sector.csv')
head(sector)

symbols_it <- sector$Symbol[sector$Sector == 'Information Technology']
names_it <- sector$Name[sector$Sector == 'Information Technology']

# Creating the 505 bivariate MTS 

mts <- list()

for (i in symbols_it) {
  
  indexes <- which(sp500$Name == i)
  returns <- scale(Delt(sp500$close[indexes]))
  volumes <- scale(Delt(sp500$volume[indexes]))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}


g <- numeric()

for (i in 1:71) {
  
  g[i] <- nrow(mts[[i]])
  
}

indexnull <- which(g != 1258)
indexnotnull <- which(g == 1258)


# Creating bivariate MTS 

mts <- list()

for (i in symbols_it[indexnotnull]) {
  
  indexes <- which(sp500$Name == i)
  returns <- scale(Delt(sp500$close[indexes]))
  volumes <- scale(Delt(sp500$volume[indexes]))
  returns <- returns[-1]
  volumes <- volumes[-1]
  mts[[i]] <- cbind(returns, volumes)
  
}

symbols_clustering <- symbols_it[indexnotnull]
names_clustering <- names_it[indexnotnull]

# Lets perform the clustering 

K <- 10
coherence2 <- listTomatrix(lapply(mts, quantile_coherence_re_im))
clustering <- pam(coherence2, K)$cluster


# Displaying clustering soluction

clusters <- list()

for (i in 1:K) {
  
  clusters[[i]] <- names_clustering[clustering == i]
  
}
