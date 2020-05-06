

# Lets check the approach generalising the Piccolos distance with p-values

# Loading the first 4 cluster from the Libras dataset 

ground_truth <- c(rep(1, 24), rep(2, 24), rep(3, 24), rep(4, 24))

S <- vector(mode = "list", length = nrow(libras))

# S is going to contain each MTS with its label

for (i in 1 : length(S)) {
  # S[[i]] <- pendigits[i,] 
  S1 <- libras[i, c(seq(1, 90, by = 2), 91)] 
  S2 <- libras[i, c(seq(2, 90, by = 2), 91)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

save(S, file = 'libras.RData')

M <- vector(mode = "list", length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- S[[i]][,seq(1, 45)]
}


R <- list()


for (i in 1 : 96) {
  R[[i]] <- t(M[[i]])
}

# First, we construct a list of matrixes 2 x 6, the first row containing the estimated coefficients for 
# a VAR(1) fitted for each one of the time series, and the second containing their standard errors 

coefs <- list() 

time <- system.time(for (i in 1 : 96) {
  v1 <- VARMA(R[[i]], p = 1)
  coefs[[i]] <- c(v1$aic, v1$bic)
})

d <- listTomatrix(coefs)
clustering <- kmeans(d, 4)$cluster
external_validation(ground_truth, clustering, method = "jaccard_index")

