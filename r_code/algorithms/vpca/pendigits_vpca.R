

########## PREPARATION OF DATASETS 


# We are going to load and prepare data  

# We start with the set of data called PenDigits 

pendigits <- read.csv('pendigits.txt', header = F)
colnames(pendigits) <- NULL


# We have to obtain 7494 MTS of length 8, each one with a label from 1 to 9 (there are 9 clusters)

S <- vector(mode = "list", length = nrow(pendigits))

# S is going to contain each MTS with its label

for (i in 1:length(S)) {
  # S[[i]] <- pendigits[i,] 
  S1 <- pendigits[i, seq(1, 17, by = 2)] 
  S2 <- pendigits[i, c(seq(2, 17, by = 2), 17)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

M <- vector(mode = "list", length = nrow(pendigits))

# M is going to contain each MTS without its label. 

for (i in 1:length(M)) {
  M[[i]] <- S[[i]][,seq(1, 8)]
}



# We can run the Fuzzy C-means algorithm by using M


# Obtaining both matrices V1 y V2, which are both 7494 x 8 matrices

V1 <- numeric()
V2 <- numeric()

for (i in 1:length(S)) {
  V1 <- rbind(V1, M[[i]][1,])
  V2 <- rbind(V2, M[[i]][2,])
}

# Dimensionality reduction via PCA

pc1 <- princomp(V1, cor = F, score = T)
summary(pc1)
pc2 <- princomp(V2, cor = F, score = T)
summary(pc2)

# In order to explain, in both instances, more than 90% of variability, we choose ps = 5 principal components 

# Computing the scores matrixes 

F1 <- pc1$scores[, 1:5]
F2 <- pc2$scores[, 1:5]

# Reconstructing the reduced sample of MTS time series

Y <- vector(mode = 'list', length = length(S))

for (i in 1 : length(S)) {
  Y[[i]] = matrix(nrow = 2, ncol = 5)
}

for (i in 1 : length(S)) {
  Y[[i]] = rbind(F1[i,], F2[i,])
}





