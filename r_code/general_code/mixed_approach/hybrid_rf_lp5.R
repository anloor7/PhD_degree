

# Applying k-means clustering with QAF and DTW distances to the dataset LP5

# Loading the dataset 

lp5 <- read.csv('lp5.data.txt', header = F, sep = '\t')
head(lp5)


# We have to obtain  164 MTS of length 15, each one with a given label (there are 5 labels)

S <- vector(mode = "list", length = 164)

for (i in 0 : 163) {
  S[[i + 1]] <- matrix(nrow = 16, ncol = 7)
  S[[i + 1]] <- lp5[ (16 * i + 1): (16 * (i + 1)),]
}

for (i in 1 : 164) {
  label <- as.numeric(S[[i]][1, 1])
  S[[i]] <- cbind(S[[i]], rep(label, 16))
  colnames(S[[i]])[8] <- 'label'
  S[[i]] <- S[[i]][2:16, 2:8]
  S[[i]] <- t(S[[i]])
  rownames(S[[i]]) <- NULL
  colnames(S[[i]]) <- NULL
}

# Removing labels

M <- vector(mode = "list", length = 164)
for (i in 1 : 164) {
  M[[i]] <- matrix(nrow = 6, ncol = 15)
  M[[i]] <- S[[i]][1:6,]
}

# Ground truth 

ground_truth <- numeric(164)

for (j in 1: 164) {
  ground_truth[j] <- S[[j]][7, 1]
}


# Hybrid approach with k-means

hybrid_distance <- function(X, Y){
  
  a <- as.vector(qaf_mts_coefs_xy(X))
  b <- as.vector(qaf_mts_coefs_xy(Y))
  c <- EuclideanDistance(a, b)
  d <- dtw_mts(X, Y)
  c
}

clustering <- km_mts(M, 5, dis = hybrid_distance)
external_validation(ground_truth, clustering)