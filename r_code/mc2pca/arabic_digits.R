

# Loading the dataset standing for test set of spoken Arabic Digits

arabic <- read.csv('test_arabic_digits.txt', header = F, sep = ' ')
column1 <- arabic$V1
indexes <- which(column1 == 'a')

# We have to obtain 2200 series with 13 variables each one, and different lengths 

n <- 2200
X <- list(n)

for (i in 1 : (n-1)) {
  
  X[[i]] <- arabic[indexes[i] : indexes[i + 1],]
}
X[[n]] <- arabic[(indexes[n]) : nrow(arabic),]

# Removing useless lines 

for (i in 1 : n) {
  useless <- which(X[[i]]$V1 == 'a')
  X[[i]] <- X[[i]][-useless,]
  rownames(X[[i]]) <- NULL
  colnames(X[[i]]) <- NULL
  X[[i]] <- apply(as.matrix(X[[i]]), c(1,2), as.numeric)
}

ground_truth <- numeric(2200)

for (i in 1 : 10) {
  ground_truth[ ((i - 1)*220 + 1) : (220*i)] <- i
}


# Lets perform mc2pca over the dataset and evaluate the results 

clustering <- mc2pca(X, 10, lambda = 0.90, niter = 300, tol = 0.01)

external_validation(clustering, ground_truth, summary_stats = T)







