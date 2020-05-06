

ecg <- read.table('ecg.txt', header = F, sep = '')
ecg <- read.table('ecg1.txt', header = F, sep = '')

ecg <- data.matrix(ecg)
colnames(ecg) <- NULL
colnames(ecg1) <- NULL

# Lets perform model based clustering on the training data

coefs <- matrix(nrow = 100, ncol = 3)
ground_truth <- ecg[,1]
x <- ecg[,2 : ncol(ecg)]

for (i in 1 : 100) {
  v1 <- arima(x[i,], order = c(1, 0, 0))
  coefs[i,] <- c(as.vector(v1$coef), v1$aic)
}

clustering <- kmeans(coefs, 2)$cluster
external_validation(ground_truth, clustering, method = "jaccard_index")
