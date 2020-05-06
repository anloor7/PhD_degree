

# Lets check the approach generalising the Piccolos distance with p-values

# Loading the first 4 cluster from the Libras dataset 

libras <- read.csv('libras.txt', header = F)
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

# First, we construct a list of vectors containing the estimated coefficients of VAR models and their values for aic and
# bic criteria 

coefs <- list() 

time <- system.time(for (i in 1 : 96) {
  v1 <- VARMA(R[[i]], p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), v1$aic, v1$bic)
})

d <- listTomatrix(coefs)
clustering <- kmeans(d, 4)$cluster
external_validation(ground_truth, clustering, method = "adjusted_rand_index")
external_validation(fuzzytocrisp(t(fcm(d, 4)$u)), ground_truth)



# Lets see how this approach performs with some simulations

coefs <- list()

phi <- matrix(c(0.8, -0.4, 0.7, 0.6), nrow = 2)
the <- matrix(c(-0.5,0,0,-0.6), nrow = 2)

for (i in 1 : 100) {
  model <- VARMAsim(100, arlags = 1, phi = phi, sigma = diag(2))$series
  v1 <- VARMA(model, p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), v1$aic, v1$bic)
}

for (i in 101 : 200) {
  model <- VARMAsim(100, malags = 1, theta = the, sigma = 4*diag(2))$series
  v1 <- VARMA(model, p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), v1$aic, v1$bic)
}

for (i in 201 : 300) {
  model <- VARMAsim(100, arlags = 1, malags = 1, phi = phi, theta = the, sigma = diag(2))$series
  v1 <- VARMA(model, p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), v1$aic, v1$bic)
}

phi <- matrix(c(0.2, -0.7, 0.2, 0.3), nrow = 2)
for (i in 301 : 400) {
  model <- VARMAsim(150, arlags = 1, malags = 1, phi = phi, theta = the, sigma = diag(2))$series
  v1 <- VARMA(model, p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), v1$aic, v1$bic)
}



ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100))
d <- listTomatrix(coefs)
clustering <- kmeans(d, 3)$cluster
external_validation(ground_truth, clustering, method = "adjusted_rand_index")


# Testing dimensionality reduction 

varmalist <- list()

phi_c1 <- matrix(c(0.6, -0.4, 0.1, 0.0, 0.5, 0.5, -0.5, 0.1, 0.1, 0.3, 0.7, 0.5,
                   0.0, 0.1, 0.1, 0.6), nrow = 4)
sigma_c1 <- matrix(c(1, 0.25, 0.10, 0.0, 0.25, 1.0, 0.25, 0.1, 0.1, 0.25, 1, 0.25,
                     0.0, 0.1, 0.25, 1), nrow = 4)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7, 0.4, 0.7, -0.5, 0.3, 0.0, 0.5, 0.7), nrow = 4)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 0.25, 1, 0.25, 0.25, 0.1, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 1 ), nrow = 4)

for (i in 1 : 100) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 101 : 200) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, malags = 1, phi = phi_c1, theta = theta_c2, sigma = sigma_c2)$series
}

sigma <- list()

for (i in 1:200) {
  sigma[[i]] <- cov(varmalist[[i]])
}

s <- cpca(sigma, lambda = 0.70)
varmalist_reduced <- list()

for (i in 1 : 200) {
  varmalist_reduced[[i]] <- varmalist[[i]] %*% s
}

coefs <- list()
for (i in 1 : 200) {
  v1 <- VARMA(varmalist_reduced[[i]], p = 1)
  coefs[[i]] <- c(as.vector(v1$coef), as.vector(v1$Sigma), v1$aic, v1$bic)
}

d <- listTomatrix(coefs)
ground_truth <- c(rep(1, 300), rep(2, 300))
clustering <- kmeans(d, 2)$cluster
external_validation(ground_truth, clustering, method = "adjusted_rand_index")

external_validation(fuzzytocrisp(t(fcm(d, 2)$u)), ground_truth)


# Testing our approach versus Baragonas's second experiment


varmalist <- list()
phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7), nrow = 3)
sigma_c2 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)
phi_c3 <- matrix(c(0.5, -0.3, 0.0, -0.3, 0.5, -0.4, 0.0, 0.4, 0.5), nrow = 3)
theta_c3 <- matrix(c(-0.5, 0.4, 0.0, 0.4, -0.5, -0.3, 0.0, 0.3, -0.5), nrow = 3)
sigma_c3 <- matrix(c(1, 0.25, 0.1, 0.25, 1, 0.25, 0.10, 0.25, 1), nrow = 3)

for (i in 1 : 100) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 101 : 200) {
  varmalist[[i]] <- VARMAsim(100, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 201 : 300) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, malags = 1, phi = phi_c1, theta = theta_c3, sigma = sigma_c3)$series
}

sigma <- list()

for (i in 1:300) {
  sigma[[i]] <- cov(varmalist[[i]])
}

s <- cpca(sigma, lambda = 0.80)
varmalist_reduced <- list()

for (i in 1 : 300) {
  varmalist_reduced[[i]] <- varmalist[[i]] %*% s
}

coefs <- list()
for (i in 1 : 300) {
  v1 <- VARMA(varmalist_reduced[[i]], p = 1)
  coefs[[i]] <- c(as.numeric(v1$coef))#, v1$aic, v1$bic)
}

d <- listTomatrix(coefs)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))
clustering <- kmeans(d, 3)$cluster
external_validation(ground_truth, clustering, method = "adjusted_rand_index")

external_validation(fuzzytocrisp(t(fcm(d, 2)$u)), ground_truth)


# Testing our approach with series of 6 variables

varmalist <- list()
phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7, -0.5, 0, 0.5, 0.3, -0.4, 0.3, 0.5, -0.2, 0.8, 0, 
                   -0.5, 0.4, -0.3, 0, 0, 0.4, 0.2, 0.5, -0.7, 0.2, 0, 0.5, 0.5, 0.5, 0, -0.6, 0.2), nrow = 6)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 0.10, 0.25, 0.25, 1, 0.25, 0.10, 0.25, 0.10, 0.10, 0.25, 1, 0.10, 0, 0.25, 0.25,
                     0.10, 0.10, 1, 0.10, 0.25, 0.10, 0.25, 0, 0.10, 1, 0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 1), nrow = 6)
theta_c2 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7, 0.5, 0.6, 0.4, -0.6, 0, 0.2, 0.5, 0.4, 0.1, 0, 0.4,
                     0.4, 0.5, 0, 0.5, 0.2, 0, 0.7, 0.1, 0, 0.3, 0.4, 0.5, 0.6, 0, 0.2, 0.1), nrow = 6)
sigma_c2 <- sigma_c1
phi_c3 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.8, 0.5, 0.0, -0.3, 0.7, -0.5, 0, 0.5, 0.3, -0.4, -0.3, 0.5, -0.2, 0.8, 0, 
                   -0.5, 0.4, -0.3, 0, 0, -0.4, 0.2, 0.5, -0.7, 0.2, 0, 0.5, 0.5, -0.5, 0, -0.6, 0.2), nrow = 6)
theta_c3 <-  matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, -0.5, 0.0, 0.3, 0.7, 0.5, 0.6, 0.4, -0.6, 0, 0.2, 0.5, 0.4, 0.1, 0, 0.4,
                      0.4, 0.5, 0, 0.5, 0.2, 0, 0.7, 0.1, 0, 0.3, 0.4, 0.5, 0.6, 0, 0.2, 0.1), nrow = 6)
sigma_c3 <- 2*sigma_c2

set.seed(1234)
for (i in 1 : 100) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series
}

for (i in 101 : 200) {
  varmalist[[i]] <- VARMAsim(100, malags = 1, theta = theta_c2, sigma = sigma_c2)$series
}

for (i in 201 : 300) {
  varmalist[[i]] <- VARMAsim(100, arlags = 1, malags = 1, phi = phi_c3, theta = theta_c3, sigma = sigma_c3)$series
}


sigma <- list()

for (i in 1:300) {
  sigma[[i]] <- cov(varmalist[[i]])
}

s <- cpca(sigma, lambda = 0.60)
varmalist_reduced <- list()

for (i in 1 : 300) {
  varmalist_reduced[[i]] <- varmalist[[i]] %*% s
}

coefs <- list()
for (i in 1 : 300) {
  v1 <- VARMA(varmalist_reduced[[i]], p = 1)
  coefs[[i]] <- c(as.numeric(v1$coef), v1$aic, v1$bic)
}

d <- listTomatrix(coefs)
ground_truth <- c(rep(1, 100), rep(2, 100), rep(3, 100))
clustering <- kmeans(d, 3)$cluster
external_validation(ground_truth, clustering, method = "adjusted_rand_index")
intCriteria(d, clustering, 'xie_beni')


external_validation(fuzzytocrisp(t(fcm(d, 3)$u)), ground_truth)
