

# Performing proof with package quantspec according to its documentation (Kley)

library(quantspec)

# Toy example: eight independent standard gaussian random 
# variables 

Y <- rnorm(8)
bn <- qRegEstimator(Y, levels = c(0.25, 0.5, 0.75))
bn

# Retrieving the attributes frequencies and parallel of the object bn

getFrequencies(bn)
getParallel(bn)

# Get the values associated with certain frequencies or levels

getValues(bn, levels = c(0.25, 0.5))

# In the present case, B = 0 bootstrap replications were 
# performed 

# Computing and plotting the frequency representations d_{32, R}^tau,
# from 32 simulated, standard normally distributed random variables 

dn <- clippedFT(rnorm(32), levels = seq(0.05, 0.95, 0.05))
plot(dn, 2 * pi * (0:64)/32, levels = c(0.25, 05))


# The quantspec package includes three demos 

# demo('sp500')
# demo('wheatprices')
# demo('qar-simulation')

# Example: SP 500

library('zoo')

plot(sp500, xlab = 'time t', ylab = '', main = '')
acf(coredata(sp500), xlab = 'lag k', ylab = '', main = '')
acf(coredata(sp500)^2, xlab = 'lag k', ylab = '', main = '')


# Computing the CR periodogram for levels 0.05, 0.5 and 0.95, all Fourier frequencies and
# and with 250 bootstrap replications determined from a moving blocks bootstrap with block length
# l = 32


CR <- quantilePG(sp500, type = "clipped", levels.1 = c(0.05, 0.5, 0.95),
                 type.boot = "mbb", B = 250, l = 32)
freq <- getFrequencies(CR)
plot(CR, frequencies = freq[freq > 0 & freq <= pi],
     ylab = expression({I[list(n, R)]^{list(tau[1], tau[2])}}(omega)))


# Simulation study: analyzing a quantile autoregressive process 

# Generating the graphical representation of a copula spectral density regarding a given QAR(1) process

csd <- quantileSD(N = 2^9, seed.init = 2581, type = "copula",
                    ts = ts1, levels.1 = c(0.25, 0.5, 0.75), R = 100, quiet = TRUE)
plot(csd, ylab = expression(f[list(q[tau[1]], q[tau[2]])](omega)))

# Computing the copula spectral density with high precision (takes time)

csd <- quantileSD(N = 2^12, seed.init = 2581, type = "copula", ts = ts1, 
                  levels.1 = c(0.25, 0.5, 0.75), R = 5000)
save(csd, file = "csd-qar1.rdata")
load("csd-qar1.rdata")
plot(csd, frequencies = 2 * pi * (1:2^8) / 2^9,
     ylab = expression(f[list(q[tau[1]], q[tau[2]])](omega))) # The parameter frequencies was used to compare this
# plot with the previous one 

# To get a first idea of how well the estimator performs, plot the smoothed CR periodogram computed for one 
# simulated QAR(1) time series of length 512 

sCR <- smoothedPG(ts1(512), levels.1 = c(0.25, 0.5, 0.75), weight = kernelWeight(W = W1, bw = 0.1))
plot(sCR, qsd = csd,
     ylab = bquote(paste(hat(G)[list(n, R)](list(tau[1], tau[2], omega)), 
                         " and ", f[list(q[tau[1]], q[tau[2]])](omega))))

# Simulation study. 5000 independent QAR(1) time series 

set.seed(2581)
ts <- ts1
N <- 128
R <- 5000
freq<-2 *pi*(1:16)/32
levels <- c(0.25, 0.5, 0.75) 
J <- length(freq)
K <- length(levels)
sims <- array(0, dim = c(4, R, J, K, K))
weight <- kernelWeight(W = W1, bw = 0.3) # Epanechnikov kernel and rather large bandwidth 

# Simulation 

for (i in 1 : R) {
Y <- ts(N)
CR <- quantilePG(Y, levels.1 = levels, type = "clipped")
LP <- quantilePG(Y, levels.1 = levels, type = "qr")
sCR <- smoothedPG(CR, weight = weight)
sLP <- smoothedPG(LP, weight = weight)
sims[1, i, , , ] <- getValues(CR, frequencies = freq)[, , , 1]
sims[2, i, , , ] <- getValues(LP, frequencies = freq)[, , , 1]
sims[3, i, , , ] <- getValues(sCR, frequencies = freq)[, , , 1]
sims[4, i, , , ] <- getValues(sLP, frequencies = freq)[, , , 1] }


# The true copula spectral density is copied to an array trueV.
# Using the arrays sims and trueV the root integrated mean squared errors are computed as follows

trueV <- getValues(csd, frequencies = freq) 
SqDev <- array(apply(sims, c(1, 2), 
                     function(x) {abs(x - trueV)^2}), dim = c(J, K, K, 4, R)) 
rimse <- sqrt(apply(SqDev, c(2, 3, 4), mean))
rimse

# We can see that the smoothed quantile periodogram posses smaller root integrated mean squared errors than 
# the quantile periodograms














