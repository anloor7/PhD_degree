

# Lets create a function to compute the distance between two time series based on the copula-based
# periodogram. We are using levels 0.25, 0.5 and 0.75

# Input
# X, Y: two univariate time series 

cr_per_distance <- function(X, Y){
  
  levels <- c(0.25, 0.5, 0.75)
  CRX <- quantilePG(X, levels.1 = levels, type = "clipped")
  CRY <- quantilePG(Y, levels.1 = levels, type = "clipped")
  sCRX <- smoothedPG(CRX, weight = weight)
  sCRY <- smoothedPG(CRY, weight = weight)
  valuesx <- getValues(sCRX)
  valuesy <- getValues(sCRY)
  valuesx <-  as.vector(abs(valuesx))
  valuesy <- as.vector(abs(valuesy))
  
  return(EuclideanDistance(valuesx, valuesy))
  
}

a <- ts2(200)
b <- ts3(200)
cr_per_distance(a, b)

X <- arima.sim(n = 1000, list(ar = c(0.8), sd = sqrt(1)))
CRX <- quantilePG(X, levels.1 = levels, type = "qr")
sCRX <- smoothedPG(CRX, weight = weight)
plot(sCRX)



X <- arima.sim(n = 1000, list(ar = c(0.8), sd = sqrt(1)))
Y <- arima.sim(n = 1000, list(ma = c(0.8), sd = sqrt(1)))
n <- length(X)
freq <- (1:n/n)


levels <- c(0.25, 0.5, 0.75)
CRX <- quantilePG(X, levels.1 = levels, type = "clipped")
CRY <- quantilePG(Y, levels.1 = levels, type = "clipped")
valuesx <- getValues(CRX)
valuesy <- getValues(CRY)
valuesx <-  as.vector(abs(valuesx))
valuesy <- as.vector(abs(valuesy))
EuclideanDistance(valuesx, valuesy)

