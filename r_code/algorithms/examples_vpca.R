

########## EN ESTE CÓDIGO, VAMOS A APLICAR EL MÉTODO DE VPCA (HE, 2020) A CONSJUNTOS SIMULADOS DE SERIES DE 
########## TIEMPO MULTIVARIANTES 


# Obtención de los conjuntos de series de tiempo multivariantes


# Vamos a simular dos conjuntos de series multivariantes,
# cada uno de los cuales contendrá 10 elementos. Es decir, tendremos K = 2 clusters verdaderos con 10 series cada uno, M = 20. El
# número de variables de cada una de las series será de N = 2 (series bivariantes). El número de instantes temporales considera
# dos será de L = 20

library(MTS) # Es un paquete para simular modelos VARMA 

# Definimos las matrices en las que serán almacenadas las series. Las matrices M1,...,M10 constituirán el primer cluster, 
# mientras que las matrices S1, ..., S10, el segundo cluster 

M <- vector(mode = 'list', length= 10)
for (i in 1:10) {
  M[[i]] = matrix(nrow = 20, ncol = 2)
}

P <- vector(mode = 'list', length= 10)
for (i in 1:10) {
  P[[i]] = matrix(nrow = 20, ncol = 2)
}


# Definimos los parámetros de las series del primer cluster 

p1 = matrix(c(0.2, -0.6, 0.3, 1), 2, 2)
sig1 = matrix(c(0.2, 0.15, 0.15, 0.2), 2, 2)
th1 = matrix(c(-0.5, 0, 0, -0.6), 2, 2)

# Simulamos las series del primer cluster 

# set.seed(2000)
for (i in 1:10) {
 # M[[i]] <- t(VARMAsim(20, arlags = c(1), malags = c(1), phi = p1, 
             #      theta = th1, sigma = sig1)$series)
 M[[i]] <- t(mvrnorm(20, mu = c(0.75, 0.75), Sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)))
}

# Definimos los parámetros de las series del segundo cluster 

p2 = matrix(c(0.2, 0.6, -0.3, 1), 2, 2)
sig2 = matrix(c(7, 0.35, 0.35, 4), 2, 2)
th2 = matrix(c(-0.1, 0, 0, 0.8), 2, 2)

# Simulamos las series del segundo cluster 


for (i in 1:10) {
  # S[[i]] <- t(VARMAsim(20, arlags = c(1), malags = c(1), phi = p2, 
    #                 theta = th2, sigma = sig2)$series)
  
 P[[i]] <- t(mvrnorm(20, mu = c(-0.75, -0.75), Sigma = matrix(c(1, -0.5, -0.5, 1), nrow = 2)))
}

# Vamos a representar conjuntamente, las diez series relativas a cada uno de los clusters. Para las representaciones vamos a usar
#, de aquí en adelante para los que queda de tesis, el paquete ggplot2, que proporciona unas representaciones visualmente muy 
# atractivas para series de tiempo 

library(ggplot2)
library(gbm)
library(grid)

# Representación de las series del primer cluster 

x1 <- seq(1, 20)
y1 <- vector(mode = 'list', length = 10)
x2 <- x1
y2 <- vector(mode = 'list', length = 10)
plotl1 <- vector(mode = 'list', length = 10)
auxiliary <- numeric()
for (i in 1:10) {
  y1 = M[[i]][1,]
  y2 = M[[i]][2,]
  df1 <- data.frame(x1 = x1, y1 = M[[i]][1,], x2 = x2, y2 = M[[i]][2,])
 auxiliary <- ggplot(df1, aes(x1 = x1, y1 = y1)) + geom_line(aes(x1, y1), col = 'red') 
 plotl1[[i]] <- auxiliary + geom_line(aes(x2, y2), col = 'blue')
}

grid.arrange(plotl1[[1]], plotl1[[2]], plotl1[[3]], plotl1[[4]], 
             plotl1[[5]], plotl1[[6]], plotl1[[7]], plotl1[[8]],
             plotl1[[9]], plotl1[[10]], top = textGrob("Primer cluster"))

# Representación de las series del segundo cluster

x1 <- seq(1, 20)
y1 <- vector(mode = 'list', length = 10)
x2 <- x1
y2 <- vector(mode = 'list', length = 10)
plotl2 <- vector(mode = 'list', length = 10)

for (i in 1:10) {
  y1 = P[[i]][1,]
  y2 = P[[i]][2,]
  df1 <- data.frame(x1 = x1, y1 = P[[i]][1,], x2 = x2, y2 = P[[i]][2,])
  plotl2[[i]] <- ggplot(df1, aes(x1 = x1, y1 = y1)) + geom_line(aes(x1, y1), col = 'red') +
    geom_line(aes(x2, y2), col = 'blue')
}

grid.arrange(plotl2[[1]], plotl2[[2]], plotl2[[3]], plotl2[[4]], 
             plotl2[[5]], plotl2[[6]], plotl2[[7]], plotl2[[8]],
             plotl2[[9]], plotl2[[10]], top = textGrob('Segundo cluster'))

# Ya tenemos definidos y representados los clusters de antemano. Ahora, vamos a intentar resolver el problema aplicando el algoritmo 
# paso por paso. El hiperparámetro K = 2 (dos clusters) va a estar fijado de antemano en todo momento 

# Calculamos las matrices V1 y V2, que serán ambas, matrices 20 x 20

V1 <- numeric()
V2 <- numeric()

for (i in 1:10) {
  V1 <- rbind(V1, M[[i]][1,])
  V2 <- rbind(V2, M[[i]][2,])
}

for (i in 1:10) {
  V1 <- rbind(V1, P[[i]][1,])
  V2 <- rbind(V2, P[[i]][2,])
}

# Calculamos las matrices de covarianzas 

sigma1 <- cov(V1)
sigma2 <- cov(V2)

# Reducimos la dimensión de las matrices V1 y V2 mediante componentes principales

pc1 <- princomp(V1, cor = T, score = T)
summary(pc1)
pc2 <- princomp(V2, cor = T, score = T)
summary(pc2)

# Tomando 5 componentes principales, podemos explicar, en ambos casos, más del 90% de la variabilidad de nuestro conjunto de
# datos. En consecuencia, construimos las matrices F1 y F2

F1 <- pc1$scores[, 1:5]
F2 <- pc2$scores[, 1:5]

# Recostruimos la muestra reducida de series de tiempo multivariantes. Y1, ..., Y20 son las series reducidas 

Y <- vector(mode = 'list', length = 20)

for (i in 1:20) {
  Y[[i]] = matrix(nrow = 2, ncol = 5)
}

for (i in 1:20) {
  Y[[i]] = rbind(F1[i,], F2[i,])
}

# Ahora, vamos a calcular la matriz S, que en nuestro caso es una matriz 10 x 10

S <- matrix(nrow = 10, ncol = 10)
constant <- 1/(2*pi*(1-(1/2))^2)
coefficients_previous <- numeric(100) # Vector con tantos elementos como entradas tiene la matriz S


# Primero fila de S

a <- 1
  for (p1 in 1:5) {
    for (n1 in 1:2) {
      for (n2 in 1:2) {
        for (p2 in 1:5){
          coefficients_previous[a] <- sqrt((n1-n2)^2 + (p1 - p2)^2)
          a <- a + 1
        }
          
      }
    }
  }
 
S_previous <- matrix(coefficients_previous, nrow = 10, ncol = 10, byrow = T) # La matriz ya casi está de la forma que nos interesa,
# excepto por el orden de las columnas 


S_previous <- S_previous[,c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]
k <- 2*(1-1/5)^2
S <- 1/(pi * k) * exp(-S_previous^2/k)

# Calculamos la matriz de distancias entre las 20 series iniciales 
# Estiramos cada serie 

y <- vector(mode = 'list', length = 10)
for (i in 1:length(Y)) {
  y[[i]] <- as.vector(Y[[i]])
  
}

distances <- matrix(nrow = 20, ncol = 20)
for (i in 1:20) {
  for (j in 1:20) {
    distances[i, j] <- sqrt((y[[i]]-y[[j]]) %*% S %*% (y[[i]]-y[[j]]))
  }
  
}

# Vamos a llevar a llevar a cabo un clustering jerárquico para esta matriz de distancias inicial 

clustering <- hclust(dist(distances))
plot(clustering)
cutree(clustering, 2) # Nos quedamos con dos clusters

# El algoritmo discrimina perfectamente aplicando un clustering jerárquico cuando los dos grupos se diferencian claramente 

# Vamos a implementar, a continuación, Fuzzy c-means 

K <- 2 # Hiperparámetro. Número de clusters 

# Seleccionamos aleatoriamente K centroides 

vector_centroids <- sample(20, 2, replace = F)


# Calculamos la matriz de distancias 

# Para ello, vamos a construir, en primer lugar, una función que calcule la distancia que estamos considerando

dist <- function(x, y){
  sqrt((x-y) %*% S %*% (x-y))
}

dist_matrix <- matrix(nrow = K, ncol = 20)
for (i in 1:20) {
  for (k in vector_centroids) {
    dist_matrix[which(vector_centroids == k), i] <- dist(y[[i]], y[[k]])
  }
  
}

rownames(dist_matrix) <- vector_centroids

# Una vez tenemos la matriz de distancias para comenzar, vamos a realizar un determinado número de iteraciones en las que
# calcularemos la membership matrix y actualizaremos el centro de los clusters 

niter <- 100
mem_matrix <- matrix(nrow = k, ncol = 20)

for (i in 1:niter) {
  
}
