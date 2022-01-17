

# Loading auxiliary functions, data and packages

library(quantspec) # Loading the quantspec package
library(ggplot2) # Loading the library ggplot2
library(fclust) # Loading the library fclust

# Function to compute the QCD-based features

quantile_quantities_re_im <- function(X, levels = seq(0.1, 0.9, 0.4)) {

  c <- ncol(X)
  d <- length(levels)
  qSPG <- smoothedPG(X, levels.1 = levels, type = 'clipped')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)


  matrix <- matrix(qSPGv, ncol = 1)
  coherence <- c(Re(matrix), Im(matrix))

}

# Function to transform a list of features to a matrix (each feature vector is a row in the matrix)

listTomatrix <- function(l){
  n <- length(l)
  s <- length(l[[1]])
  m <- matrix(nrow = n, ncol = s)
  for (i in 1 : n) {
    m[i,] <- l[[i]]
  }
  m
}



load('mts_50.RData') # Loading the 50 first companies of S&P 500 index
load('vector_names_50.RData') # Loading the corresponding names


mts <- mts_50
vector_names <- vector_names_50


# Paypal starts on 2015-07-06, so we are going to consider all series from that date

mts_reduced <- list()
start <- 1258 - 654 + 1

for (i in 1 : 50) {

  if (i != 47) {

    mts_reduced[[i]] <- mts[[i]][start : 1258 ,]

  } else {

    mts_reduced[[i]] <- mts[[i]]

  }

}


 # Retaining the top 20 companies

vector_names_top_20 = c('AAPL', 'MSFT', 'AMZN', 'GOOGL', 'GOOG',
                 'FB', 'TSLA', 'BRK.B', 'V', 'JNJ',
                 'WMT', 'JPM', 'MA', 'PG', 'UNH',
                 'DIS', 'NVDA', 'HD', 'PYPL', 'BAC')

mts_top_20 = list()

for (i in 1 : 20) {

  index = which(vector_names_top_20[i] == vector_names)
  mts_top_20[[i]] = mts_reduced[[index]]

}


# Computing the dissimilarity matrix based on d_{QCD} and PCA

coherence2 <- listTomatrix(lapply(mts_top_20, quantile_quantities_re_im)) # QCD-features
p_comp <- prcomp(coherence2) # PCA
n_p_comp <- ceiling(0.12*min(dim(p_comp$x))) # Selection of principal components
p_comp_matrix <- data.frame(p_comp$x[, 1 : n_p_comp])
dis_matrix = proxy::dist(p_comp_matrix) # Computation of distance matrix


# Representing a 2-dimensional sccaling plot of the companies

mds <- cmdscale(dis_matrix, 2, list. = T)
mds$GOF
df <- data.frame(cbind(mds$points), vector_names_top_20)

# Distance matrix for 2d reduced space

plot_mds <- ggplot(df, aes(x = X1, y = X2)) + geom_point(size = 2,
                                                         col = 'blue') +
  xlab('Coordinate 1') + ylab('Coordinate 2') +
  geom_text(aes(label = vector_names_top_20), size = 4,
            # position = position_jitter(width = 0, height = 0.0),
            nudge_y = c(rep(0.01, 18), -0.01, 0.01),
            nudge_x = 0.00,
            check_overlap = F) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5,  size = 18),
        legend.position = '') + ggtitle('Top 20 companies of S&P 500')


# Performing fuzzy C-means for C = 6 and m = 1.8

m <- 4
clustering_pca <- FKM(p_comp_matrix, k = 6, m = 1.8, RS = 300)
clustering_pca$k
table_fuzzy <- round(clustering_pca$U, 4)
rownames(table_fuzzy) <- vector_names_top_20
colnames(table_fuzzy) <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6')
table_fuzzy # Membership matrix

