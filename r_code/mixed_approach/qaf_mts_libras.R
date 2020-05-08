

# Lets perfom clustering via multivariate QAF with Libras dataset


# Loading and preparing the data 

libras <- read.csv('libras.txt', header = F)
colnames(libras) <- NULL
S <- vector(mode = "list", length = nrow(libras))

for (i in 1 : length(S)) {
  S1 <- libras[i, c(seq(1, 90, by = 2), 91)] 
  S2 <- libras[i, c(seq(2, 90, by = 2), 91)]
  S[[i]] <- rbind(as.numeric(S1), as.numeric(S2))
}

M <- vector(mode = "list", length = length(S))
for (i in 1:length(S)) {
  M[[i]] <- t(S[[i]][,seq(1, 45)])
}

# Applying k-means algorithm 

M <- lis
gamma <- listTomatrix(lapply(M, qaf_mts_coefs_xy))
