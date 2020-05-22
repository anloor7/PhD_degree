

# Cluster analysis of dataset Libras via Quantile Coherence 

load('libras.RData')
load('gt_libras.RData')

# Performin k-means

coherences <- listTomatrix(lapply(libras, quantile_coherence))
clustering <- kmeans(coherences, 15)$cluster
external_validation(ground_truth_libras, clustering)


