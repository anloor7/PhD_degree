

# We are going to construct the barplots for the paper 

# Loading the data 

title_size = 15.5
small_title_size = 14
text_size = 13.5
legend_title_size = title_size - 2
x_size = text_size + 1

files_rdata <- list.files(pattern = "*.RData")
l_rdata <- length(files_rdata)

for (i in 1:l_rdata) {
  
  load(files_rdata[i])
  
}



# Barplots for normal innovations



# Boxplots for scenario 1, normal innovations

howmany_qc2 <- add_5c(count(howmany_clusters(clustering_qc2_1_1000_n)))
howmany_w <- add_5c(count(howmany_clusters(clustering_w_1_1000_n)))
howmany_gcc <- add_5c(count(howmany_clusters(clustering_gcc_1_1000_n)))
howmany_kst <- add_5c(count(howmany_clusters(clustering_kst_1_1000_n)))
howmany_mah <- add_5c(count(howmany_clusters(clustering_mah_1_1000_n)))
howmany_pca <- add_5c(count(howmany_clusters(clustering_pca_1_1000_n)))

df_a <- rbind(howmany_qc2, howmany_w, howmany_gcc, howmany_kst, howmany_mah, howmany_pca)
metric <- c(rep('QCD', nrow(howmany_qc2)), rep('W', nrow(howmany_w)), rep('GCC', nrow(howmany_gcc)),
            rep('J', nrow(howmany_kst)), rep('M', nrow(howmany_mah)), rep('PCA', nrow(howmany_pca)))

df <- cbind(df_a, metric)

n1 <- barplot_nclusters(df)
plot_n1 <- n1 + ggtitle('Scenario 1') + theme(plot.title = element_text(hjust = 0.5, size = title_size), 
                                              legend.text = element_text(size = legend_size),
                                              legend.title = element_text(size = legend_title_size),
                                              axis.text = element_text(size = text_size))


# Boxplots for scenario 2, normal innovations 

howmany_qc2 <- add_4c(count(howmany_clusters(clustering_qc2_2_1000_n)))
howmany_w <- add_4c(count(howmany_clusters(clustering_w_2_1000_n)))
howmany_gcc <- add_4c(count(howmany_clusters(clustering_gcc_2_1000_n)))
howmany_kst <- add_4c(count(howmany_clusters(clustering_kst_2_1000_n)))
howmany_mah <- add_4c(count(howmany_clusters(clustering_mah_2_1000_n)))
howmany_pca <- add_4c(count(howmany_clusters(clustering_pca_2_1000_n)))

df_a <- rbind(howmany_qc2, howmany_w, howmany_gcc, howmany_kst, howmany_mah, howmany_pca)
metric <- c(rep('QCD', nrow(howmany_qc2)), rep('W', nrow(howmany_w)), rep('GCC', nrow(howmany_gcc)),
            rep('J', nrow(howmany_kst)), rep('M', nrow(howmany_mah)), rep('PCA', nrow(howmany_pca)))

df <- cbind(df_a, metric)

n2 <- barplot_nclusters(df)
plot_n2 <- n2 + ggtitle('Scenario 2') + theme(plot.title = element_text(hjust = 0.5,  size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              legend.title = element_text(size = legend_title_size),
                                              axis.text = element_text(size = text_size))


# Barplots for scenario 3, normal innovations


howmany_qc2 <- add_3c(count(howmany_clusters(clustering_qc2_3_1000_n)))
howmany_w <- add_3c(count(howmany_clusters(clustering_w_3_1000_n)))
howmany_gcc <- add_3c(count(howmany_clusters(clustering_gcc_3_1000_n)))
howmany_kst <- add_3c(count(howmany_clusters(clustering_kst_3_1000_n)))
howmany_mah <- add_3c(count(howmany_clusters(clustering_mah_3_1000_n)))
howmany_pca <- add_3c(count(howmany_clusters(clustering_pca_3_1000_n)))

df_a <- rbind(howmany_qc2, howmany_w, howmany_gcc, howmany_kst, howmany_mah, howmany_pca)
metric <- c(rep('QCD', nrow(howmany_qc2)), rep('W', nrow(howmany_w)), rep('GCC', nrow(howmany_gcc)),
            rep('J', nrow(howmany_kst)), rep('M', nrow(howmany_mah)), rep('PCA', nrow(howmany_pca)))

df <- cbind(df_a, metric)

barplot_nclusters(df)

n3 <- barplot_nclusters(df)
plot_n3 <- n3 + ggtitle('Scenario 3') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              legend.title = element_text(size = legend_title_size),
                                              axis.text = element_text(size = text_size))


grid.arrange(plot_n1, plot_n2, plot_n3)








# Barplots for student innovations



# Barplots for scenario 1, student innovations

howmany_qc2 <- add_5c(count(howmany_clusters(clustering_qc2_1_1000_t)))
howmany_w <- add_5c(count(howmany_clusters(clustering_w_1_1000_t)))
howmany_gcc <- add_5c(count(howmany_clusters(clustering_gcc_1_1000_t)))
howmany_kst <- add_5c(count(howmany_clusters(clustering_kst_1_1000_t)))
howmany_mah <- add_5c(count(howmany_clusters(clustering_mah_1_1000_t)))
howmany_pca <- add_5c(count(howmany_clusters(clustering_pca_1_1000_t)))

df_a <- rbind(howmany_qc2, howmany_w, howmany_gcc, howmany_kst, howmany_mah, howmany_pca)
metric <- c(rep('QCD', nrow(howmany_qc2)), rep('W', nrow(howmany_w)), rep('GCC', nrow(howmany_gcc)),
            rep('J', nrow(howmany_kst)), rep('M', nrow(howmany_mah)), rep('PCA', nrow(howmany_pca)))

df <- cbind(df_a, metric)

t1 <- barplot_nclusters(df)
plot_t1 <- t1 + ggtitle('Scenario 1') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              legend.title = element_text(size = legend_title_size),
                                              axis.text = element_text(size = text_size))


# Barplots for scenario 2, student innovations 

howmany_qc2 <- add_4c(count(howmany_clusters(clustering_qc2_2_1000_t)))
howmany_w <- add_4c(count(howmany_clusters(clustering_w_2_1000_t)))
howmany_gcc <- add_4c(count(howmany_clusters(clustering_gcc_2_1000_t)))
howmany_kst <- add_4c(count(howmany_clusters(clustering_kst_2_1000_t)))
howmany_mah <- add_4c(count(howmany_clusters(clustering_mah_2_1000_t)))
howmany_pca <- add_4c(count(howmany_clusters(clustering_pca_2_1000_t)))

df_a <- rbind(howmany_qc2, howmany_w, howmany_gcc, howmany_kst, howmany_mah, howmany_pca)
metric <- c(rep('QCD', nrow(howmany_qc2)), rep('W', nrow(howmany_w)), rep('GCC', nrow(howmany_gcc)),
            rep('J', nrow(howmany_kst)), rep('M', nrow(howmany_mah)), rep('PCA', nrow(howmany_pca)))

df <- cbind(df_a, metric)

t2 <- barplot_nclusters(df)
plot_t2 <- t2 + ggtitle('Scenario 2') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              legend.title = element_text(size = legend_title_size),
                                              axis.text = element_text(size = text_size))

grid.arrange(plot_t1, plot_t2)

