

# We are going to construct the boxplots for the paper 

# Loading the data 

title_size = 15.5
small_title_size = 14
text_size = 13.5
legend_size = 17

files_rdata <- list.files(pattern = "*.RData")
l_rdata <- length(files_rdata)

for (i in 1:l_rdata) {
  
  load(files_rdata[i])
  
}


# Boxplots for normal innovations 

# Boxplots for scenario 1, normal innovations


value <- c(ari_qc2_1_1000_n, larsen_qc2_1_1000_n, loo_qc2_1_1000_n, jaccard_qc2_1_1000_n,
           ari_w_1_1000_n, larsen_w_1_1000_n, loo_w_1_1000_n, jaccard_w_1_1000_n,
           ari_gcc_1_1000_n, larsen_gcc_1_1000_n, loo_gcc_1_1000_n, jaccard_gcc_1_1000_n,
           ari_kst_1_1000_n, larsen_kst_1_1000_n, loo_kst_1_1000_n, jaccard_kst_1_1000_n,
           ari_mah_1_1000_n, larsen_mah_1_1000_n, loo_mah_1_1000_n, jaccard_mah_1_1000_n,
           ari_pca_1_1000_n, larsen_pca_1_1000_n, loo_pca_1_1000_n, jaccard_pca_1_1000_n)

index <- c(rep(c(rep('ARI', 100), rep('LA', 100), rep('LOO1NN', 100), rep('JI',100)), 6))
metric <- c(rep('QCD', 400), rep('W', 400), rep('GCC', 400), rep('J', 400), rep('M', 400), rep('PCA', 400))

df <- data.frame(value = value, index = index, metric = metric)

n1 <- box_plot_paper_qcd(df)
plot_n1 <- n1 + ggtitle('Scenario 1') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              axis.text = element_text(size = text_size),
                                              strip.text.x = element_text(size = small_title_size))



# Boxplots for scenario 2, normal innovations 


value <- c(ari_qc2_2_1000_n, larsen_qc2_2_1000_n, loo_qc2_2_1000_n, jaccard_qc2_2_1000_n,
           ari_w_2_1000_n, larsen_w_2_1000_n, loo_w_2_1000_n, jaccard_w_2_1000_n,
           ari_gcc_2_1000_n, larsen_gcc_2_1000_n, loo_gcc_2_1000_n, jaccard_gcc_2_1000_n,
           ari_kst_2_1000_n, larsen_kst_2_1000_n, loo_kst_2_1000_n, jaccard_kst_2_1000_n,
           ari_mah_2_1000_n, larsen_mah_2_1000_n, loo_mah_2_1000_n, jaccard_mah_2_1000_n,
           ari_pca_2_1000_n, larsen_pca_2_1000_n, loo_pca_2_1000_n, jaccard_pca_2_1000_n)

index <- c(rep(c(rep('ARI', 100), rep('LA', 100), rep('LOO1NN', 100), rep('JI',100)), 6))
metric <- c(rep('QCD', 400), rep('W', 400), rep('GCC', 400), rep('J', 400), rep('M', 400), rep('PCA', 400))

df <- data.frame(value = value, index = index, metric = metric)

n2 <- box_plot_paper_qcd(df)
plot_n2 <- n2 + ggtitle('Scenario 2') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              axis.text = element_text(size = text_size),
                                              strip.text.x = element_text(size = small_title_size))


# Boxplots for scenario 3


value <- c(ari_qc2_3_1000_n, larsen_qc2_3_1000_n, loo_qc2_3_1000_n, jaccard_qc2_3_1000_n,
           ari_w_3_1000_n, larsen_w_3_1000_n, loo_w_3_1000_n, jaccard_w_3_1000_n,
           ari_gcc_3_1000_n, larsen_gcc_3_1000_n, loo_gcc_3_1000_n, jaccard_gcc_3_1000_n,
           ari_kst_3_1000_n, larsen_kst_3_1000_n, loo_kst_3_1000_n, jaccard_kst_3_1000_n,
           ari_mah_3_1000_n, larsen_mah_3_1000_n, loo_mah_3_1000_n, jaccard_mah_3_1000_n,
           ari_pca_3_1000_n, larsen_pca_3_1000_n, loo_pca_3_1000_n, jaccard_pca_3_1000_n)

index <- c(rep(c(rep('ARI', 100), rep('LA', 100), rep('LOO1NN', 100), rep('JI',100)), 6))
metric <- c(rep('QCD', 400), rep('W', 400), rep('GCC', 400), rep('J', 400), rep('M', 400), rep('PCA', 400))

df <- data.frame(value = value, index = index, metric = metric)


n3 <- box_plot_paper_qcd(df)
plot_n3 <- n3 + ggtitle('Scenario 3') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              axis.text = element_text(size = text_size),
                                              strip.text.x = element_text(size = small_title_size)) +
  coord_cartesian(ylim = c(-0.15, 1)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1))

grid.arrange(plot_n1, plot_n2, plot_n3)


























# Boxplots for student innovations


# Boxplots for scenario 1, student innovations


value <- c(ari_qc2_1_1000_t, larsen_qc2_1_1000_t, loo_qc2_1_1000_t, jaccard_qc2_1_1000_t,
           ari_w_1_1000_t, larsen_w_1_1000_t, loo_w_1_1000_t, jaccard_w_1_1000_t,
           ari_gcc_1_1000_t, larsen_gcc_1_1000_t, loo_gcc_1_1000_t, jaccard_gcc_1_1000_t,
           ari_kst_1_1000_t, larsen_kst_1_1000_t, loo_kst_1_1000_t, jaccard_kst_1_1000_t,
           ari_mah_1_1000_t, larsen_mah_1_1000_t, loo_mah_1_1000_t, jaccard_mah_1_1000_t,
           ari_pca_1_1000_n, larsen_pca_1_1000_n, loo_pca_1_1000_n, jaccard_pca_1_1000_n)

index <- c(rep(c(rep('ARI', 100), rep('LA', 100), rep('LOO1NN', 100), rep('JI',100)), 6))
metric <- c(rep('QCD', 400), rep('W', 400), rep('GCC', 400), rep('J', 400), rep('M', 400), rep('PCA', 400))

df <- data.frame(value = value, index = index, metric = metric)

t1 <- box_plot_paper_qcd(df)
plot_t1 <- t1 + ggtitle('Scenario 1') + theme(plot.title = element_text(hjust = 0.5, size = title_size),
                                              legend.text = element_text(size = legend_size),
                                              axis.text = element_text(size = text_size),
                                              strip.text.x = element_text(size = small_title_size))


# Boxplots for scenario 2, student innovations 

value <- c(ari_qc2_2_1000_t, larsen_qc2_2_1000_t, loo_qc2_2_1000_t, jaccard_qc2_2_1000_t,
           ari_w_2_1000_t, larsen_w_2_1000_t, loo_w_2_1000_t, jaccard_w_2_1000_t,
           ari_gcc_2_1000_t, larsen_gcc_2_1000_t, loo_gcc_2_1000_t, jaccard_gcc_2_1000_t,
           ari_kst_2_1000_t, larsen_kst_2_1000_t, loo_kst_2_1000_t, jaccard_kst_2_1000_t,
           ari_mah_2_1000_t, larsen_mah_2_1000_t, loo_mah_2_1000_t, jaccard_mah_2_1000_t,
           ari_pca_2_1000_t, larsen_pca_2_1000_t, loo_pca_2_1000_t, jaccard_pca_2_1000_t)

index <- c(rep(c(rep('ARI', 100), rep('LA', 100), rep('LOO1NN', 100), rep('JI',100)), 6))
metric <- c(rep('QCD', 400), rep('W', 400), rep('GCC', 400), rep('J', 400), rep('M', 400), rep('PCA', 400))

df <- data.frame(value = value, index = index, metric = metric)

t2 <- box_plot_paper_qcd(df)
plot_t2 <- t2 + ggtitle('Scenario 2') + theme(plot.title = element_text(hjust = 0.5, size = title_size,),
                                              legend.text = element_text(size = legend_size),
                                              axis.text = element_text(size = text_size),
                                              strip.text.x = element_text(size = small_title_size))

grid.arrange(plot_t1, plot_t2)
