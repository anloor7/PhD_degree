

# First comparison: boxplot and times


# Box plot

# Loading the data 

file_names <- as.list(dir(pattern = "boxplot_*"))
lapply(file_names,load,.GlobalEnv)



# Making the graph 

values <- c(boxplot_qaf_pam, boxplot_qaf_km, boxplot_qc1_pam, boxplot_qc1_km, boxplot_qc2_pam, boxplot_qc2_km, boxplot_dtw1_pam,
            boxplot_dtw2_pam, boxplot_pdc_pam, boxplot_kst_pam, boxplot_kst_km, boxplot_pca_pam, boxplot_pca_km, 
            boxplot_mc2pca_km)
nam <- c('QAF PAM', 'QAF KM', 'QC1 PAM', 'QC1 KM', 'QC2 PAM', 'QC2 KM', 'DTW1 PAM', 'DTW2 PAM', 'PDC PAM', 'KST PAM',
         'KST KM', 'PCA PAM', 'PCA KM', 'MC2PCA KM')

df <- data.frame(x = rep(0, 280), y = values)

a <- seq(1, 280, by = 20)
for (i in a) {
  df$x[i : (i + 19)] <- rep(nam[which(a == i)], 20)
}

p <- ggplot(df, aes(x = x, y = y, fill = x))
p + geom_boxplot() + xlab('Method') + ylab('ARI') +
  theme(axis.text.x = element_blank(), legend.title = element_blank(),
        axis.ticks.x = element_blank()) 

# Time 

qaf_pam_time <- 3.25
qaf_km_time <- 1.12
qc1_pam_time <- 37.02
qc1_km_time <- 13.44
qc2_pam_time <- 38.10
qc2_km_time <- 16.69
dtw1_pam_time <- 736.8
dtw2_pam_time <- 90
pdc_pam_time <- 12.96
kst_pam_time <- 513
kst_km_time <- 538
pca_pam_time <- 3.33
pca_km_time <- 3.29
mc2pca_km_time <- 73.14 


# Second comparison: boxplot and times

# Box plot

# Loading the data 

file_names <- as.list(dir(pattern = "*he"))
lapply(file_names,load,.GlobalEnv)



# Making the graph 

values <- c(boxplot_qaf_pam_he, boxplot_qaf_km_he, boxplot_qc1_pam_he, boxplot_qc1_km_he, boxplot_qc2_pam_he, boxplot_qc2_km_he, boxplot_dtw1_pam_he,
            boxplot_dtw2_pam_he, boxplot_pdc_pam_he, boxplot_kst_pam_he, boxplot_kst_km_he, boxplot_pca_pam_he, boxplot_pca_km_he, 
            boxplot_mc2pca_km_he)
nam <- c('QAF PAM', 'QAF KM', 'QC1 PAM', 'QC1 KM', 'QC2 PAM', 'QC2 KM', 'DTW1 PAM', 'DTW2 PAM', 'PDC PAM', 'KST PAM',
         'KST KM', 'PCA PAM', 'PCA KM', 'MC2PCA KM')

df <- data.frame(x = rep(0, 280), y = values)

a <- seq(1, 280, by = 20)
for (i in a) {
  df$x[i : (i + 19)] <- rep(nam[which(a == i)], 20)
}

p <- ggplot(df, aes(x = x, y = y, fill = x))
p + geom_boxplot() + xlab('Method') + ylab('ARI') +
  theme(axis.text.x = element_blank(), legend.title = element_blank(),
        axis.ticks.x = element_blank()) 





# Time 


qaf_pam_time_he <- 1.67
qaf_km_time_he <- 0.58
qc1_pam_time_he <- 32.76
qc1_km_time_he <- 12.75
qc2_pam_time_he <- 31.95
qc2_km_time_he <- 10.24
dtw1_pam_time_he <- 537
dtw2_pam_time_he <- 36.98
pdc_pam_time_he <- 4.46
kst_pam_time_he <- 3672
kst_km_time_he <- 408 
pca_pam_time_he <- 4.06
pca_km_time_he <- 3.92
mc2pca_km_time_he <- 51
  
  
  