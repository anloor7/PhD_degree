

# Visualizing a MTS with ggplot2. Example

phi_c1 <- matrix(c(0.6, -0.4, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3, 0.7, -0.5, 0, 0.5, 0.3, -0.4, 0.3, 0.5, -0.2, 0.8, 0, 
                   -0.5, 0.4, -0.3, 0, 0, 0.4, 0.2, 0.5, -0.7, 0.2, 0, 0.5, 0.5, 0.5, 0, -0.6, 0.2), nrow = 6)
sigma_c1 <- matrix(c(1, 0.25, 0.1, 0.25, 0.10, 0.25, 0.25, 1, 0.25, 0.10, 0.25, 0.10, 0.10, 0.25, 1, 0.10, 0, 0.25, 0.25,
                     0.10, 0.10, 1, 0.10, 0.25, 0.10, 0.25, 0, 0.10, 1, 0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 1), nrow = 6)

series <- VARMAsim(100, arlags = 1, phi = phi_c1, sigma = sigma_c1)$series # Object MTS
df_series <- as.matrix(series) # Converting the object to a data frame

df <- NULL
temp_df <- NULL
n <- nrow(df_series)
c <- ncol(df_series)
for(i in 1 : c){
  temp_df <- data.frame(x = 1:n, y = df_series[,i], col = rep(i , n))
  df <- rbind(df,temp_df)} 

ggplot(df, aes(x = x, y = y, group = col, colour = factor(col))) + geom_line() + xlab('Time') +
  ylab('Series') + theme(legend.position = "none")



