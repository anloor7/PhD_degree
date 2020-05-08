

# This is a function to plot a MTS via ggplot2

# Input:
# series: matrix or dataframe containing MTS object (rows being time observations and columns being variables)

mts_plot <- function(series){
  df_series <- as.matrix(series)
  
  df <- NULL
  temp_df <- NULL
  n <- nrow(df_series)
  c <- ncol(df_series)
  for(i in 1 : c){
    temp_df <- data.frame(x = 1:n, y = df_series[,i], col = rep(i, n))
    df <- rbind(df,temp_df)} 
  
  ggplot(df, aes(x = x, y = y, group = col, colour = factor(col))) + geom_line(size = 0.8) + xlab('Time') +
    ylab('Series') + theme(legend.position = "none")
}

