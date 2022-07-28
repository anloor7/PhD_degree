

prediction_matrix <- function(series, max_lag){
  
  series_length <- length(series)
  
  # Creating a dataset which contains, in each row, the vectors subjected
  # to be predicted
  
  dataset_predictions <- matrix(0, nrow = series_length - max_lag, ncol = max_lag)
  starting_vector <- 1 : max_lag
  
  for (i in (1 : (series_length - max_lag))) {
  
    dataset_predictions[i,] <- series[starting_vector + (i - 1)]
    
  }
  
  return(dataset_predictions)
  
}
