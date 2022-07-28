

global_model_function <- function(series, max_lag = 1, model_caret = 'lm',...) {
  
  # Creating the lag-embedded matrix
  
  l <- length(series)
  
  list_matrices <- list()
  
  for (i in 1 : l) {
    series_length <- length(series[[i]])
    list_matrices[[i]] <- matrix(0, nrow = series_length - max_lag, ncol = max_lag + 1)
    
    for (j in 1 : max_lag) {
      
      list_matrices[[i]][,j] <- series[[i]][j : (j + series_length -  max_lag - 1)]
      
    }
    
    list_matrices[[i]][, (max_lag + 1)] <- series[[i]][-c(1 : max_lag)]
    
  }
  
  matrix_total <- do.call('rbind', list_matrices)
  
  
  # Constructing a Machine Learning model over the lag-embedded matrix
  # via the caret library
  
  last_col <- ncol(matrix_total)
  train_explanatory <- data.frame(matrix_total[, 1 : (last_col - 1)])
  train_response <- matrix_total[, last_col]
  
  trControl <- trainControl(method = "none")
  if (ncol(train_explanatory) == 1) {colnames(train_explanatory) <- 'X1'}
  model <- train(train_explanatory, train_response, method = model_caret,
                 trControl = trControl,...)
  model
  
}