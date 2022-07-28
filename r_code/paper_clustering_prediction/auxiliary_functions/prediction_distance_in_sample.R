

prediction_distance_in_sample <- function(series_test, series_matrix, vector_coefs, mean_vector) {
  
  predictions <- series_matrix %*% t(t(vector_coefs)) + t(t(mean_vector))
  return(mean(abs(series_test - predictions)))
  
}
