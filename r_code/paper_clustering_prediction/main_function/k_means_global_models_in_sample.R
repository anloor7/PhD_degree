

# This is a function to compute the k-means approach through global models and
# the in-sample error

library(caret)

k_means_global_models_in_sample <- function(series, K, max_lag = 1, niter = 150, model_caret = 'lm',
                                  init = NULL) {
  
  n_series <- length(series)
  series_lengths <- lapply(series, length)
  
  # Computing the matrices associated to each series (later these matrices will be
  # used to compute the predictions)
  
  prediction_matrices <- lapply(series, prediction_matrix, max_lag = max_lag)
  
  # Random inicialization of centroids
  
  if (is.null(init)) {
    
    sample_numbers <- sample(1 : K)
    sample_numbers <- c(sample_numbers, 
                        sample(K, n_series - K, replace = T)) # Random assignment of series to clusters
    initial_models <- list()
    
    for (i in 1 : K) {
      
      indexes_k <- which(sample_numbers == i)
      initial_models[[i]] <- global_model_function(series[indexes_k], max_lag = max_lag)
      
    }
    
  } else {
    
    initial_models <- list()
    
    for (i in 1 : K) {
      
      indexes_k <- which(init == i)
      initial_models[[i]] <- global_model_function(series[indexes_k], max_lag = max_lag)
      
    }
    
  }
  
  
  
  # Repeat
  
  indexes_assignment <- 0
  old_assignment <- 1
  models <- initial_models
  iter <- 0
  objective_function <- numeric()
  clustering_solution_list <- list()
  centroids_list <- list()
  mae_iter <- numeric()
  mae_vector_list <- list()
  
  while (iter < niter & sum(indexes_assignment == old_assignment) != length(indexes_assignment) & if(iter > 7) {objective_function[iter - 6] > min(c(objective_function[iter - 5], 
                                                                                                                                                     objective_function[iter - 4],
                                                                                                                                                     objective_function[iter - 3],
                                                                                                                                                     objective_function[iter - 2],
                                                                                                                                                     objective_function[iter - 1]))} else {T}) {
    
    
    # Assigning each series to its closest centroid
    
    old_assignment = indexes_assignment
    indexes_assignment <- numeric(n_series)
    store_distances <- numeric(n_series)
    
    
    for (j1 in 1 : n_series){
      
      distances <- numeric()
      
      for(j2 in 1 : K){
        
        vector_coefficients <- summary(models[[j2]])$coef[,1][-1]
        if (length(vector_coefficients) != max_lag) {
          
          vector_coefficients <- c(vector_coefficients, 
                                   rep(0, max_lag - length(summary(models[[j2]])$coef[,1][-1])))
        }
          
        distances[j2] <- prediction_distance_in_sample(series[[j1]][-(1 : max_lag)],
                                                      prediction_matrices[[j1]],
                                                      vector_coefficients,
                                                      rep(summary(models[[j2]])$coef[,1][1], series_lengths[[j1]] - max_lag))
        
      }
      
      indexes_assignment[j1] <- which.min(distances)
      store_distances[j1] <- min(distances)
      
    }
    
    u_cluster_vector <- unique(indexes_assignment)
    
    if (length(u_cluster_vector) != K) {
      n_clusters <- length(unique(indexes_assignment))
      
      for (i in 1 : n_series) {
        
        indexes_assignment[i] <- which(indexes_assignment[i] == 
                                         u_cluster_vector)
      }
      
      K <- n_clusters
      
    }
    
    
    # Recomputation of centroids (models)
    
    
    for (p in 1 : K) {
      
      indexes_k <- which(indexes_assignment == p)
      models[[p]] <- global_model_function(series[indexes_k], max_lag = max_lag,
                                           model_caret = model_caret)
      
    }
    
    
    iter <- iter + 1
    
    objective_function[iter] <- mean(store_distances)
    clustering_solution_list[[iter]] <- indexes_assignment
    centroids_list[[iter]] <- models
    mae_iter[iter] <- mean(store_distances)
    mae_vector_list[[iter]] <- store_distances
    
  }
  
  clustering_solution_list <- clustering_solution_list[-iter]
  clustering_solution_list <- c(list(sample_numbers), clustering_solution_list)
  
  
  # Computing the value of the objective function
  
  
  # distances_min <- numeric()
  # for (j1 in 1 : n_series){
  
  #  index <- indexes_assignment[j_1]
  #  distances_min[j1] <- prediction_distance(training_series[[j1]], test_series[[j1]],
  #                                         models[[index]], max_lag = max_lag)
  
  
  # }
  
  if (iter > 7) {
  
  if (objective_function[iter - 6] <= min(c(objective_function[iter - 5], 
                                                       objective_function[iter - 4],
                                                       objective_function[iter - 3],
                                                       objective_function[iter - 2],
                                                       objective_function[iter - 1]))) {
    
    
    
    list_return <- list(clustering = clustering_solution_list[[iter - 6]],
                        centroids = centroids_list[[iter - 6]],
                        niter = iter,
                        mae = mae_iter[iter - 6],
                        mae_vector = mae_vector_list[[iter - 6]],
                        objective_function = objective_function,
                        clustering_solution_list = clustering_solution_list,
                        centroids_list = centroids_list)
    
  } else { 
    
    list_return <- list(clustering = indexes_assignment,
                        centroids = models,
                        niter = iter,
                        mae = mean(store_distances),
                        mae_vector = store_distances,
                        objective_function = objective_function, 
                        clustering_solution_list = clustering_solution_list,
                        centroids_list = centroids_list)
    
  }
    
  } else {
    
    list_return <- list(clustering = indexes_assignment,
                        centroids = models,
                        niter = iter,
                        mae = mean(store_distances),
                        mae_vector = store_distances,
                        objective_function = objective_function,
                        clustering_solution_list = clustering_solution_list,
                        centroids_list = centroids_list)
    
  }
  
  
  return(list_return)
  
  
  
  
}

