

# Function to convert a list into a matrix



# This is a function to transform a list of vectors of the same length into a matrix

listTomatrix <- function(l){
  n <- length(l)
  s <- length(l[[1]])
  m <- matrix(nrow = n, ncol = s)
  for (i in 1 : n) {
    m[i,] <- l[[i]]
  }
  m
}


# Function for computing the marginal probabilities

marginal_probabilities <- function(series, categories) {
  
  series <- factor(series, levels = levels(categories))
  series_length <- length(series) # Series length
  n_cat <- length(categories) # Number of categories in the dataset
  
  
  # Computing the marginal probabilities in the series
  
  marginal_probabilities <- numeric()
  
  for (i in 1 : n_cat) {
    
    count_i <- sum(series == categories[[i]])
    marginal_probabilities[i] <- count_i/series_length
    
  }
  
  return(marginal_probabilities)
  
}


# Functions for computing the Cramer's v-based features

# Auxiliary function to compute p_i_j_k

p_i_j_k_function <- function(series, i_cat, j_cat, k) {
  
  series_length <- length(series)
  a <- series[(k + 1) : series_length]
  b <- series[1 : (series_length - k)]
  
  
  number <- series_length - k
  count <- numeric(number)
  
  for (i in 1 : number) {
    
    if (a[i] == i_cat & b[i] == j_cat) {
      
      count[i] <- 1
      
    } else {
      
      count[i] <- 0
      
    }
    
  }
  
  
  return(sum(count)/(series_length-k))
  
}


# Main function

features_cramer_v_function <- function(series, categories, max_lag = 1) {
  
  series <- factor(series, levels = levels(categories))
  series_length <- length(series) # Series length
  n_cat <- length(categories) # Number of categories in the dataset
  
  
  # Computing the marginal probabilities in the series
  
  marginal_probabilities <- numeric()
  
  for (i in 1 : n_cat) {
    
    count_i <- sum(series == categories[[i]])
    marginal_probabilities[i] <- count_i/series_length
    
  }
  
  # Computing the term s2p
  
  s2p <- sum(marginal_probabilities^2)
  
  # Computing the factors inside the square root of Kruskal's tau
  
  features <- numeric()
  
  count <- 1
  for (i in 1 : n_cat) {
    
    for (j in 1 : n_cat) {
      
      for (k in 1 : max_lag) {
        
        p_i_j_k <- p_i_j_k_function(series = series, i_cat = categories[i], j_cat = categories[j], k = k)
        numerator <- p_i_j_k - marginal_probabilities[i]*marginal_probabilities[j]
        denominator <- marginal_probabilities[i]*marginal_probabilities[j]
        features[count] <- (numerator^2/denominator)
        count <- count + 1
        
      }
      
    }
    
  }
  
  features[is.na(features)] <- 0
  return(features)
  
  
  
}


# Functions to compute the Cohen's k-based features

# Auxiliary function to compute p_i_j_k

p_i_j_k_function <- function(series, i_cat, j_cat, k) {
  
  series_length <- length(series)
  a <- series[(k + 1) : series_length]
  b <- series[1 : (series_length - k)]
  
  
  number <- series_length - k
  count <- numeric(number)
  
  for (i in 1 : number) {
    
    if (a[i] == i_cat & b[i] == j_cat) {
      
      count[i] <- 1
      
    } else {
      
      count[i] <- 0
      
    }
    
  }
  
  
  return(sum(count)/(series_length-k))
  
}


# Main function

features_cohen_function <- function(series, categories, max_lag = 1) {
  
  series <- factor(series, levels = levels(categories))
  series_length <- length(series) # Series length
  n_cat <- length(categories) # Number of categories in the dataset
  
  
  # Computing the marginal probabilities in the series
  
  marginal_probabilities <- numeric()
  
  for (i in 1 : n_cat) {
    
    count_i <- sum(series == categories[[i]])
    marginal_probabilities[i] <- count_i/series_length
    
  }
  
  # Computing the term s2p
  
  s2p <- sum(marginal_probabilities^2)
  
  # Computing the factors inside the square root of Kruskal's tau
  
  features <- numeric()
  
  count <- 1
  for (i in 1 : n_cat) {
    
    for (k in 1 : max_lag) {
      
      p_i_j_k <- p_i_j_k_function(series = series, i_cat = categories[i], j_cat = categories[i], k = k)
      numerator <- p_i_j_k - marginal_probabilities[i]^2
      denominator <- 1 - s2p
      features[count] <- numerator/denominator
      count <- count + 1
      
    }
    
  }
  
  features[is.na(features)] <- 0
  return(features)
  
  
  
}


# Function to compute the binarization-based features

features_categorical_series <- function(series, l = 1, levels_dataset) {
  
  l_series <- length(series) # Length of each series
  n_cat <- length(levels_dataset) # Number of categories
  
  
  # Feature extraction for each UTS
  
  # Creating a matrix containing the series of indicators 
  
  
  matrix_series <- matrix(0, nrow = n_cat, ncol = l_series)
  
  for (j in 1 : n_cat) {
    
    matrix_series[j,] <- as.numeric(series == levels_dataset[j])
    
  }
  
  series_indicators <- matrix_series
  
  
  
  
  # Creating the vector of features
  
  features_autocorrelations <- matrix(0, nrow = n_cat, ncol = l)
  features_cross_correlations <- matrix(0, nrow = n_cat*(n_cat-1)/2, ncol = 2*l + 1)
  
  
  for (i in 1 : n_cat) {
    
    features_autocorrelations[i,] <- TSA::acf(series_indicators[i,], lag.max = l, plot = F)$acf[,,1]
    
  }
  
  features_autocorrelations[is.na(features_autocorrelations)] <- 0
  
  
  count <- 1
  for (i in 1 : n_cat) {
    
    for(j in i : n_cat) {
      
      if (i != j) { 
        
        features_cross_correlations[count,] <- stats::ccf(series_indicators[i,], 
                                                          series_indicators[j,], 
                                                          lag.max = l, plot = F)$acf[,,1]
        
        count <- count + 1
        
      }
      
    }
    
    
  }
  
  features_cross_correlations[is.na(features_cross_correlations)] <- 0
  
  
  
  
  # Constructing the whole vector of features
  
  whole_vector_features <- c(c(features_autocorrelations), c(features_cross_correlations))
  return_list <- list(whole_features = whole_vector_features, features_cross_correlations = features_cross_correlations)
  
  return(return_list)
  
  
  
  
  
  
  
  
  
}


# Function to simulate NDARMA processes



# This is a function for simulating an NDARMA(p, q) process
# Input parameters:
# marg_prob: marginal probabilities concerning the series of noises
# probs_multin: a vector of length p + q + 1 containing the probabilities concerning the multinomial distribution
# p: value of p
# q: value of q
# series_length: length of the series that we want to simulate

ndarma_sim_function <- function(marg_prob, probs_multin, p, q, series_length) {
  
  n_cat <- length(marg_prob) # Number of categories
  
  # First step (i): constructing the p values X_t-1, ..., X_t-p
  
  p_first <- sample(0 : (n_cat-1), size = p, replace = T, prob = marg_prob)
  
  # First step (ii): constructing the q + 1 values epsilon_t, epsilon_t-1, ..., epsilon_t-q
  
  q_1_first <- sample(0 : (n_cat-1), size = q + 1, replace = T, prob = marg_prob)
  
  # Joining the two vectors
  
  first_parameters <- c(p_first, q_1_first)
  
  
  # Creating the series
  
  series <- numeric()
  
  
  for (i in 1 : series_length) {
    
    epsilon_t <- sample(0 : (n_cat-1), size = 1, replace = T, prob = marg_prob)
    
    if (q != 0) {
      
      first_parameters[(p+1) : (p+q+1)] <- c(epsilon_t, first_parameters[(p+1) : (p+q)])
      
    } else {
      
      first_parameters[p + 1] <- epsilon_t
      
    }
    
    mult_vector <- as.numeric(rmultinom(1, 1, prob = probs_multin))
    series[i] <- sum(mult_vector*first_parameters)
    
    if (p == 1) {
      
      first_parameters[1 : p] <- c(series[i])
      
    } else if (p > 1) {
      
      first_parameters[1 : p] <- c(series[i], first_parameters[1 : (p-1)])
      
    }
    
    
    
  }
  
  return(series)
  
  
}



