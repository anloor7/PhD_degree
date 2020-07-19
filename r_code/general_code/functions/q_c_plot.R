

# This is a function to return plots of quantile coherency (real part) and quantile coherency (imaginary part) given 
# two UTS components of a given MTS (both components can be the same)

# Input parameters:
# X,Y: two UTS making the MTS of interest 
# q1, q2: levels of probability that we want to take into consideration

q_quantities_plot <- function(X, Y, q1, q2){
  
  
  Z <- cbind(X, Y)
  l <- nrow(Z)
  qSPG <- smoothedPG(Z, levels.1 = q1, levels.2 = q2, type = 'qr')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  
  m <- qSPGv[, 2, 1, 1, 1, 1]
  real_part <- Re(m)
  imaginary_part <- -Im(m)
  df <- data.frame(x = freq, real = real_part, imaginary = imaginary_part)
  
  plot_real <- ggplot(df, aes(x = freq/(2*pi), y = real)) + geom_line(col = 'red') + 
    xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
    ggtitle('Quantile cospectrum')
  plot_imaginary <- ggplot(df, aes(x = freq/(2*pi), y = imaginary)) + geom_line(col = 'blue') + 
    xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
    ggtitle('Quantile quadrature spectrum')
  grid.arrange(plot_real, plot_imaginary)
  
}


q_c_plot <- function(X, Y, q1 = 0.5, q2 = 0.5){
  
  
  Z <- cbind(X, Y)
  l <- nrow(Z)
  qSPG <- smoothedPG(Z, levels.1 = q1, levels.2 = q2, type = 'qr')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getCoherency(qSPG, frequencies = freq)
  
  m <- qSPGv[, 2, 1, 1, 1, 1]
  real_part <- Re(m)
  imaginary_part <- Im(m)
  df <- data.frame(x = freq, real = real_part, imaginary = imaginary_part)
  
  plot_real <- ggplot(df, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
    xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
    ggtitle('Quantile coherency (real part)')
  plot_imaginary <- ggplot(df, aes(x = freq/(2*pi), y = imaginary)) + geom_line(col = 'blue') + 
    xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
    ggtitle('Quantile coherency (imaginary part)')
  grid.arrange(plot_real, plot_imaginary)
  
}


# This is a function to returns plots of quantile cospectrum and quantile quadrature spectrum, either real 
# or imaginary part for any given pair of UTS within a MTS

# Input parameters:
# X: a MTS
# part: 'Re' for real part and 'Im' for imaginary part
# q1, q2: levels of probability that we want to take into consideration


q_quantities_plot_mts <- function(X, part = 'Re', q1 = 0.5, q2 = 0.5){
  
  
  c <- ncol(X)
  qSPG <- smoothedPG(X, levels.1 = q1, levels.2 = q2, type = 'qr')
  freq <- getFrequencies(qSPG) # Fourier frequencies
  qSPGv <- getValues(qSPG, frequencies = freq)
  l <- length(freq)
  qSPGc <- getCoherency(qSPG, frequencies = freq)
  
  
  if (part == 'Re'){
    
    k <- 1
    re_plots <- list()
    
    for (i in (1 : c)) {
      
      for (j in (1 : c)){
        
        
        real_part <- Re(qSPGv[,i,1,j, 1, 1])
        df <- data.frame(x = freq, real = real_part)
        re_plots[[k]] <-  ggplot(df, aes(x = freq/(2*pi), y = real)) + geom_line(col = 'red') + 
          xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
          ggtitle('Quantile coherency (real part)')
        k <- k + 1
        
      }
    }
  
    do.call("grid.arrange", c(re_plots, ncol = c))
    
  } 
  
  
  if (part == 'Im') {
    
    k <- 1
    im_plots <- list()
    
    for (i in (1 : c)) {
      
      for (j in (1 : c)){
        
        im_part <- Im(qSPGv[,i,1,j, 1, 1])
        df <- data.frame(x = freq, im = im_part)
        im_plots[[k]] <-  ggplot(df, aes(x = freq/(2*pi), y = im)) + geom_line(col = 'blue') + 
          xlab(TeX('$\\frac{\\omega}{2\\pi}$')) + ylab('') +
          ggtitle('Quantile coherency (imaginary  part)')
        k <- k + 1
        
      }
    }
    
    do.call("grid.arrange", c(im_plots, ncol = c))
    
    
  }
  
  
  
}
