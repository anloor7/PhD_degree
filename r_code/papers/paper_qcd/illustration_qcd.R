
 
# Lets go to ilustrate QCD in three different scenarios regarding bivariate series


setwd('/Users/angellopezoriona/Library/Mobile Documents/com~apple~CloudDocs/academic_life/PhD/papers/papers_2020/paper_qcs/stored_results_copia_seguridad')

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}




files_rdata <- list.files(pattern = "*.RData")
l_rdata <- length(files_rdata)

for (i in 1:l_rdata) {
  
  load(files_rdata[i])
  
}


# White noise processes

set.seed(1234)
x1 <- rnorm(20000)
x2 <- rnorm(20000)
P1 <- cbind(x1, x2)
# q_quantities_plot_mts(P1, part = 'Re', q1 = 0.1, q2 = 0.9)




# VARMA process 

l <- 20000
phi <- matrix(c(0.8, -0.4, 0.7, 0.6), ncol = 2)
sigma <- diag(2)
epsilon <- rmvnorm(l,  sigma = sigma) # Normal
P2 <- varma(l, k = 2, VAR = phi, innov = epsilon)
# q_quantities_plot_mts(P2, part = 'Re', q1 = 0.1, q2 = 0.9)




# QVAR processes

l <- 20000
add <- 150
l <- l + add
Xt <- list()
Xt[[1]] <- c(0, 0)

for (i in 2 : (l + 1)) {
  
  theta <- 0*diag(2)
  u12 <- runif(1)
  u21 <- runif(1)
  theta[1, 2] <- 1.5*(u12 - 0.5)
  theta[2, 1] <- 1.5*(u21 - 0.5)
  
  theta0 <- numeric()
  theta0[1] <- qnorm(u12, 0)
  theta0[2] <- qnorm(u21, 0)
  
  Xt[[i]] <- theta %*% Xt[[i - 1]] + theta0
  
}

cluster1 <- listTomatrix(Xt)
P3 <- cluster1[((add+1):l),]
# q_quantities_plot_mts(P3, part = 'Re', q1 = 0.1, q2 = 0.9)



# Plots




# Plot 1

title_size = 15.5
text_size = 13.5
legend_size = 15
x_size = text_size + 1
q1 <- 0.5
q2 <- 0.1

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df11 <- data.frame(x = freq, real = real_part)
plot1 <- ggplot(df11, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Quantile cospectrum'),
       subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.1$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                                       size = title_size),
                                                                             axis.text = element_text(size = text_size),
                                                                             axis.title = element_text(size = x_size),
                                  plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df11, file = 'df11.RData')

c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df12 <- data.frame(x = freq, real = real_part)
plot1 <- plot1 + geom_line(data = df12, aes(x = x/(2*pi), y = real), col = 'blue')
# save(df12, file = 'df12.RData')

c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df13 <- data.frame(x = freq, real = real_part)
plot1 <- plot1 + geom_line(data = df13, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df13, file = 'df13.RData')

S1 <- rep('P1', 10001)
S2 <- rep('P2', 10001)
S3 <- rep('P3', 10001)
S123 <- c(S1, S2, S3)
  
df_total <- rbind(df11, df12, df13)
df_total$Process <- S123

plot1_legend <- ggplot(df_total, aes(x = x, y = real, color = Process)) + geom_line() +
  scale_colour_manual(values=c("red","blue","green")) +
  scale_fill_discrete(name = "Process") +
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  ggtitle(TeX('Quantile cospectrum ($\\tau_1=0.5$, $\\tau_2=0.1$)')) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "top",
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size))

shared_legend <- extract_legend(plot1_legend)

# Plot 2 


q1 <- 0.5
q2 <- 0.1

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df21 <- data.frame(x = freq, real = real_part)
plot2 <- ggplot(df21, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Quantile quadrature spectrum'),
       subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.1$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                           size = title_size),
                                                                 axis.text = element_text(size = text_size),
                                                                 axis.title = element_text(size = x_size),
                                                                 plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df21, file = 'df21.RData')  


c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df22 <- data.frame(x = freq, real = real_part)
plot2 <- plot2 + geom_line(data = df22, aes(x = x/(2*pi), y = real), col = 'blue') 
# save(df22, file = 'df22.RData')


c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df23 <- data.frame(x = freq, real = real_part)
plot2 <- plot2 + geom_line(data = df23, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df23, file = 'df23.RData')

plot12 <- grid.arrange(arrangeGrob(plot1, plot2, ncol = 2))

plot12 <- grid.arrange(
             shared_legend, nrow = 2,
             arrangeGrob(plot1, plot2, nrow = 1, ncol = 2), ncol = 1,
             heights = c(1, 10))



# Plot 3

q1 <- 0.5
q2 <- 0.5

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df31 <- data.frame(x = freq, real = real_part)
plot3 <- ggplot(df31, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') + labs(title = TeX('Quantile cospectrum'),
                                               subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.5$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                                                                   size = title_size),
                                                                                                         axis.text = element_text(size = text_size),
                                                                                                         axis.title = element_text(size = x_size),
                                                                                                         plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df31, file = 'df31.RData')

c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df32 <- data.frame(x = freq, real = real_part)
plot3 <- plot3 + geom_line(data = df32, aes(x = x/(2*pi), y = real), col = 'blue')
# save(df32, file = 'df32.RData')

c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df33 <- data.frame(x = freq, real = real_part)
plot3 <- plot3 + geom_line(data = df33, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df33, file = 'df33.RData')



# Plot 4 


q1 <- 0.5
q2 <- 0.5

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df41 <- data.frame(x = freq, real = real_part)
plot4 <- ggplot(df41, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Quantile quadrature spectrum'),
       subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.5$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                           size = title_size),
                                                                 axis.text = element_text(size = text_size),
                                                                 axis.title = element_text(size = x_size),
                                                                 plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df41, file = 'df41.RData')  
  
  
c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df42 <- data.frame(x = freq, real = real_part)
plot4 <- plot4 + geom_line(data = df42, aes(x = x/(2*pi), y = real), col = 'blue') 
# save(df42, file = 'df42.RData')


c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df43 <- data.frame(x = freq, real = real_part)
plot4 <- plot4 + geom_line(data = df43, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df43, file = 'df43.RData')

plot34 <- grid.arrange(plot3, plot4, nrow = 1)





# Plot 5


q1 <- 0.5
q2 <- 0.9

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df51 <- data.frame(x = freq, real = real_part)
plot5 <- ggplot(df51, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Quantile cospectrum'),
       subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.9$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                           size = title_size),
                                                                 axis.text = element_text(size = text_size),
                                                                 axis.title = element_text(size = x_size),
                                                                 plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df51, file = 'df51.RData')

c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df52 <- data.frame(x = freq, real = real_part)
plot5 <- plot5 + geom_line(data = df52, aes(x = x/(2*pi), y = real), col = 'blue')
# save(df52, file = 'df52.RData')

c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- Re(qSPGv[,1,1,2, 1, 1])
df53 <- data.frame(x = freq, real = real_part)
plot5 <- plot5 + geom_line(data = df53, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df53, file = 'df53.RData')



# Plot 6


q1 <- 0.5
q2 <- 0.9

c <- ncol(P1)
qSPG <- smoothedPG(P1, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df61 <- data.frame(x = freq, real = real_part)
plot6 <- ggplot(df61, aes(x = x/(2*pi), y = real)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Quantile quadrature spectrum'),
       subtitle = TeX(('$\\tau_1=0.5$, $\\tau_2=0.9$'))) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                           size = title_size),
                                                                 axis.text = element_text(size = text_size),
                                                                 axis.title = element_text(size = x_size),
                                                                 plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 

# save(df61, file = 'df61.RData')


c <- ncol(P2)
qSPG <- smoothedPG(P2, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df62 <- data.frame(x = freq, real = real_part)
plot6 <- plot6 + geom_line(data = df62, aes(x = x/(2*pi), y = real), col = 'blue') 
  
# save(df62, file = 'df62.RData')


c <- ncol(P3)
qSPG <- smoothedPG(P3, levels.1 = q1, levels.2 = q2, type = 'qr')
freq <- getFrequencies(qSPG) # Fourier frequencies
qSPGv <- getValues(qSPG, frequencies = freq)
l <- length(freq)
real_part <- -Im(qSPGv[,1,1,2, 1, 1])
df63 <- data.frame(x = freq, real = real_part)
plot6 <- plot6 + geom_line(data = df63, aes(x = x/(2*pi), y = real), col = 'green') +
  theme(axis.text = element_text(size = text_size))
# save(df63, file = 'df63.RData')

plot56 <- grid.arrange(plot5, plot6, nrow = 1)


# Plot 7

c <- ncol(P1)
cs <- crossSpectrum(P1, spans = 1000)
freq <- cs$freq
re_per <- Re(cs$Pxy)
im_per <- -Im(cs$Pxy)
l <- length(freq)
text_size = 12
df71 <- data.frame(x = freq, per = re_per)
df71_i <- data.frame(x = freq, per = im_per)
plot71 <- ggplot(df71, aes(x = x, y = per)) + geom_line(col = 'red') + 
xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional cospectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                                                                size = title_size),
                                                                                                      axis.text = element_text(size = text_size),
                                                                                                      axis.title = element_text(size = x_size),
                                                                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-1, 1, by = 0.1))
# save(df71, file = 'df71.RData')

plot71_i <- ggplot(df71, aes(x = x, y = im_per)) + geom_line(col = 'red') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional quadrature spectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                size = title_size),
                                                      axis.text = element_text(size = text_size),
                                                      axis.title = element_text(size = x_size),
                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-1, 1, by = 0.1))
# save(df71, file = 'df71.RData')

trad_plot_1 <- grid.arrange(plot71, plot71_i, nrow = 1)


c <- ncol(P2)
cs <- crossSpectrum(P2, spans = 1000)
freq <- cs$freq
re_per <- Re(cs$Pxy)
im_per <- -Im(cs$Pxy)
l <- length(freq)
df72 <- data.frame(x = freq, per = re_per)
df72_i <- data.frame(x = freq, per = im_per)
plot72 <- ggplot(df72, aes(x = x, y = per)) + geom_line(col = 'blue') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional cospectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                size = title_size),
                                                      axis.text = element_text(size = text_size),
                                                      axis.title = element_text(size = x_size),
                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) +
  scale_y_continuous(limits = c(-0.1, 4.8), breaks = seq(0, 4, by = 1))
# save(df71, file = 'df71.RData')
# save(df72, file = 'df72.RData')
plot72_i <- ggplot(df72_i, aes(x = x, y = per)) + geom_line(col = 'blue') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional quadrature spectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                size = title_size),
                                                      axis.text = element_text(size = text_size),
                                                      axis.title = element_text(size = x_size),
                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) 
# save(df71, file = 'df71.RData')

trad_plot_2 <- grid.arrange(plot72, plot72_i, nrow = 1)

c <- ncol(P3)
cs <- crossSpectrum(P3, spans = 1000)
freq <- cs$freq
re_per <- Re(cs$Pxy)
im_per <- -Im(cs$Pxy)
l <- length(freq)
df73 <- data.frame(x = freq, per = re_per)
df73_i <- data.frame(x = freq, per = im_per)
plot73 <- ggplot(df73, aes(x = x, y = per)) + geom_line(col = 'green') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional cospectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                size = title_size),
                                                      axis.text = element_text(size = text_size),
                                                      axis.title = element_text(size = x_size),
                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-1, 1, by = 0.1))
# save(df71, file = 'df71.RData')
# save(df72, file = 'df72.RData')
plot73_i <- ggplot(df73_i, aes(x = x, y = per)) + geom_line(col = 'green') + 
  xlab(TeX('$\\omega/2\\pi$')) + ylab('') +
  labs(title = TeX('Traditional quadrature spectrum')) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                size = title_size),
                                                      axis.text = element_text(size = text_size),
                                                      axis.title = element_text(size = x_size),
                                                      plot.subtitle =  element_text(hjust = 0.5, size = title_size)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-1, 1, by = 0.1))
# save(df71, file = 'df71.RData')

trad_plot_3 <- grid.arrange(plot73, plot73_i, nrow = 1)

plot789 <- grid.arrange(shared_legend, nrow = 2,
                        arrangeGrob(trad_plot_1, trad_plot_2, trad_plot_3, nrow = 3, ncol = 1), ncol = 1,
                        top = textGrob("", 
                                     vjust = 1.8, just = 'center', 
                                     gp = gpar(col = "black", fontsize = 14)),
                        heights = c(0.5, 10))


# plot12 <- grid.arrange(arrangeGrob(plot1, plot2, ncol = 2),
#                        shared_legend, nrow = 2, heights = c(10, 1))

final_plot <- grid.arrange(plot12, plot34, plot56, nrow = 3)






