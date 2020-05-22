

# Simulating from scratch 

# Each variable of the series has a time-varying conditional variance and can be modeled as a univariate GARCH
# process 

vectoruspec <- list(omega = 0.5, alpha1 = 0.5, beta1 = 0.5, gamma1 = 0.5, ma1 = -0.4) # Parameters of eGARCH model 
uspec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm", mean.model = list(armaOrder=c(0, 1), include.mean = F), 
                    fixed.pars = vectoruspec)



model <- multispec(replicate(4,uspec)) # Replication of eGarch model for the four variables 

vectormvspec <- list(dcca1 = 0.4, dccb1 = 0.4) # Parameters of DCC model 
mvspec <- dccspec(model, dccOrder = c(1, 1), model = "DCC", distribution = # Defining DCC model 
                    "mvnorm", fixed.pars = vectormvspec)


sim <- dccsim(mvspec, n.sim = 300, preQ = diag(4), Qbar = diag(4), Nbar = diag(4))
simulations <- slot(sim, 'msim')[[4]][[1]] # Simulated series 

