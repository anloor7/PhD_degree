


apple <- read.csv('apple.csv')

s <- 4
l <- 1259

returns_apple <- Delt(apple$close, type = 'log')[s:l]
volume_apple <- Delt(apple$volume)[s:l]

mts_apple <- cbind(returns_apple, volume_apple)
mts_plot(mts_apple)

ibm <- read.csv('ibm.csv')

returns_ibm <- Delt(ibm$close, type = 'log')[s:l]
volume_ibm <- Delt(ibm$volume)[s:l]

mts_ibm <- cbind(returns_ibm, volume_ibm)


google <- read.csv('google.csv')

returns_google <- Delt(google$close, type = 'log')[s:l]
volume_google <- Delt(google$volume)[s:l]

mts_google <- cbind(returns_google, volume_google)



amazon <- read.csv('amazon.csv')

returns_amazon <- Delt(amazon$close, type = 'log')[s:l]
volume_amazon <- Delt(amazon$volume)[s:l]

mts_amazon <- cbind(returns_amazon, volume_amazon)


bax <- read.csv('bax.csv')

returns_bax <- Delt(bax$close, type = 'log')[s:l]
volume_bax <- Delt(bax$volume)[s:l]

mts_bax <- cbind(returns_bax, volume_bax)

bdx <- read.csv('bdx.csv')

returns_bdx <- Delt(bdx$close, type = 'log')[s:l]
volume_bdx <- Delt(bdx$volume)[s:l]

mts_bdx <- cbind(returns_bdx, volume_bdx)


cxo <- read.csv('cxo.csv')

returns_cxo <- Delt(cxo$close, type = 'log')[s:l]
volume_cxo <- Delt(cxo$volume)[s:l]

mts_cxo <- cbind(returns_cxo, volume_cxo)

cop <- read.csv('cop.csv')

returns_cop <- Delt(cop$close, type = 'log')[s:l]
volume_cop <- Delt(cop$volume)[s:l]

mts_cop <- cbind(returns_cop, volume_cop)


mtss <- list(mts_apple, mts_ibm, mts_google, mts_amazon, 
             mts_bax, mts_bdx, mts_cxo, mts_cop)





# Clustering 

coherence2 <- listTomatrix(lapply(mtss, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, 2)$cluster



# Añadiendo dos compañías de restaurantes

mcd <- read.csv('mcd.csv')

returns_mcd <- Delt(mcd$close, type = 'log')
volume_mcd <- Delt(mcd$volume, type = 'log')

mts_mcd <- cbind(returns_mcd, volume_mcd)


cmg <- read.csv('cmg.csv')

returns_cmg <- Delt(cmg$close, type = 'log')
volume_cmg <- Delt(cmg$volume, type = 'log')

mts_cmg <- cbind(returns_cmg, volume_cmg)


ko <- read.csv('ko.csv')

returns_ko <- Delt(ko$close, type = 'log')
volume_ko <- Delt(ko$volume, type = 'log')

mts_ko <- cbind(returns_ko, volume_ko)

pep <- read.csv('pep.csv')

returns_pep <- Delt(pep$close, type = 'log')
volume_pep <- Delt(pep$volume, type = 'log')

mts_pep <- cbind(returns_pep, volume_pep)


pfe <- read.csv('pfe.csv')

returns_pfe <- Delt(pfe$close, type = 'log')
volume_pfe <- Delt(pfe$volume, type = 'log')

mts_pfe <- cbind(returns_pfe, volume_pfe)


prgo <- read.csv('prgo.csv')

returns_prgo <- Delt(prgo$close, type = 'log')
volume_prgo <- Delt(prgo$volume, type = 'log')

mts_prgo <- cbind(returns_prgo, volume_prgo)


mtss <- list(mts_pfe, mts_prgo, mts_cxo, mts_cop)


coherence2 <- listTomatrix(lapply(mtss, quantile_coherence_re_im))
dis_matrix <- proxy::dist(coherence2, EuclideanDistance)  
clustering <- pam(dis_matrix, 2)$cluster

