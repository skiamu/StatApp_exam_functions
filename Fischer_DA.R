# funzione per la Fischer Discriminant Analysis
# 
# INPUT:
#       X = dataframe in cui ogni riga è un'osservazione delle caratteristiche
#           n righe p colonne
#       group = factor vector con le label che indicano la popolazione di 
#               appartenenza per ogni osservazioni in X
#       n_dp = number of discriminant components   
#       print.resul = 
#       print.plot =
#       X.new = dataframe con le caratteristiche a cui si vuole attribuire una
#               classe col classificatore do fischer       
# OUTPUT:
# 
# 
library(dplyr)

Fisher_DA <- function(X, 
                      group, jittering = T,
                      n_dp = 2,
                      print.result = T,
                      print.plot = T,
                      X.new){
   set.seed(280787)
   sd_jitt = 0.025
   # smuovere un po i dati se le misure sono state fatte con arronondamneto
   if(jittering){
      X <- X + cbind(rnorm(n, sd = sd_jitt))
   }
   # numerosità campionaria
   n <- dim(X)[1]
   # dimensione delle caratteristiche
   p <- dim(X)[2]
   # numerosità campionaria per ogni gruppo
   ng <- table(group)
   # numero di gruppi
   g <- length(ng)
   # media dentro ogni gruppo (medie vettore colonna)
   Mg <- matrix(ncol = g,nrow = p)
   for(i in 1:g){
      Mg[,i] <- colMeans( X[which(group == levels(group)[i]),] )
   }
   # grandmean
   M <- colMeans(X)
   # inizializzo S a 0
   Sp <- matrix(0, ncol = p, nrow = p)
   # calcolo una stima di sigma, Spool
   for(i in 1:g){
      # covarianza dentro gruppo i
      Sg <- (ng[i] - 1) / (n-g) * cov ( X[which(group == levels(group)[i]),] ) 
      Sp <- Sp +  Sg
   }
   # ATTENZIONE: verificare se si deve dividere per n o n-g
   B <- matrix(0, ncol = p, nrow = p)
   for(i in 1:g){
      B <- B + ng[i]/n * (Mg[,i] - M)%*% t((Mg[,i] - M))
   }
   
   val.Sp <- eigen(Sp)$val
   vec.Sp <- eigen(Sp)$vec
   invSp.2 <- matrix(0, nrow = p, ncol = p)
   # calcolo la matrice sigma^-1/2
   for(i in 1:p){
      invSp.2 <- invSp.2 + 1 / sqrt(val.Sp[i]) * vec.Sp[,i] %*% t(vec.Sp[,i])
   }
   # calcolo decomposizione spettrale per la matrice sigma^-1/2 * B * sigma^-1/2
   spec.dec <- eigen(invSp.2 %*% B %*% invSp.2)
   # matrice le cui colonne sone le discriminant direction a1,.....ap
   A <- invSp.2 %*% spec.dec$vectors
   # calcolo gli scores nelle prime n_dp discriminant components, ossia proietto
   # le carattaristiche nello spazio delle prime n_dp discriminant component
   scores <- as.matrix(X) %*% A[,1:n_dp]
   
   ########## CLASSIFICATION ###########
   # coordinate canoniche dei vettori delle medie nei gruppi nello spazio delle
   # prime n_dp discriminant component
   cc.M <- t(Mg) %*% A[,1:n_dp]
   
   # calcola la classe con il classificatore di fisher dei dati nel training set
   # per poi calcolare una stima dell'APER. x assegnato alla popolazione k se la
   # distanza della proiezione di x nello spazio delle distriminant component
   # e più vicina alla media mu_k dove mu_k è a sua volta la proiezione in 
   # questo spazio del vettore delle medie dei gruppi
   f.class <- classify(scores = scores, cc.M = cc.M)
   if(print.result){
      print(table(classe.vera = group, classe.attr = f.class))
   }
   # compute the error of the fischer classifier
   errors <- (f.class != group)
   APERf   <- sum(errors) / n
   if(print.result){
      print("------------------")
      cat("APER = ", APERf)
   }
   
   ########### ATTRIBUTION ################
   cc.new <- as.matrix(X.new) %*% A[,1:n_dp]
   new.class <- classify(cc.new,cc.M)
   
   ########### PLOTTING ################
   if(print.plot){
      # extract the first two canonical component because we want to plot on a plane
      cc1 <- scores[,1]
      cc2 <- scores[,2]
      color.group <- group
      levels(color.group) <- rainbow(g)
      
      x11()
      plot(cc1, cc2, main='Fisher discriminant analysis', xlab='first canonical coordinate',
           ylab='second canonical coordinate', pch=20, col=as.character(color.group))
      legend(min(cc1), min(cc2)+2, legend = levels(group), fill= rainbow(g), cex=.7)
      
      for(i in 1:g){
         points(cc.M[i,1], cc.M[i,2], pch=4, col = rainbow(g)[i] , lwd=2, cex=1.5)
      }
      x.cc  <- seq(min(cc1),max(cc1),len=200)
      y.cc  <- seq(min(cc2),max(cc2),len=200)
      xy.cc <- expand.grid(cc1=x.cc, cc2=y.cc)
      
      z <- matrix(ncol = g,nrow = length(sqrt(rowSums(scale(xy.cc,cc.M[1,],scale=FALSE)^2))))
      for(i in 1:g){
         z[,i] <- sqrt(rowSums(scale(xy.cc,cc.M[i,],scale=FALSE)^2))
      }
      zi.cc <- NULL
      for(i in 1:g){
         for(j in 1:length(z[,i])){
            zi.cc[j] <- z[j,i] - min(z[j,-i])
         }
         contour(x.cc, y.cc, matrix(zi.cc, 200), levels=0, drawlabels=F, add=T)
      }
   }
   
   return(new.class)
} # end function


classify <- function(scores,cc.M){
   n <- dim(scores)[1]
   f.class <- factor(levels = levels(group))
   dist.m <- matrix(0, nrow = n, ncol = g)
   for(i in 1:n){
      for(j in 1:g){
         dist.m[i,j] <- norm(scores[i,] - cc.M[j,],type = "2")^2
      }
      f.class[i] <- levels(group)[which.min(dist.m[i,])]
   }
   return(f.class)
}# end function classify
