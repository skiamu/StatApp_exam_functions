########### MANOVA ONE-WAY ##########
#
#
# funzione per lo studio della MANOVA ONE-WAY
#
#
# INPUT: 
#       Y = dataframe della risposta, la risposta è vettoriale 
#       treat = fattore 1, stessa lunghezza di Y con livelli
#               per ogni osservazione
#       int.diff = T se si vogliono gli intervalli di confidenza per le differenze
#                  delle medie nei gruppi del fattore 1
#      int.mean = T se si vuole gli intervalli di confidenza per le medie
#                 nei gruppi
#                 
# NOTE : manca da fare il check sull'ugluaglianza delle covarianze
# e mettere a posto l'ultimo plot
#
oneway_manova <- function(Y,
                          treat,
                          int.diff = T,
                          int.mean = T,
                          plot.result = F,
                          alpha = 0.05,
                          print.result = T){
   
   # dimensione della risposta
   p <- dim(Y)[2]
   # numerosità in ciascun gruppo (livello)
   ng <- table(treat)
   # numerosità campionaria
   n <- dim(Y)[1]
   # numero di livelli del fattore
   g <- length(levels(treat))
   
   ####### DATA EXPLORATION #################
   
   colori <- NULL
   for(i in 1:g){
      colori[which(treat == levels(treat)[i])] <- rainbow(g)[i]
   }
   if(plot.result){
      x11()
      pairs(Y, col = colori, pch=16)
   }
   # plotto la variabilità per ogni livello del fattore delle componenti
   # della risposta. Questo plot è da guardare trasversalmente osservando
   # se ci sono differenze significative tra box dello stesso colore, in quanto
   # sotto H0 queste differenze non sono statisticamente significative
   if(plot.result){
      x11()
      par(mfrow = c(1,g), las = 2)
      for(i in 1:g){
         idx <- which(treat == levels(treat)[i])
         boxplot(Y[idx,],main = levels(treat)[i], ylim = range(Y), col = rainbow(p))
      }
   }
   # plotto la variabilità per ogni componente della risposta tra i livelli
   # del fattore. Va guardato riquadro per riguadro, sotto H0 non c'è differenza
   # tra i box dentro un riquadro, immaginando di rifiutare H0 posso 
   # farmi un'ideadi quale componente provova il rifiuto individuando
   # qualle dove i box sono più diversi.
   if(plot.result){
      x11()
      par(mfrow = c(1,p), las = 2)
      for(i in 1:p){
         boxplot(Y[,i] ~ treat, main = colnames(Y)[i], ylim = range(Y), col = rainbow(g))
      }
   }
   
   ########## CHECK ASSUMPTIONS ################
   # controlla la normalità p-dimensionale della risposta
   Ps <- NULL
   for(i in 1:g){
      Ps <- c(Ps,mcshapiro.test(as.matrix(Y[which(treat == levels(treat)[i]),]))$p)
   }
   if(print.result){
      if(any(Ps < alpha)){
         warning("gaussianità non supportata dai dati")
      }
   }
   # controlla che la matrice di covarianza della risposta sia sta stessa
   # dentro ogni gruppo del fattore
   
   
   ########## MODEL FITTING ################
   
   # se rifiuto questo test ho evidenza per dire che il trattamento ha 
   # effetto sul vettore della risposta media
   fit <- manova(as.matrix(Y) ~ treat)
   z <- summary.manova(fit,test="Wilks")
   print("###################### SUMMARY MANOVA ############################")
   print(summary.manova(fit,test="Wilks"))
   
   # Si vuole capire quale componente del vettore risposta ha causato
   # il rifiuto del test (se rifiuto c'è stato). Faccio p ANOVA
   print("################### SUMMARY ANOVA SULLE COMPONENTI ################")
   print(summary.aov(fit))
   
   # Si vuole capire quale livello inflenza effettivamnte la media delle
   # componenti della risposta. Costruisco intervalli di confidenza
   # per le differenze tre medie per ogni compenente del vettore 
   # delle risposta
   
   # stima della matrice covarianza
   W <- summary.manova(fit)$SS$Residuals
   # grand mean
   M <- sapply(Y,mean)
   # matrice delle medie nei gruppi, la media è un vettore riga
   Mg <- matrix(nrow = g, ncol = p)
   for(i in 1:g){
      Mg[i,] <- sapply(Y[which(treat == levels(treat)[i]),], mean)
   }
   treat_lev <- levels(treat)
   
   if(int.diff){
      # numero di differenze tra medie (g-1) * g / 2 per ogni componente
      k <- p*g*(g-1)/2
      # quantile della t student per Bonferroni CI
      qT <- qt(1-alpha/(2*k), n-g)
      # CI è una lista che contiene in ogni entrata ha gli intervalli per le
      # differenze tra gruppi, per ogni componente del vettore delle risposte.
      # Quindi è lungo quanto le varie differenze che si possono formare
      CI <- vector("list",g*(g-1)/2)
      nomi.diff <- vector("character",g*(g-1)/2)
      count <- 0
      for(i in 1:(g-1)){
         for(j in (i+1):g){
            count <- count + 1
            nomi.diff[count] <- paste0(treat_lev[i],"-",treat_lev[j])
            inf <- Mg[i,]-Mg[j,] - qT * sqrt( diag(W)/(n-g) * (1/ng[i]+1/ng[j]) )
            sup <- Mg[i,]-Mg[j,] + qT * sqrt( diag(W)/(n-g) * (1/ng[i]+1/ng[j]) )
            CI[[count]] <- cbind(inf,sup)
         }
      }
      names(CI) <- nomi.diff
      print("######## CI DIFFERENZE MEDIE NEI GRUPPI PER COMPONENTI ############")
      cat("la correzione di Bonferroni usata è k = ", k, "\n")
      cat("alpha = ",alpha, "\n")
      print(CI)
   }
   if(int.mean){
      k <- p * g
      qT <- qt(1-alpha/(2*k),n-g)
      CI.mean <- vector("list",g)
      for(i in 1 :g){
         inf <- Mg[i,] - qT * sqrt(diag(W)/(n-g)* (1/ng[i]))
         sup <- Mg[i,] + qT * sqrt(diag(W)/(n-g)* (1/ng[i]))
         CI.mean[[i]] <- cbind(inf,sup)
      }
      names(CI.mean) <- treat_lev
      print("######## CI  MEDIE NEI GRUPPI PER COMPONENTI ############")
      cat("la correzione di Bonferroni usata è k = ", k, "\n")
      cat("alpha = ", alpha, "\n")
      print(CI.mean)
      
   }
   
   return(list(CI = myifelse(int.diff,CI,"no intervalli differenze"),
               CI.mean = myifelse(int.mean,CI.mean,"no intervalli medie"),
               summary.manova = z,
               p.gauss = Ps))
   
} # end function



# INPUT :
#        Y = matrice nxp della risposta
#        X = dataframe nx2 dei fattori



twoway_manova <- function(Y,X,alpha = 0.05,print.result = T,with.interaction = F){
   
   # primo fattore
   F1 <- X[,1]
   # secondo fattore
   F2 <- X[,2]
   # dimensione della risposta
   p <- dim(Y)[2]
   # numero osservazioni dentro ogni gruppo del primo fattore
   ng <- table(F1)
   # numero osservazioni dentro ogni gruppo del secondo fattore
   nb <- table(F2)
   # numero livelli primo fattore
   g <- length(levels(F1))
   # numero livelli secondo fattore
   b <- length(levels(F2))
   # numero osservazioni dentro ogni rettangolo
   if(length(F1) %% g*b == 0){
      n <- length(F1) %% g*b
   }else{
      warning("unbalanced design")
   }
   # livelli fattore 1
   treat_lev1 <- levels(F1)
   # livelli fattore 2
   treat_lev2 <- levels(F2)
   ################# check assumption gaussianity ###################
   P.gauss <- vector("double",g*b)
   count <- 0
   for(i in 1:g){
      for(j in 1:b){
         count <- count + 1
         idx <- which(F1 == treat_lev1[i] & F2 == treat_lev2[j])
         P.gauss[count] <- mcshapiro.test(Y[idx,])$p
      }
   }
   
   if(with.interaction && n != 1){
      fit <- manova(as.matrix(Y) ~ F1 + F2 + F1:F2)
      z <- summary.manova(fit, test="Wilks")
      print("############ SUMMARY TWO-WAY MANOVA WITH INTERACTION ###########")
      print(z)
   }else{
      fit <- manova(as.matrix(Y) ~ F1 + F2)
      z <- summary.manova(fit, test="Wilks")
      print("############ SUMMARY TWO-WAY MANOVA ADDITIVE MODEL ###########")
      print(z)
   }
   
   print("################### SUMMARY MARGINAL ANOVA #########################")
   print(summary.aov(fit))
   
   # INTERVALLI DI CONFIDENZA
   # numero di intervalli da controllare
   k <- g*(g-1)/2*p + b*(b-1)/2*p
   # gradi di libertà residui (da usare negli intervalli)
   dof.res <- fit$df.residual
   # inizializzazione
   IC <- vector("list",p)
   # quantile della t student per Bonferroni CI
   qT <- qt(1-alpha/(2*k), dof.res)
   # stima della matrice covarianza: diversa dal caso della manova ma non si sa il motivo
   W <- t(fit$residuals) %*% fit$residuals 
   nomi.Y <- colnames(Y)
   for(t in 1:p){
      IC.F1 <- data.frame(nrow = g*(g-1)/2, ncol = 2)
      IC.F2 <- data.frame(nrow = b*(b-1)/2, ncol = 2)
      
      # calcolo medie dentro i livelli per i due fattori
      Mg <- by(Y[,t],F1,mean)
      Mb <- by(Y[,t],F2,mean)
      count <- 0
      nomi <- vector("character",g*(g-1)/2)
      for(i in 1:(g-1)) {
         for(j in (i+1):g) {
            count <- count + 1
            nomi[count] <- paste(treat_lev1[i],"-",treat_lev1[j])
            IC.F1[count,] <-
               c(Mg[i]-Mg[j] - qt(1-alpha/(2*k), dof.res) * 
                    sqrt( diag(W)[t] * ( 1/ng[i] + 1/ng[j] )/dof.res),
                 Mg[i]-Mg[j] + qt(1-alpha/(2*k), dof.res) * 
                    sqrt( diag(W)[t] * ( 1/ng[i] + 1/ng[j] )/dof.res))
         }
      }
      rownames(IC.F1) <- nomi
      colnames(IC.F1) <- c("LB","UB")
      count <- 0
      nomi <- vector("character",b*(b-1)/2)
      for(i in 1:(b-1)) {
         for(j in (i+1):b) {
            count <- count + 1
            nomi[count] <- paste(treat_lev2[i],"-",treat_lev2[j])
            # inserisco separatamente perche non posso fare un vettore misto
            IC.F2[count,] <-
               c(Mb[i]-Mb[j] - qt(1-alpha/(2*k), dof.res) * 
                    sqrt( diag(W)[t] * ( 1/nb[i] + 1/nb[j] )/dof.res),
                 Mb[i]-Mb[j] + qt(1-alpha/(2*k), dof.res) * 
                    sqrt( diag(W)[t] * ( 1/nb[i] + 1/nb[j] )/dof.res))
         }
      }
      rownames(IC.F2) <- nomi
      colnames(IC.F2) <- c("LB","UB")
      IC[[t]] <- list(IC.F1 = IC.F1,IC.F2 = IC.F2)
   }
   names(IC) <- nomi.Y
   print("############### IC DIFFERENZA MEDIE SULLE COMPONENTI ###############")
   cat("Si è usata una correzione di Bonferroni k = ",k, "\n")
   for(i in IC){print(i)}
   
   
   return(list(IC = IC,
               fit = fit,
               P.gauss = P.gauss))
   
   
   
}