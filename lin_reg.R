# function for running linear regression in an compact way
# 
# INPUT:
#      X = dataframe con i regressori
#      Y = vettore della risposta
#      formula = string with the regression formula (Y ~ X1 + ...)
#      list.A = lista con le matrici per i test da fare sui parametri
#      list.b = lista con il termine noto b del test A*beta = b
#      X0.new = dataframe con le nuove osservazioni su cui fare predizione
#      interval =  NULL or c("confidence") or c("confidence","prediction"),
#                  if it's NULL no prediction is done
#      alpha.pred = alpha to be uded in prediction, if we want simultaneous interval
#                   just include the Bonferroni correction
#      remove = T se si vuole rimuove gli outlier dallla design matrix
#      sim.band = TRUE if you want a simultaneous prediction and confidence band,
#                 it will be plotted againt every regressor
#      pointwise = TRUE if you want a pontwise prediction and confidance band,
#                  it will be plotted againt every regressor
#                  
# OUTPUT:
# 
# 
# REMARK : se non voglio intervalli inferenza sui parametri lascio list.A NULL
#          se non voglio predizione du nuove osservazioni lascio interval NULL
#          se voglio predirre simultaneo su dati in input applico correzioni di
#          Bonferroni sull'alpha.pred
library(MASS)
library(car)
library(rgl)
path <- "/home/andrea/StatApp/StatApp_test"
source(paste(path,"PCA/PCA_function.R",sep = "/"))
source(paste(path,"outlier.R",sep = "/"))

lin_reg <- function(Y, 
                    X,
                    formula,
                    list.A = NULL,
                    list.b = NULL,
                    X0.new = NULL,
                    interval = NULL,
                    alpha = 0.05,
                    alpha.pred = 0.05,
                    print.plot = F,
                    print.result = T,
                    print.band = F,
                    PCA = F,
                    print.plot.DIAGN = F,
                    sim.band = F,
                    pointwise = F,
                    remove = F){
   # numerosità campionaria
   n <- dim(X)[1]
   # numero di regressori
   r <- dim(X)[2]
   # rimuove gli outlier solamente se chiesto, il grafico con la distanza di
   # malanobis lo fa sempre
   # find_outlier(X, remove = remove, plot = F)
   
   ######## DATA EXPLORATION ###########
   if(print.plot){
      # scatter plot per vedere le marginali
      x11()
      plot(data.frame(Y,X))
      # se due regressori plotta in 3D i punti per vedere la relazione
      if(r == 2){
         open3d()
         points3d(x = X[,1], y = X[,2], z = Y, size = 4,  aspect = T)
         box3d()
         axes3d()
      }
      
   }
   
   
   if(PCA){
      # calcolo gli scores che potrebbero servirmi dopo, e plotto il risultato solo
      # se rischiesto
      pc <- PC(X, print.results = F, graph = T)
   }
   
   ########## MODEL FITTING ############
   fit <- lm(as.formula(formula), data = X)
   # stampa il summary della regressione
   if(print.result){print(summary(fit))}
   # stampa per diagniostica
   if(print.plot.DIAGN){
      x11()
      plot(fit)
   }
   # plot 3D con osservazioni e piano fittato solo se 2 regressori
   if(print.plot){
      if(r == 2){
         griglia <- expand.grid(x=range(X[,1]),y=range(X[,2]))
         names(griglia) <- colnames(X)
         points3d(x=X[,1], y=X[,2], z=fitted(fit), size=4)
         points3d(x=X[,1], y=X[,2], z=Y, size=4,col = "Red")
         surface3d(range(X[,1]), 
                   range(X[,2]), 
                   matrix(predict(fit, griglia),2,2),alpha = 0.5)
      }
   }
   # stampare test di normalità sui residui, test quantitativo di quanto
   # osservato qualitativamente dai grafici
   if(print.result){shapiro.test(rstandard(fit))}
   
   ########### TEST SUI PARAMETRI MODELLO ############
   if(!is.null(list.A)){
      stopifnot(length(list.A) == length(list.b))
      # inizializzo le liste che conterranno i CI
      IC.sim <- IC.Bonf <- vector("list",length(list.A))
      # vettori dei parametri stimati del modello
      B <- coefficients(fit)
      for(i in 1:length(list.A)){
         A <- list.A[[i]]
         b <- list.b[[i]]
         print("#######################")
         print(linearHypothesis(fit, A, b))
         # build confidence intervals
         CI.sim <- matrix(0,nrow = dim(A)[1],ncol = 2)
         CI.Bonf <- matrix(0,nrow = dim(A)[1],ncol = 2)
         for(j in 1:dim(A)[1]){
            a <- A[j,]
            LB.sim <- a%*%B - sqrt(a%*%vcov(fit)%*%a)*sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            UB.sim <- a%*%B + sqrt(a%*%vcov(fit)%*%a)*sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            LB.Bonf <- a%*%B - sqrt(a%*%vcov(fit)%*%a) * qt(1-alpha/(2*dim(A)[1]),n-r-1)
            UB.Bonf <- a%*%B + sqrt(a%*%vcov(fit)%*%a) * qt(1-alpha/(2*dim(A)[1]),n-r-1)
            CI.sim[j,] <- c(LB.sim,UB.sim)
            CI.Bonf[j,] <- c(LB.Bonf,UB.Bonf)
         }
         # queste liste contengono intervalli di confidenza per le combinazioni
         # lineari dei parametri beta specificate nelle matrici A per ogni A nella
         # lista list.A
         IC.sim[[i]] <- CI.sim
         IC.Bonf[[i]] <- CI.Bonf
      }
      print("intervalli di confidenza di Bonf per i parametri specificati nelle matrici A dei test")
      print(IC.Bonf)   
      
   }
   
   ########## PREVISIONE ##################
   if(!is.null(interval)){# se interval è null non faccio predizione
      if(length(interval) == 2){
         names(X0.new) <- colnames(X)
         # predizione + CI per risposta media
         Y.conf <- predict(fit, X0.new, interval='confidence', level=1-alpha.pred) 
         # predizione + CI per nuova osservazione
         Y.pred <- predict(fit, X0.new, interval='prediction', level=1-alpha.pred)
         new.response <- list(Y.conf, Y.pred)
      }
      if(length(interval) == 1){
         # in output restituisco un dataframe dove nella prima colonna ho il valore
         # fittato e nella altre due l' il confidence o prediction interval
         print(alpha.pred)
         new.response <- predict(fit, X0.new, interval=interval, level=1-alpha.pred)
      }
   }
   # costuizione di bande di confidenza simultanee da plottare,
   if(print.band){
      # punti della griglia
      m = 100
      # inizializzazione matrice griglia, attenzione che la prima colonna è 
      # quella di 1 per beta0
      Z0 <- matrix(1,nrow = m, ncol = r+1)
      # costruire dataframe per la predizione, ogni colonna i è una griglia
      # di m punti nel range della i-esimo regressore
      for(i in 2:(r+1)){
         Z0[,i] <- seq(range(X[,i-1])[1],range(X[,i-1])[2], length = m)
      }
      # se voglio banda simultanea
      if(sim.band){
         # matrice che conterrà il punto fittato nella prima colonna e lower and
         # upper bound della fascia di confidenza simultanea delle altre 2 colonne
         Conf <- matrix(0,nrow = m, ncol = 3)
         Pred <- matrix(0,nrow = m, ncol = 3)
         # stima di sigma^2
         s2_err <- sum(residuals(fit)^2)/fit$df.residual
         for(j in 1:dim(Z0)[1]){
            z0 <- Z0[j,]
            LB.conf <- z0%*%B - sqrt(z0%*%vcov(fit)%*%z0) * sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            UB.conf <- z0%*%B + sqrt(z0%*%vcov(fit)%*%z0) * sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            LB.pred <- z0%*%B - sqrt(s2_err + z0%*%vcov(fit)%*%z0) * sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            UB.pred <- z0%*%B + sqrt(s2_err + z0%*%vcov(fit)%*%z0) * sqrt((r+1)*qf(1-alpha,r+1,n-r-1))
            Conf[j,] <- c(z0%*%B, LB.conf, UB.conf)
            Pred[j,] <- c(z0%*%B, LB.pred, UB.pred)
         }
         
         # plotta le bande di confidenza per ogni marginale, ossia risposta vs
         # regressore, per ogni regressore
         for(i in 1:r){# scorro su tutti i regressori
            x11()
            # plotto regressore i vs risposta
            plot(data.frame(X[,i],Y), xlab = colnames(X)[i],
                 ylab = "response", las = 1, xlim=range(X[,i]), ylim=range(Pred))
            # devo prendere la conolla i+1-esima perche la prima di Z0 è
            # quella di 1
            lines(Z0[,i+1], Conf[,1], lty=2, col='blue', lwd=3) # valori fittati
            lines(Z0[,i+1], Conf[,2], lty=2, col='red', lwd=2) # lower bound
            lines(Z0[,i+1], Conf[,3], lty=2, col='red', lwd=2) # upper bound
            
            lines(Z0[,i+1], Pred[,2], lty=3, col='gold', lwd=2) # lower bound
            lines(Z0[,i+1], Pred[,3], lty=3, col='gold', lwd=2) # upper bound
            title(main = "banda confidenza e predizione simultanea",
                  sub = paste("livello confidenza = ", 1-alpha, sep = " "))
            legend('bottomright',c('conf', 'pred','fit'), lty=1, col=c('red','gold','blue'))
         }
      }
      if(pointwise){# voglio banda one-at-time, per avere la banda come 
         # output mi basta mettere una griglia continua in input in X0.new
         # tolgo la colonna di 1 che usando predict non mi serve
         Z0   <- data.frame(Z0[,2:(r+1)])
         names(Z0) <- colnames(X)
         Conf <- predict(fit, Z0, interval='confidence')
         Pred <- predict(fit, Z0, interval='prediction')
         for(i in 1:r){
            x11()
            plot(data.frame(X[,i],Y), xlab = colnames(X)[i],
                 ylab = "response", las = 1, xlim=range(X[,i]), ylim= range(Pred))
            lines(Z0[,i], Conf[,'fit'], lty=2, col='blue',lwd=3)
            lines(Z0[,i], Conf[,'lwr'], lty=2, col='red', lwd=2)
            lines(Z0[,i], Conf[,'upr'], lty=2, col='red', lwd=2)
            
            lines(Z0[,i], Pred[,'lwr'], lty=3, col='gold', lwd=2)
            lines(Z0[,i], Pred[,'upr'], lty=3, col='gold', lwd=2)
            title(main = "banda confidenza e predizione point-wise",
                  sub = paste("livello confidenza = ", 1-alpha, sep = " "))
            legend('bottomright',c('conf', 'pred','fit'), lty=1, col=c('red','gold','blue'))
         }
      }
   }
   
   
   # cosa sto restituendo?
   # model = output di lm, pssia tutte le info sulla regressione
   # PC = tutte le info sulla PCA sui regressori
   # IC = intervalli di confidenza simultanei sui parametri per ogni test fatto
   # prediction = intervalli di confidenza o di predizione su nuove osservazioni
   # prediction.sim = intervalli di confidenza o predizioni su nuove ossrvazioni
   #                  SIMULTANEI
   return(list(model = fit,
               PC = myifelse(PCA,pc,"no PCA"),
               IC = myifelse(!is.null(list.A),list(IC.Bonf = IC.Bonf,IC.sim = IC.sim),"no inferenza"),
               prediction = myifelse(!is.null(interval),new.response,"no predizione"),
               myifelse(sim.band,list(Conf.sim = Conf,Pred.sim = Pred),"no predizione simulatanea")))
   
}# end function


myifelse <- function(condition,x,y){
   if(condition)
      return(x)
   else
      return(y)
   
}# end function













