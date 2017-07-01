############ ONE-WAY ANOVA ###################
#
# funzione per lo studio dell'ANOVA ONE-WAY

#
# INPUT: 
#       X = dataframe come le seguenti colonne
#           1) Y = vettore delle risposte
#           2) treat = factor vector della stessa lunghezza di Y 
#                      con i trattamenti somministrati alla stat unit
#                      per cui si è registrata la risposta
# 
# OUTPUT:
#       ICrange = dataframe con gli intervalli di confidenza per tutte
#                 le differenze tra i tau_i. Una volta dimostrato che il
#                 fattore ha effetto sulla risposta media mi serve per
#                 capire quale livello/i è/sono quelli determinanti
# 
# 
oneway_anova <- function(X,plot.result = T, alpha = 0.05, print.result = T,var.int = F){
   
   # estraggo dal dataframe vettore delle risposte e dei trattamenti
   Y <- X[,1]; treat <- X[,2]
   
   
   # This is a case of one-way ANOVA: one variable (Y) observed 
   # over g levels
   n       <- length(treat)      # total number of obs.
   ng      <- table(treat)       # number of obs. in each group
   treat_lev   <- levels(treat)      # levels of the treatment
   g       <- length(treat_lev)     # number of levels (i.e., of groups)
   
   ##### ASSUMPTION VERIFICATION  #####
   # lo posso fare a priori oppure a posteriori sui residui. Le hp da 
   # verificare sono le seguenti:
   # 1) normalità dentro ciascun gruppo
   # 2) uguaglianza delle g matrici di covarianza
   
   # verifico 1)
   Ps <- NULL
   for(i in 1:g){
      Ps <- c(Ps,shapiro.test(Y[ treat == treat_lev[i] ])$p)
   }
   if(any(Ps < alpha)){# se almeno un gruppo ha p.value sotto alpha
      # indici dei gruppi in cui l'ipotesi di gaussianità non è vera
      idx.noGauss <- which(Ps < alpha)
      for(i in 1:length(idx.noGauss)){
         warning("Gaussian hp not verified for groups: ", treat_lev[idx.noGauss[i]],
                 " p.value = ", Ps[idx.noGauss[i]])
      }
   }
   # verifico 2)
   # test of homogeneity of variances
   # H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
   # H1: there exist i,j s.t. sigma.i!=sigma.j
   # sono felice di accettare, come nel caso di gaussianità, l'ipotesi nulle è 
   # invertita rispetto al caso delle scoperte
   cov.test <- bartlett.test(Y, treat)
   if(cov.test$p.value < alpha){
      warning("same covariance hp not verified, p.value = ",cov.test$p.value )
   }
   
   ####### ANOVA FITTING ######
   fit <- aov(Y ~ treat)
   print("################ SUMMARY ANOVA ONE-WAY #########################")
   print(summary(fit))
   
   # se p.value basso significa che rifiutiamo H0 ossia c'è evidenza per
   # affermare che i trattementi hanno effetto sulla media della risposta.
   # Sono interessato a sapere quale tra i g trattamenti ha provocato il rifiuto
   # dell'hp nulla. Devo costruire intervalli di confidenza per la differenza 
   # tra medie nei gruppi
   
   # numero di differenze che posso costruire con un vettore a g componenti.
   # Aumento di uno se voglio IC anke per la varianza
   k <- ifelse(var.int, g*(g-1)/2 + 1, g*(g-1)/2)
   # grand mean
   Media   <- mean(Y)
   # la funzione tapply applica una funzione dentro ogni gruppo,
   # il primo argomento è la risposta il secondo il vettore coi fattori
   Mediag  <- tapply(Y, treat, mean)
   # lo posso anke prendere dal summari di aov
   SSres <- sum(residuals(fit)^2)
   # stima della varianza sigma^2
   S <- SSres/(n-g)
   
   IC <- data.frame(nrow = g*(g-1)/2,ncol = 2)
   nomi <- vector("character",g*(g-1)/2)
   count <- 0
   for(i in 1:(g-1)) {
      for(j in (i+1):g) {
         count <- count + 1
         nomi[count] <- paste(treat_lev[i],"-",treat_lev[j])
         IC[count,] <-
            c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * 
                 sqrt( S * ( 1/ng[i] + 1/ng[j] )),
              Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * 
                 sqrt( S * ( 1/ng[i] + 1/ng[j] )))
      }
   }
   rownames(IC) <- nomi
   colnames(IC) <- c("LB","UB")
   if(var.int){
      IC.var <- c(LB = (n-g) * S / qchisq(1-alpha/(2*k),n-g),
                  UB = (n-g) * S / qchisq(alpha/(2*k),n-g))
   }
   
   # stampo gli intervalli di confidenza per ogni differenza delle medie
   print("############## IC DIFFERENZE MEDIE (tau_i - tau_j) #################")
   cat("correzione di Bonferroni k = ", k, "\n")
   print(IC)
   if(var.int){cat("IC per la varianza : ", IC.var, "\n")}
   
   
   # STIMA MEDIE CON IL MODELLO ANOVA
   tau <- Mediag - Media
   M.hat <- Media + tau
   
   
   
   return(list( IC = IC,
                fit = fit, 
                assunzioni = list(gauss = Ps,var = cov.test),
                M.hat = M.hat,
                tau = tau))
   
} # end function

############ TWO-WAY ANOVA ###############
#
# INPUT :
#       X = dataframe con le seguenti colonne
#           Y = vettore numerico delle risposte
#           F1 = vettore di fattori (fattore 1)
#           F2 = vettore di fattori (fattore 2)
#       plot.result = TRUE se si vogliono plot
#       print.result = TRUE se si voglioni i risultati numerici
#       formula = formula per fittare il modello di ANOVA
#             
# OUTPUT:
#       IC:F1 = dataframe con gli intervalli di confidenza simultanei delle 
#               differenze di medie nel primo fattore
#       IC.F2 = dataframe con gli intervalli di confidenza simultanei delle 
#               differenze di medie nel secondo fattore
#       Pcov = p.value del bartlett test
#       Ps = p. value del shapiro test
#              
# OSS : i CI in output mi servono, una volta dimostrato che uno dei due fattori o 
#       entrambi sono significativi per la media della risposta, 
#       per capire effettivamente quale livelli di un fattore determina
#       un effetto significativo sulle risposta media
#       

twoway_anova <- function(X,
                         plot.result = T,
                         alpha = 0.05, 
                         print.result = T,
                         formula = NULL,
                         var.int = F){
   # vettore delle risposte
   Y <- X[,1]
   # primo fattore
   F1 <- X[,2]
   # secondo fattore
   F2 <- X[,3]
   # numerosità del campione
   n <- length(Y)
   # numero osservazioni dentro ogni gruppo del primo fattore
   ng <- table(F1)
   # numero osservazioni dentro ogni gruppo del secondo fattore
   nb <- table(F2)
   # numero livelli primo fattore
   g <- length(levels(F1))
   # numero livelli secondo fattore
   b <- length(levels(F2))
   
   ##### ASSUMPTION VERIFICATION  #####
   # lo posso fare a priori oppure a posteriori sui residui. Le hp da 
   # verificare sono le seguenti:
   # 1) normalità dentro ciascun quadrato
   # 2) uguaglianza delle g matrici di covarianza
   Ps <- NULL
   enough.data <- T
   for(i in 1:length(levels(F1))){
      for(j in 1:length(levels(F2))){
         idx <- intersect(which(F1 == levels(F1)[i]),which(F2 == levels(F2)[j]))
         if(length(Y[idx]) > 3)
            Ps <- c(Ps,shapiro.test(Y[idx])$p)
         else{
            warning("non ho abbastanza dati nel quadrato ",levels(F1)[i],
                    " di F1 ", levels(F2)[i]," di F2 per controllare gaussianità")
            enough.data <- F
         }
      }
   }
   if(any(Ps < alpha)){# se almeno un quadrato ha p.value sotto alpha
      # indici dei gruppi in cui l'ipotesi di gaussianità non è vera
      idx.noGauss <- which(Ps < alpha)
      for(i in 1:length(idx.noGauss)){
         warning("Gaussian hp not verified for groups: ", i,
                 " p.value = ", Ps[idx.noGauss[i]])
      }
   }
   
   # verifico 2)
   # test of homogeneity of variances
   # H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
   # H1: there exist i,j s.t. sigma.i!=sigma.j
   # sono felice di accettare, come nel caso di gaussianità, l'ipotesi nulle è 
   # invertita rispetto al caso delle scoperte
   if(enough.data){
      w <- bartlett.test(Y,F1:F2)
      p.bartlet <- w$p.value
      print("############# BARTLETT TEST ##################")
      print(w)
   }
   if(p.bartlet < alpha){
      warning("non posso assumere sigma uguali")
   }
   
   
   ###### fitting the models #############
   
   # comincio con il modello completo per vedere se interazione significativa
   fit <- aov(Y ~ F1 + F2 + F1:F2, data = X)
   print("################## AOV SUMMARY ##############################")
   print(summary(fit))
   # faccio un test simultaneo sui due fattori
   M <- mean(Y) # media intero campione
   Mg <- tapply(Y, F1, mean) # media dentro fattore 1
   Mb <- tapply(Y, F2, mean) # media dentro fattore 2
   SS_F1 <- sum(ng %*% (Mg - M)^2)              # or from the summary: 1.53    
   SS_F2 <- sum(nb %*% (Mb  - M)^2)              # or from the summary: 66.70
   SSres   <- sum((fit$residuals)^2) # or from the summary: 16.37
   S <- SSres/fit$df.residual
   Ftot      <- ( (SS_F1 + SS_F2) / ((g-1)+(b-1)))/(SSres / (fit$df.residual))
   Ptot      <- 1 - pf(Ftot, (g-1)+(b-1), fit$df.residual) # attention to the dgf!
   if(print.result){
      print("################ TEST SIMULTANEO SUI DUE FATTORI #################")
      cat("test simultaneo sui due  fattori, p.value = ", Ptot,"\n")
   }
   
   
   # INTERVALLI DI CONFIDENZA SULLE DIFFERENZE DI MEDIE
   # gli intervalli sono costruiti sul modello dato in input in formula, fa 
   # differenza perchè cambiano i gradi di liberà nella stima della sigma
   fit <- aov(formula,data = X)
   treat_lev1 <- levels(F1)
   treat_lev2 <- levels(F2)
   # gradi di libertà che saranno usati nel quantile delle t
   dof.res <- fit$df.residual
   # la correzione di Bonf è pari al numero di differenze che si possono fare.
   # Nel caso specifichi in input anke l'intervallo per la varianza si aggiunge 1
   k <- ifelse(var.int, g*(g-1)/2 + b*(b-1)/2 + 1, g*(g-1)/2 + b*(b-1)/2)
   IC.F1 <- data.frame(nrow = g*(g-1)/2,ncol = 2)
   IC.F2 <- data.frame(nrow = b*(b-1)/2, ncol = 2)
   nomi.g <- vector("character",g*(g-1)/2)
   nomi.b <- vector("character",b*(b-1)/2)
   count <- 0
   
   for(i in 1:(g-1)) {
      for(j in (i+1):g) {
         count <- count + 1
         nomi.g[count] <- paste(treat_lev1[i],"-",treat_lev1[j])
         IC.F1[count,] <-
            c(Mg[i]-Mg[j] - qt(1-alpha/(2*k), dof.res) * 
                 sqrt( S * ( 1/ng[i] + 1/ng[j] )),
              Mg[i]-Mg[j] + qt(1-alpha/(2*k), dof.res) * 
                 sqrt( S * ( 1/ng[i] + 1/ng[j] )))
      }
   }
   rownames(IC.F1) <- nomi.g
   colnames(IC.F1) <- c("LB","UB")
   count <- 0
   for(i in 1:(b-1)) {
      for(j in (i+1):b) {
         count <- count + 1
         nomi.b[count] <- paste(treat_lev2[i],"-",treat_lev2[j])
         IC.F2[count,] <-
            c(Mb[i]-Mb[j] - qt(1-alpha/(2*k), dof.res) * 
                 sqrt( S * ( 1/nb[i] + 1/nb[j] )),
              Mb[i]-Mb[j] + qt(1-alpha/(2*k), dof.res) * 
                 sqrt( S * ( 1/nb[i] + 1/nb[j] )))
      }
   }
   rownames(IC.F2) <- nomi.b
   colnames(IC.F2) <- c("LB","UB")
   
   if(var.int){
      IC.var <- c(LB = dof.res * S / qchisq(1-alpha/(2*k),dof.res),
                  UB = dof.res * S / qchisq(alpha/(2*k),dof.res))
   }
   print("################### IC DIFFERENZE MEDIE #############################")
   print(summary(fit))
   cat("la correzione di Bonferroni usata è k = ", k,"\n")
   print(IC.F1)
   print(IC.F2)
   if(var.int){cat("IC per la varianza : ", IC.var, "\n")}
   
   
   # STIMA MEDIE CON IL MODELLO TWO-WAY ANOVA
   # grand mean
   X.grand <- mean(Y)
   # vettore dei tau
   tau <- Mg - X.grand
   # vettore dei beta 
   beta <- Mb - X.grand
   M.hat <- data.frame(matrix(nrow = g,ncol = b))
   rownames(M.hat) <- treat_lev1
   colnames(M.hat) <- treat_lev2
   for(i in 1:g){
      for(j in 1:b){
         M.hat[i,j] <-  X.grand + tau[i] + beta[j]
      }
   }
   
   
   
   return(list(IC.F1 = IC.F1,
               IC.F2 = IC.F2,
               p.bartlet = p.bartlet,
               p.gauss = Ps,
               fit = fit,
               M.hat = M.hat,
               param = list(tau = tau,beta = beta)))
   
   
   
}# end function 
