# test per l'uguaglianza delle medie di due popolazioni normali indipendenti.
# Il fatto che mi spinge a scegliere questo test rispetto ad un paired è che 
# magari ho diversa numerosità campionaria n1 != n2.
# 
# INPUT :
#        X1 = dataframe popolazione 1
#        X2 = dataframe popolazione 2
#        print.plot = TRUE se si vogliono i plot
#        delta0 = ipotesi nulle del test
#        A = matrice con le righe le combinazioni linear su cui si vuole
#            calcolare l'intervallo di coonfidenza per la differenza delle medie
# 
# 
# OUTPUT:
#        IC = lista che contiene gli intervalli di confidenza per la differenza
#             tra le medie delle due popolazioni sia simulatanei che di
#             Bonferroni
# 
# NOTE:
# 
two.pop_test <- function(X1, X2, delta0, alpha = 0.05, print.plot = T, A = NULL){
   
   # numerosità popolazione 1
   n1 <- dim(X1)[1]
   # numerosità popolazione 2
   n2 <- dim(X2)[1]
   # dimensione popolazione 1
   p1 <- dim(X1)[2]
   # dimensione popolazione 2
   p2 <- dim(X2)[2]
   if(p1 != p2){
      stop("population dimensions don't agree")
   }
   p <- p1
   # se le due popolazioni sono univariate allora faccio semplicemnte un t-test
   # per il confronto delle medie delle due popolazioni assumendo uguale varianza
   if(p == 1){
      p.gauss <- c(shapiro.test(X1[,1])$p.value,shapiro.test(X2[,1])$p.value)
      z <- t.test(x = X1, y = X2, var.equal = T,conf.level = 1-alpha)
      print("###### TEST UNIVARIATO UGUAGLIANZA MEDIE POP NORM INDIP ###########")
      print(z)
      return(list(test = z, p.gauss = p.gauss))
   }
   # calcolo la media campionaria per i due campioni
   M1 <- sapply(X1,mean)
   M2 <- sapply(X2,mean)
   # calcolo matrice di covarianza campionaria per i due campioni
   S1  <-  cov(X1)
   S2  <-  cov(X2)
   # calcolo la matrice di covarianza pool: uguale per hp nei due campioni
   S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
   
   ######## checking the assumptions ###########
   # 1) normality
   check_gaussianity(X1, alpha,1)
   check_gaussianity(X2, alpha, 2)
   #2) uguaglianze delle matrici S1 e S2
   
   ######## test and CI ###############
   # inversa della matrice S pool
   Sinv <- solve(S)
   # valore calcolato della statistica
   T2 <- n1 * n2 /(n1 + n2) * (M1 - M2 - delta0) %*% Sinv %*% (M1 - M2 - delta0)
   # quantile della distribuzione
   cfr.fisher <- (n1+n2-2)*p / (n1+n2-1-p) * qf(1-alpha,p,n1+n2-1-p)
   # per correggere con bonferroni uso p perche costruisco tanti CI
   # quante componenti hi nel vettore quindi p
   cfr.t <- qt(1 - alpha/(2 * p), n1 + n2 - 2)
   # calcolo p.value del test
   p.value <- 1 - pf( (n1+n2-1-p) / ((n1+n2-2)*p) * T2, p, n1+n2-1-p )
   print("############## TWO INDIPENDENT NORMAL POPULATION TEST ##########")
   cat("@ conf level alpha = ",alpha,ifelse(p.value<alpha,"we reject H0","we cannot reject H0"),"\n")
   cat("T2 test p-value = ", p.value,"\n")
   cat("statistica calcolata T2.0 = ",T2, "\n")
   cat("quantile della F * coeff = ", cfr.fisher, "\n")
   
   # INTERVALLI CONFIDENZA DIFFERENZA MEDIE DUE POPOLAZIONI
   # intervalli di confidenza pr la differenza tra medie delle due popolazioni sia
   # simultanei che bonferroni
   IC.sim <- cbind(M1 - M2 -sqrt(cfr.fisher*diag(S)*(1/n1+1/n2)),
                   M1 - M2 +sqrt(cfr.fisher*diag(S)*(1/n1+1/n2)))
   colnames(IC.sim)<-c("inf","sup")
   IC.Bonf <- cbind(M1 - M2 - cfr.t * sqrt(diag(S)*(1/n1+1/n2)),
                    M1 - M2 + cfr.t * sqrt(diag(S)*(1/n1+1/n2)))
   colnames(IC.Bonf)<-c("inf","sup")
   IC <- list(simult = IC.sim, Bonf = IC.Bonf)
   print("############# IC DIFFERENZA MEDIE, BONFERRONI E T2 #################")
   cat("intervalli di confidenza al livello alpha =",alpha, " per differenza medie componenti","\n")
   print(IC)
   
   # INTERVALLI DI CONFIDENZA PER COMBINAZIONI LINEARI
   if(!is.null(A)){
      IC.sim.2 <- cbind(A%*%(M1 - M2) -sqrt(cfr.fisher*diag(A%*%S%*%t(A))*(1/n1+1/n2)),
                        A%*%(M1 - M2) +sqrt(cfr.fisher*diag(A%*%S%*%t(A))*(1/n1+1/n2)))
      colnames(IC.sim.2)<-c("inf","sup")
      k <- dim(A)[1]
      cfr.t <- qt(1 - alpha/(2 * k), n1 + n2 - 2)
      IC.Bonf.2 <- cbind(A%*%(M1 - M2) - cfr.t * sqrt(diag(A%*%S%*%t(A))*(1/n1+1/n2)),
                         A%*%(M1 - M2) + cfr.t * sqrt(diag(A%*%S%*%t(A))*(1/n1+1/n2)))
      colnames(IC.Bonf.2)<-c("inf","sup")
      IC.2 <- list(simult = IC.sim.2,Bonf = IC.Bonf.2)
      print("############# IC LINEAR COMB DIFFERENZA MEDIE, BONFERRONI E T2 #################")
      cat("intervalli di confidenza al livello alpha =",alpha, " per differenza medie componenti","\n")
      cat("correzione di bonferroni k = ",k, "\n")
      print(IC.2)
      
      return(list(componenti = IC, linear_comb = IC.2))
      
   }
   
   
   
   return(IC)
}# end function




check_gaussianity <- function(X, alpha,i){
   # i indica il gruppo su cui sto controllando l'hp
   
   # check for gaussianity
   mctest <- mcshapiro.test(X)
   if(mctest$pvalue < alpha)
      warning("gaussian hp not supported by data for pop ", i, " ,shapiro p-value = ",mctest$pvalue,"\n")
   else
      cat("You are assuming gaussianity, there's no evidence to say the contrary, shapiro p-value = ",mctest$pvalue,"\n")
   
} # end function check_gaussianity
