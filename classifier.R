load("/home/andrea/StatApp/StatApp_test/inference/mcshapiro.test.RData")


# INPUT: 
#       X = vettore delle caratteristiche
#       group = vettore di fattori dei gruppi
#       jittering = T se si vuole scuotere un po i dati
#       prior = vettore delle prior probability
#       mis.cost = costi di misclassificazione
#       X0.new = dataframe con le nuove osservazioni da predirre
#       kind = "LDA" o "QDA"
classifier <- function(X,
                       group,
                       jittering = F,
                       prior = NULL,
                       mis.cost = NULL,
                       print.result = T,
                       print.plot = T,
                       X0.new = NULL,
                       kind = "LDA",
                       sd_jitt = NULL){
   
   
   # numerosità campionaria
   n <- dim(X)[1]
   # dimensione delle caratteristiche
   p <- dim(X)[2]
   if(p < 2){stop("p = 1, usare altro")}
   # numerosità campionaria per ogni gruppo
   ng <- table(group)
   # numero di gruppi
   g <- length(ng)
   # media dentro ogni gruppo (medie vettore colonna)
   Mg <- matrix(ncol = g,nrow = p)
   for(i in 1:g){
      Mg[,i] <- colMeans( X[which(group == levels(group)[i]),] )
   }
   if(length(group) != n){stop("errore dimensione vettore dei fattori")}
   if(!is.null(prior) &&  length(prior) != g){stop("errore dimensione prior")}
   if(!is.null(mis.cost) &&  length(mis.cost) != g){stop("errore dimensione costi misclass")}
   
   ###### VERIFICA ASSUNZIONI 
   p.gauss <- vector("double",g)
   for(i in 1:g){
      # indice elementi del gruppi i
      idx <- which(group == levels(group)[i])
      p.gauss[i] <- shapiro.test(X[idx,])$p.value
   }
   if(p == 1 && kind == "LDA"){
      bar.test <- bartlett.test(X[,1],group)
   }
   
   set.seed(280787)
   # smuovere un po i dati se le misure sono state fatte con arronondamneto
   if(jittering){
      X <- X + cbind(rnorm(n, sd = sd_jitt))
   }
   if(is.null(prior)){# stimo le prior dalle frequenze del campione
      prior <- ng / n
      prior.sample = T
   }else{
      prior.sample = F
   }
   
   # se dei costi di misclassificazione sono dati in input allora includili
   # nelle prior per usare i comando lda e qda.
   # Fare attenzione che effettivamente il primo costo deve andare a moltiplicare
   # la prima prior e non debbano essere scambiati
   if(!is.null(mis.cost)){
      prior <- prior * mis.cost / sum(mis.cost * prior)
      stopifnot(sum(prior) == 1)
      print("#################### AVVISO MODIFICHE PRIOR #####################")
      cat("nel calcolo delle prior modificate sto moltiplicando la prior del
          gruppo ",levels(group)[1], " con il costo ", mis.cost[1],"\n")
   }
   if(kind == "LDA"){
      fit <- lda(group ~ .,data = X, prior = as.numeric(prior))
   }else{
      fit <- qda(group ~ .,data = X, prior = as.numeric(prior))
   }
   print(fit)
   
   pred <- predict(fit)
   ######### PERFORMANCE DEL CLASSIFICATORE
   misc <- table(classe.vera=group, classe.allocata=pred$class)
   print("################# CLASSIFIER PERFORMANCES ###########################")
   print(misc)
   if(prior.sample){
      errori <- pred$class != group
      APER <- sum(errori) / n
   }else{
      APER <- 0
      # tieni conto nel calcolo dell'APER delle priorse non stimate dal campione
      # e quindi date in input
      for(i in 1:g){
         APER <- APER + (sum(misc[i,-i])/sum(misc[i,])) * 
            prior[i]
      }
   }
   
   if(print.result){
      cat("classifier's APER = ", APER, "\n")
   }
   
   ######### PREDICTION NEW OBSERVATION
   if(!is.null(X0.new)){
      new.pred <- predict(fit, X0.new)
      print("################### PREDICTION #####################")
      print(new.pred)
   }
   if(print.plot && p == 2){
      # extract the first two canonical component because we want to plot on a plane
      cc1 <- X[,1]
      cc2 <- X[,2]
      # ricordarsi che le medie sono nelle colonne
      cc.M <- Mg
      color.group <- group
      levels(color.group) <- rainbow(g)
      
      x11()
      # plotto le caratteristiche nel piani
      plot(cc1, cc2, main= kind, xlab=colnames(X)[1],
           ylab=colnames(X)[2], pch=20, col=as.character(color.group))
      legend(min(cc1), min(cc2)+2, legend = levels(group), fill= rainbow(g), cex=.7)
      
      for(i in 1:g){
         points(cc.M[1,i], cc.M[2,i], pch=4, col = rainbow(g)[i] , lwd=2, cex=1.5)
      }
      x.cc  <- seq(min(cc1),max(cc1),len=200)
      y.cc  <- seq(min(cc2),max(cc2),len=200)
      xy.cc <- expand.grid(cc1=x.cc, cc2=y.cc)
      names(xy.cc) <- colnames(X)
      z  <- predict(fit, xy.cc)$post    
      for(i in 1:g){
         zi <- z[,i] - z[,-i]
         contour(x.cc, y.cc, matrix(zi, 200), levels=0, drawlabels=F, add=T)
      }
   }
   
   return(list(pvalue.gauss = p.gauss,
               pvalue.bart = myifelse(p == 1 && kind == "LDA",bar.test,"no bartlett test"),
               model.fit = fit,
               APER = APER,
               sample.pred = pred,
               new.pred = myifelse(!is.null(X0.new),new.pred,"no new predictions"),
               prior = prior))
   
}# end function


myifelse <- function(condition,x,y){
   if(condition)
      return(x)
   else
      return(y)
   
}# end function
