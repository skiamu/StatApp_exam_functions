# this function draws an ellipse:
#1) sample = TRUE --> it means that matrix Sigma is the sample variance. In this case
#                     the elliptic region is a confidence region for the mean vector.
#                     Depending on small n or large n the functions computes the
#                     confidence regiob @ level alpha
#2) sample = FALSE --> it means that matrix Sigma is the variance matrix of
#                      a normal bivariate population. In this case the elliptic
#                      region is a probability reagion @ level 1-alpha 
#
#
# REMARK: the probability elliptic region is centered in mu (true mean)
#         while the confidence elliptic region for the means is a stochastic
#         region in R^2 that (1-alpha)% of the times contains the true mean and it's
#         centered on the sample mean x_bar.
plot_ellipse <- function (mu,
                          Sigma,
                          alpha = 0.05,
                          sample = T,
                          n = NULL,
                          large_n = F,
                          title_plot = NULL,
                          print.result = T){
   # INPUT:
   #      mu = mean vector (2-D)
   #      Sigma = 2x2 positive def matrix
   #      alpha = 
   #      sample = T if Sigma is the sample varianze  F if it's the population's
   #               variance
   #      n = sample cardinality (numerosità campionaria)
   #      large_n = if TRUE asymptotic region
   #      title_plot = title of ellipse plot (rejection region, confidence region ecc) 
   # OUTPUT: none
   
   library(car) 
   # Sigma must be 2x2
   if(all(dim(Sigma) != c(2,2))){
      stop("dimensioni matrice errate per plot 2-D")
   }
   # applico teorema spettarle: if sample = TRUE I havt to apply the spectral
   # theorem to the matrix Sigma/n which is the covariance matrix of the vector
   # X_bar
   w <- eigen(myifelse(sample,Sigma/n,Sigma)) 
   # estraggo autovalori
   lambda1 <- w$values[1]
   lambda2 <- w$values[2]
   # estraggo autovettori
   e1 <- w$vectors[,1]
   e2 <- w$vectors[,2]
   # dimensione popolazione
   p <- 2 
   # quantile distribuzione: se campionario e non large_n allora è una F,
   # se campionario e large_n oppure non campionario è una chi-quadro
   qnt <- myifelse(sample && !large_n, qf(1-alpha,p,n - p), qchisq(1-alpha,2))
   # raggio ellisse
   r <- myifelse(sample && !large_n,sqrt((n-1)*p/(n-p))* 
                    sqrt(qnt),sqrt(qnt))
   S <- myifelse(sample,Sigma/n,Sigma)
   
   # stampa caratteristiche ellisse
   max.semiax <- sqrt(lambda1) * r
   min.semiax <- sqrt(lambda2) * r 
   if(print.result){
      print("--------------------")
      cat("(lambda1, lambda2) = ",c(lambda1,lambda2),"\n",sep = "  ")
      print("--------------------")
      print("autovettori = ")
      print(w$vectors)
      print("--------------------")
      cat("semiasse maggiore = ", max.semiax,"\n")
      cat("semiasse minore = ",min.semiax,"\n")
      cat("raggio = ",r,"\n")
   }
   # ylim <- c(mu[2] - 2 * sqrt(S[2,2]), mu[2] + 2 * sqrt(S[2,2]))
   # xlim <- c(mu[1] - 2 * sqrt(S[1,1]), mu[1] + 2 * sqrt(S[1,1]))
   ylim <- c(mu[2] - max.semiax,mu[2] + max.semiax)
   xlim <- c(mu[1] - max.semiax, mu[1] + max.semiax)
   # plotting
   x11()
   if(is.null(title_plot)){
      title_plot <- ifelse(sample,ifelse(large_n,paste("Asyntotic Confidence Region for mu @ level alpha ",alpha,sep = ""),paste("Confidence Region  for mu @ level alpha ",alpha,sep = "")),paste("Probability Region @ level alpha ",alpha,sep = ""))                                                                     
   }
   plot(mu[1],mu[2],
        ylim= ylim,
        xlim = xlim,
        col='blue',pch=19,xlab='X.1',ylab='X.2',asp=1,
        main = title_plot)
   ellipse(center=mu, shape= S, radius=r, col = 'blue')
   
   m = e1[2]/e1[1]
   abline(a=mu[2] - m * mu[1], b = m, lty=2, col='grey')
   abline(a=mu[2] + mu[1]/m,b=-1/m,col='grey',lty=2)
   abline(v=0,col='grey',lty=2)
   abline(h=0,col='grey',lty=2)
}# end function

# implementazione funzione myifelse
myifelse <- function(condition,x,y){
   if(condition)
      return(x)
   else
      return(y)
}


