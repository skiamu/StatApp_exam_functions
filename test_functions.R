# function for T2 Hotelling test of the mean vector.
# 
# This function performs the T2 test on the mean vector. The test is the 
# following:  H0 : mu = mu0
#             HA : mu != mu0
# If requested (C != NULL), it performs a simultaneous F - test on k linear
# combinations of the mean vector
# 
# INPUT:
#      X = dataframe 
#      large_n = if TRUE asymptotic analysis, otherwihe gaussian hp
#      mu0 = null hp for the T2 test 
#      S = sample covariance matrix
#      X_bar = sample mean
#      alpha = first type error
#      C = matrix whose rows are linear combianation coefficient
#         for the F-test
#      intervals = "Bonf" if bonferroni CI on vector components
#                  "T2" if T2 CI on vector component
# OUTPUT:
#      p.value = p-value of the T2 test or the asyntotic test (if n large)
# load function for checking gaussianity and functions for confidence region

T2.test <- function( X = NULL,
                     mu0 = NULL,
                     large_n = F,
                     S = NULL,
                     X_bar = NULL,
                     alpha = 0.05,
                     C = NULL,
                     intervals = "Bonf"
){
   
   # if not given as input, compute sample mean and sample covariance
   if(is.null(X_bar)) X_bar <- colMeans(X)
   if(is.null(S)) S <- cov(X)
   # sample cardinality    # population dimension
   n <- dim(X)[1];         p <- dim(X)[2]
   # max and min values within the dataframe, it's a global variable
   Range <<- range(X)
   # if small n we relay on gaussian hp, we need to check it
   if(!large_n) {
      p.gauss <- check_gaussianity_test(X,alpha)
   }
   if(!is.null(mu0)){
      # computed statistic
      x.T2 <- n * (X_bar - mu0)%*%solve(S)%*%(X_bar - mu0)
      cfr.qnt <- ifelse(large_n,qchisq(1-alpha,p),((n-1)*p)/(n-p)*qf(1-alpha,p,n-p))
      #compute the p_value
      p.value <- ifelse(large_n,1-pchisq(x.T2, p),1 - pf((n - p) / (p * (n - 1)) * x.T2, p, n-p))
      # print the result
      print("###################### T2 TEST ##################################")
      cat("@ conf level alpha = ",alpha,ifelse(p.value<alpha,"we reject H0","we cannot reject H0"),"\n")
      cat("T2 test p-value = ", p.value,"\n")
      cat("statistica calcolata T2.0 = ",x.T2,"\n")
      cat("quantile distribuzione (con eventauale coefficiente moltiplicativo) = ",cfr.qnt,"\n")
      # plot the rejection region if X is bidimensional
      if(p==2){
         # load file with function "plot_ellipse"
         source("/home/andrea/StatApp/StatApp_test/inference/plot_ellipse.R")
         # plot rejection region and give information about the ellipse
         plot_ellipse(mu0, S, alpha = alpha, sample = T, n = n,
                      large_n = ifelse(large_n,T,F),
                      title_plot = paste("Rejection region @ level ", alpha, sep = ""))
         # add the null hp, center of the ellipse
         points(mu0[1], mu0[2], pch = 19, cex = 1.5, lwd = 2, col ='blue')
         # add the sample mean, if this point is from too far from mu0 we
         # reject the null hp. "far" in the sense of mahalanobis distance
         points(X_bar[1], X_bar[2], pch = 4, cex = 1.5, lwd = 2, col ='red')
         legend("topright",legend = c("null hp","sample mean"),col = c("blue","red"),
                pch = c(19,4))
         # add observation
         points(X[,1],X[,2])
      }
      
      # compute simultaneous confidence intervals for the components
      # of the mean vector mu, if we've rejected H0 we want to see if
      # the sample mean of some components don't belong to its confidence 
      # interval
      IC <- ConfidenceRegion(X,large_n = F,alpha = alpha,to.do = intervals,print.plot = F)
      
   }
   # se C in input non è nulla sgnifica che voglio fare test simultaneo su
   # combinazioni lineari dele vettore delle medie
   if(!is.null(C)){
      # numero di combinazioni lineari da testare
      k <- dim(C)[1]
      X_bar.C <- C%*%X_bar
      S.C <- C%*%S%*%t(C)
      T2.0 <- n * t(X_bar.C)%*%solve(S.C)%*%(X_bar.C)
      cfr.fisher <- ((n-1)*k)/(n-k)*qf(1-alpha,k,n-k)
      
      # p.value del test
      p.value.sim <- 1-pf( T2.0 * (n-k)/((n-1)*k), k, n-k )
      print("################ TEST SIMULTANEO ##########################")
      cat("@ conf level alpha = ",alpha,ifelse(p.value.sim<alpha,"we reject H0","we cannot reject H0"),"\n")
      cat("p.value.sim = ", p.value.sim,"\n")
      cat("statistica calcolata T2.0 = ",T2.0,"\n")
      cat("quantile distribuzione (con eventauale coefficiente moltiplicativo) = ",cfr.fisher,"\n")
      # calcolo intervalli di confidenza simultanei per le combinazioni lineari
      cfr.fisher <- sqrt(((n-1)*k)/(n-k) * qf(1-alpha,k,n-k))
      IC.sim <- cbind(LB = X_bar.C - cfr.fisher * sqrt(diag(S.C)/n),
                      mean = X_bar.C,
                      UB = X_bar.C + cfr.fisher * sqrt(diag(S.C)/n))
      print("################ CI T2-SIMULTANEI #############################")
      cat("livello del test alpha = ", alpha, "\n")
      print("intervalli di confidenza simultanei (T2) per le combinazioni C")
      print(IC.sim)
   }
   
   return(list(p.value = myifelse(!is.null(mu0),p.value,"no test componenti"),
               IC = myifelse(!is.null(mu0),IC,"no IC for vector components"),
               IC_linear_comb = myifelse(!is.null(C),IC.sim,"no IC for linear comb of vector components"),
               p.gauss = p.gauss))
   
} # end t2.test function

# definition of the function for checking gaussianity
check_gaussianity_test <- function(X, alpha){
   
   # check for gaussianity
   mctest <- mcshapiro.test(X)
   if(mctest$pvalue < alpha)
      warning("WARNING: gaussian hp not supported by data, shapiro p-value = ",mctest$pvalue,"\n")
   else
      cat("You are assuming gaussianity, there's no evidence to say the contrary, shapiro p-value = ",mctest$pvalue,"\n")
   p.gauss <- mctest$pvalue
   return(p.gauss)
} # end function check_gaussianity_test

myifelse <- function(condition,x,y){
   if(condition)
      return(x)
   else
      return(y)
   
}# end function

