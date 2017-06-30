# function "ConfidenceRegion" does the following things:
# 1) computes confidence regions for the true mean of a normal population
# 2) computes confidence regions for the true mean of any population (n large)
# 3) simultaneous confidence intervals (F)
# 4) Bonferroni simultaneous confidence intervals (t)
#
#  REMARK: 1)when the gaussian hp is used it is always check using mcshapiro test, the corresponding
#            p-value is printed.
#          2) call the function with the following different values of argument
#             to.do:
#             "CR" = confidence region
#             "T2" = F simulataneous conf intervals
#             "Bonf" = Bonferroni simultaneous conf intervals
# # INPUT:
#      X = dataframe 
#      large_n = if TRUE asymptotic analysis, otherwihe gaussian hp
#      S = sample covariance matrix
#      X_bar = sample mean vector
#      alpha = type 1 error
#      k = number of Bonferroni intrervals (maybe removable)
#      to.do = string for telling the function to do different things
#              "CR" --> confidence region @ level alpha
#              "T2" --> F simultaneous confidence intervals
#              "Bonf" --> Bonferroni confidence intervals
#      A = matrix whose rows are the vector where we want to compute
#          the simultaneous confidence intervals.
#          Example: A[i,] = e_i if you want the CI for the i-th component
#      plot = TRUE if you want the CI plot              
# OUTPUT:
#      IC = matrix with the simultaneous conf intervals
#  
ConfidenceRegion <- function(X = NULL,
                             large_n = F,
                             S = NULL,
                             X_bar = NULL,
                             alpha = 0.05,
                             k = NULL,
                             to.do = NULL,
                             A = NULL,
                             print.plot = F){
   
   
   # if not given as input, compute sample mean and sample covariance
   if(is.null(X_bar)){X_bar <- colMeans(X)}  
   if(is.null(S)){S <- cov(X)}  
   if(is.null(to.do)){stop("to.do not specified, i don't know what to do")}
   # sample cardinality        # population dimension
   n <- dim(X)[1];             p <- dim(X)[2]
   # range dataframe value, it's a global variable
   # if A, the matrix of the directions for sim conf intervals, is not specified; by default let's
   # compute the interval for the mean components
   if(is.null(A))  A <- diag(p)
   if(!large_n){# n is small, relay on Gaussian hp
      # check if the gaussian hp is supported by data (function defined below)
      check_gaussianity_test(X, alpha)
      # let's start inference analysis
      switch(to.do,
             CR = {# case  "confidence region for the mean"
                
                # plot the confidence region @ level 1-alpha for the mean
                plot_ellipse(X_bar, S, alpha = alpha, sample = T, n = n)
                # add the projection of the ellipse along the axis: in other
                # words compute the F simultaneous intervals and plot them
                IC <- simult_CI(X_bar, S, alpha = alpha, A, n, plot = F)
                segments(IC[1,1],0,IC[1,2],0,lty=1,lwd=2,col='red')
                segments(0,IC[2,1],0,IC[2,2],lty=1,lwd=2,col='red')
             },
             T2 = {# case "simulatneous confidence interval"
                # IC is a matrix dim(A)[1]x2, each row is a CI (function defined below)
                IC <- simult_CI(X_bar, S, alpha = alpha, A, n, print.plot = print.plot)
                # print the results
                print("################# T2 CONFIDENCE INTERVALS ###############")
                cat("Simultaneus confidence interval @ level alpha = ",alpha,"\n")
                print(IC)
                return(IC)
             },
             Bonf = {# case "Bonferroni confidence interval"
                
                # IC is a matrix dim(A)[1]x2, each row is a CI
                IC <- simult_CI(X_bar, S, alpha = alpha, A, n, Bonf = T, print.plot = print.plot)
                # print the results
                print("################# BONFERRONI CONFIDENCE INTERVALS ###############")
                cat("Bonferroni confidence interval @ level alpha = ",alpha,"\n")
                print(IC)
                return(IC)
             },
             {
                print("default on switch \n")
             }
      ) # end switch
   }
   else{# n is large, build an asyntotic confidence region for the mean
      
      if(p == 2 && print.plot){
         # plot the confidence region @ level 1-alpha for the mean
         plot_ellipse(X_bar, S, alpha = alpha, sample = T, n = n,large_n = T)
      }
      # costruzione intervalli confidenza chi-squared simultanei usando
      # l'hp di n grande e on quella di gaussianità
      k <- dim(A)[1] # numero CI
      CI <- data.frame(matrix(nrow = k,ncol = 2))
      cfr.chisq <- qchisq(1-alpha,p)
      for(i in 1:k){
         a <- A[i,]
         CI[i,] <- c(a%*%X_bar-sqrt(cfr.chisq)*sqrt(a%*%S%*%a/n),
                     a%*%X_bar+sqrt(cfr.chisq)*sqrt(a%*%S%*%a/n))
      }
      colnames(CI) <- c("LB","UB")
      print("############# CHI-SQUARED ASYNTOTIC CI ######################")
      cat("I seguenti CI sono per combinazione lineari a*mu ","\n")
      cat("intervalli di confidenza al livello alpha = ", alpha, "\n")
      print(CI)
      
   }
   
   
} # end function mean_inference

# compute simultaneous confidence interval for linear combnations of mu
simult_CI <- function(X_bar, S, alpha = 0.05, A, n, Bonf = FALSE,print.plot = FALSE){
   # INPUT:
   #       X_bar = sample mean
   #       S = sample covariance
   #       1-alpha = confidence level
   #       A = matrix whose ROWS ARE the direction for the linear comb
   #       n = sample cardinality
   #       Bonf = if TRUE Bonferroni interval otherwise F intervals
   #       plot = if TRUE plots the simultanueous intervals
   # OUTPUT:
   #       IC = matrix whose row are confidence intervals
   #  
   
   # population dimension     # number of simultaneous interval
   p <- dim(S)[2];            k <- dim(A)[1]
   # initialization: IC will contain the simulataneous intervals
   IC <- matrix(nrow = k, ncol = 3)
   # coefficient that appears in the IC formula, it depends on Bonf or F
   csi <- ifelse(Bonf,qt(1 - alpha/(2*k),n-1), sqrt(p*(n-1)/(n-p) * qf(1-alpha, p, n-p)))
   for(i in 1:k){
      # select the direction for the linear combination
      a <- A[i,]
      IC[i,] <- c(inf = t(a)%*% X_bar - sqrt(t(a)%*%S%*%a / n) * csi,
                  M = a%*%X_bar,
                  sup = t(a)%*% X_bar + sqrt(t(a)%*%S%*%a / n) * csi)
   }
   # plot the intervals
   if(print.plot){plot_intervals(IC, k)}  
   return(IC)
   
}# end simult_CI

# definition function "plot_intervals"
plot_intervals <- function(IC, k){
   #
   # INPUT :
   #        IC = confidence intervals matrix
   #        X_bar = vectof of sample means
   #        k = number of confidence intervals
   #   
   
   
   # open graphic device
   x11()
   # plot just the framework so as to use "segments" and "points"
   matplot(1:k, 1:k, pch='', ylim = range(IC), xlab='Variables',ylab='T2 for a component', 
           main='Simultaneous T2 conf. int. for the components')
   for(i in 1:k) segments(i,IC[i,1],i,IC[i,3],lwd=3,col=i)
   # add sample means
   points(1:k, IC[,2], pch=16, col=1:k)
   
} # end function plot_intervals

check_gaussianity_test <- function(X, alpha){
   
   # check for gaussianity
   mctest <- mcshapiro.test(X)
   if(mctest$pvalue < alpha)
      warning("WARNING: gaussian hp not supported by data, shapiro p-value = ",mctest$pvalue,"\n")
   else
      cat("You are assuming gaussianity, there's no evidence to say the contrary, shapiro p-value = ",mctest$pvalue,"\n")
}# end function "check_gaussianity_test"

