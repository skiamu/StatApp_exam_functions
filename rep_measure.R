# function for hypothesisi testing on a repeated measure dataset
# 
# INpUT: X = dataframe with observation in the rows and the columns are 
#            different experimental conditions on the same statistical unit
#        C = contrast matrix    
#        delta0 = null hp
#        
rep_measure <- function(X, C, delta0, alpha = 0.05, plot = T ){
   
   # dataframe dimensions
   n <- dim(X)[1];   q <- dim(X)[2]
   # check the contrast matrix dimension
   if(!all(dim(C) == c(q-1, q))){
      stop("dimensione contrast matrix non conformi")
   }
   # sample mean
   X_bar <- colMeans(X)
   # sample covariance
   S <- cov(X)
   # I suppose data drom a gaussian distribution, check this assumption
   check_gaussianity(X, alpha)
   
   ####### T2 test for C * mu #######
   # H0 : C * mu = delta.0
   # HA : C * mu != delta.0
   
   Md <- C %*% X_bar
   Sd <- C %*% S %*% t(C)
   if(det(Sd)<1e-4){stop("matrice Sd non invertibile")}
   inv.Sd <- solve(Sd)
   # test statistic
   T2 <- n * t(Md - delta0) %*% inv.Sd %*% (Md - delta0)
   # test p-value
   p.value <- 1 - pf(T2 * (n - q + 1) / ((n - 1) * (q - 1)),
                     q - 1,
                     n - q + 1)
   # print the results
   print("#######################################")
   cat("@ conf level alpha = ",alpha,ifelse(p.value<alpha,"we reject H0","we cannot reject H0"),"\n")
   cat("T2 test p-value = ", p.value,"\n")
   
   ####### test for the components of C * mu ##########
   k     <- q - 1   # number of increments (i.e., dim(C)[1])
   cfr.t <- qt(1-alpha/(2*k),n-1)
   # bonferroni CI @ level alpha for (C*mu)_j   j=1..k
   IC.BF <- cbind( Md - cfr.t*sqrt(diag(Sd)/n) , Md, 
                   Md + cfr.t*sqrt(diag(Sd)/n))
   
   cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) 
   # simultaneous F CI @ level alpha for (C*mu)_j   j=1..k
   IC.T2 <- cbind( Md - sqrt(cfr.fisher*diag(Sd)/n) , Md,
                   Md + sqrt(cfr.fisher*diag(Sd)/n))
   # print the results
   print("######################################")
   cat("Simultaneus confidence interval @ level alpha = ",alpha,"\n")
   print(IC.T2)
   print("######################################")
   cat("Bonferroni confidence interval @ level alpha = ",alpha,"\n")
   print(IC.BF)
   if(plot){
      library(plotrix) 
      x11()
      par(mfrow = c(1,2))
      plotCI(1:k, IC.T2[,2], ui = IC.T2[,3], li = IC.T2[,1])
      title(main = "simultaneous CI")
      for(i in 1:k){abline(h = delta0[i])}
      plotCI(1:k, IC.BF[,2], ui = IC.BF[,3], li = IC.BF[,1])
      title(main = "Bonferroni CI")
      for(i in 1:k){abline(h = delta0[i])}
      
   }
   
   
} # end function rep_measure







# definition of the function for checking gaussianity
check_gaussianity <- function(X, alpha){
   
   load("mcshapiro.test.RData")
   # check for gaussianity
   mctest <- mcshapiro.test(X)
   if(mctest$pvalue < alpha)
      warning("WARNING: gaussian hp not supported by data, shapiro p-value = ",mctest$pvalue,"\n")
   else
      cat("You are assuming gaussianity, there's no evidence to say the contrary, shapiro p-value = ",mctest$pvalue,"\n")
   
} # end function check_gaussianity

