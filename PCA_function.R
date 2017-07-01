# source: https://www.ime.usp.br/~pavan/pdf/MAE0330-PCA-R-2013
#
# R function to perform principal component analysis using either spectral or singular 
# value decomposition performed by the functions eigen and svd, respectively. 
# Output from the function consists of an object containing the eigenvalues, 
# eigenvectors, scores, summary of output, and the standard deviations of the new components. 
# The default settings are to use the correlation matrix (i.e., scaled=T)
# for spectral decomposition, remove all observations with missing information, 
# and to print a summary of the results. To perform singular value decomposition 
# the argument for method should be changed to "svd". By specifying the argument graph=T, 
# the function will provide two graphs: (1) a screeplot and (2) a biplot of the first two components. 
# If there are missing values within the dataset, the user has two options: 
# (1) rm.na=T (default) removes all observations with missing information or 
# (2) rm.na=F replaces missing values with the mean of the variable.
# 
# call the function in this way: 
#        pc_data1<-PC(X,method="eigen",scaled=T,graph=F,rm.na=T,print.results=T)
PC <- function(X,
               method="eigen",
               scaled=T,
               graph=F,
               rm.na=T,
               print.results=T
){
   if (any(is.na(X))){
      tmp<-X
      if(rm.na==T){X<-na.omit(data.frame(X));X<-as.matrix(X)}
      else{X[is.na(X)] = matrix(apply(X, 2, mean, na.rm = TRUE),
                                ncol = ncol(X), nrow = nrow(X), byrow = TRUE)[is.na(X)]}}
   else{tmp<-X}
   
   if(method=="eigen"){
      if(scaled==1){X1=cor(X);X2=scale(X)}
      else{X1=cov(X);X2=scale(X,scale=F)}
      total.var<-sum(diag(cov(X2)))
      values<-eigen(X1)$values;vectors<-eigen(X1)$vectors;sdev=sqrt(values)}
   
   if(method=="svd"){
      if(sum(scaled,center)>1){X2<-scale(X)}else{
         if(scaled==1){X2=scale(X,center=F)}else{
            if(center==1){X2=scale(X,scale=F)}else{X2=X}}}
      total.var<-sum(diag(cov(X2)))
      var<-nrow(X2)-1
      vectors<-svd(X2)$v;sdev=svd(X2)$d/sqrt(var);values<-sdev*sdev}
   
   prop.var<-rep(NA,ncol(X));cum.var<-rep(NA,ncol(X));scores<-X2%*%vectors
   namex<-as.character(1:ncol(X));scorenames<-rep(NA,ncol(X))
   for(i in 1:(ncol(X))){
      scorenames[i]<-do.call(paste,c("PC",as.list(namex[i]),sep=""))}
   colnames(scores)<-scorenames
   rownames(vectors)<-colnames(X);colnames(vectors)<-scorenames
   for(i in 1:ncol(X)){prop.var[i]<-var(scores[,i])/total.var}
   for(i in 1:ncol(X)){cum.var[i]<-sum(prop.var[1:i])}
   importance<-t(matrix(c(sdev,prop.var,cum.var),ncol=3))
   importance<-as.table(importance)
   colnames(importance)<-scorenames
   rownames(importance)<-c("Standard Deviation","Proportion of Variance","CumulativeProportion")
   z<-list(values=values,vectors=vectors,scores=scores,importance=importance
           ,sdev=sdev)
   if(graph==1){
      x11()
      biplot(scores[,1:2],vectors[,1:2], main="Biplot of Data",xlab=do.call
             (paste,c("PC1 (",as.list(round(z$importance[2,1]*100,2)),"%)",sep=""))
             ,ylab=do.call(paste,c("PC2(",as.list(round(z$importance[2,2]*100,2)),"%)",sep="")), cex=0.7)
      x11()
      screeplot(z,type='l',main='Screeplot of Components')
      abline(1,0,col='red',lty=2)
      # add the barplots for easing interpretation
      x11()
      p <- dim(X2)[2]
      par(mar = c(1,4,0,2), mfrow = c(p,1))
      for(i in 1:p) barplot(vectors[,i], ylim = c(-1, 1))
      # Explained variance with barplots
      x11()
      layout(matrix(c(2,3,1,3),2,byrow=T))
      plot(princomp(X2), las=2, main='Principal components', ylim=c(0,ifelse(scaled,6,max(princomp(X2)$sdev^2))))
      barplot(sapply(data.frame(X2),sd)^2, las=2, main='Original Variables', ylim=c(0,ifelse(scaled,6,max(princomp(X2)$sdev^2))), ylab='Variances')
      plot(cumsum(sdev^2)/sum(sdev^2), type='b', axes=F, xlab='number of components', 
           ylab='contribution to the total variance', ylim=c(0,1))
      abline(h=1, col='blue')
      abline(h=0.8, lty=2, col='blue')
      box()
      axis(2,at=0:10/10,labels=0:10/10)
      axis(1,at=1:ncol(X2),labels=1:ncol(X2),las=2)
      # boxplot original variables vs principal components
      x11()
      layout(matrix(c(1,2),2))
      boxplot(X2, las=2, col='red', main='Original variables')
      boxplot(scores, las=2, col='red', main='Principal components')
      }
   if(print.results==T){
      if(method=="eigen"){print("PCA Analysis Using Spectral Decomposition")}
      if(method=="svd"){print("PCA Analyis Using Singular Value Decomposition")}
      if (any(is.na(tmp))){
         if(rm.na==T){print("Warning:One or more rows of data were omitted from analysis")}
         if(rm.na==F){print("Warning: Mean of the variable was used for Missing values")}}
      print(importance)}
   z<-list(values=values,vectors=vectors,scores=scores,importance=importance
           ,sdev=sdev)
   
}#End Function