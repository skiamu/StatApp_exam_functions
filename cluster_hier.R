# funzione per il cluster gerarchico agglomerativo
# 
# INPUT : 
#        X = dataframe con le osservazioni da clusterizzare nelle colonne
#        distance = vettore con le distanze da usare nel calcolo della matrice
#                   delle distanze [ vettore di stringhe]
#        linkage = vettore con il metodo da usare per unire due cluster durante
#                  l'algoritmo [vettore di stringhe]           
#        k = numero di cluster che si vuole formare, prima si chiama la funzione
#            con k = NULL, si decide quanti cluster formare e poi si richiama
#            la funzione decidendo k
#            
# OUTPUT :
#        X.distance = lista con gli output della funzione dist per ogni 
#                     coppia distance/linkage
#        X.clust =  lista con gli output della funzione hclust  per ogni combinazione della
#                  coppia distance/linkage
#        coph.coeff = coefficiente cofanetico
#        label = lista con gli outputr della funzione cutree per ogni distance/linkage        
# 
# 
# OSS :




Clustering <- function(X,
                       distance,
                       linkage,
                       print.plot = T,
                       k = NULL){
   # numero di distanze
   n_d <- length(distance)
   # numero di linkage
   n_l <- length(linkage)
   # plotta uno scatter plot per vedere se si individua qualche gruppo nelle 
   # marginali
   # if(print.plot){
   #    x11()
   #    plot(X)
   # }
   # lista che contiene le matrici delle distanze per ogni distanza in input
   X.distance <- vector("list",n_d)
   names(X.distance) <- distance
   for(i in 1:n_d){
      X.distance[[i]] <- dist(X, method = distance[i])
   }
   # lista che conterrà i clustering per ogni distanza in input per ogni linkage
   X.clust <- vector("list",n_d * n_l)
   nomi <- vector("character", n_d * n_l)
   coph.coeff <- vector("double", n_d * n_l)
   counter <- 0
   for(i in 1:n_d){
      for(j in 1:n_l){
         counter <- counter + 1
         nomi[counter] <- paste(distance[i],linkage[j],sep = ".")
         X.clust[[counter]] <- hclust(X.distance[[i]], method = linkage[j])
         coph <- cophenetic(X.clust[[counter]])
         coph.coeff[counter] <- cor(X.distance[[i]],coph)
      }
   }
   # dai nomi alla lista coi cluster
   names(X.clust) <- names(coph.coeff) <- nomi
   print("#################################################################")
   print(coph.coeff)
   # print the dendrogram
   if(print.plot){
      counter <- 0
      for(i in 1:n_d){
         x11()
         par(mfrow=c(1,n_l))
         for(j in 1:n_l){
            counter <- counter + 1
            plot(X.clust[[counter]], main=nomi[counter], hang=-0.1, xlab='', 
                 labels=F, cex=0.6, sub='')
         }
      }
   }
   X.label <- vector("list",n_d * n_l)
   names(X.label) <- nomi
   if(!is.null(k)){
      counter <- 0
      for(i in 1:n_d){
         for(j in 1:n_l){
            counter <- counter + 1
            X.label[[counter]] <- cutree(X.clust[[counter]],k)
            if(print.plot){
               x11()
               plot(X,col = X.label[[counter]]+1,main = nomi[counter])
            }
         }
      }
   }   
   return(list(distanze = X.distance,
               cluster = X.clust,
               coph.coeff = coph.coeff,
               label = X.label))
   
}# end function
