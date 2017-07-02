




run_kmeans <- function(X,how.choose.k = F){
   
   # plotta i dati per avere un'idea su che k scegliere per l'algoritmo
   # k-means
   x11()
   plot(X)
   # leggi in input il numero di cluster
   k <- readinteger()
   # runna k-means
   result.k <- kmeans(X, centers=k)
   
   x11()
   plot(X, col = result.k$cluster+1)
   
   if(how.choose.k){
      b <- w <- NULL
      for(k in 1:10){
         result <- kmeans(X, k)
         w <- c(w, sum(result$wit))
         b <- c(b, result$bet)
      }
      
      x11()
      matplot(1:10, b/(w+b), pch='', xlab='clusters', ylab='between/tot', main='Choice of k', ylim=c(0,1))
      lines(1:10, b/(w+b), type='b', lwd=2)
      
      x11()
      matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
      lines(1:10, w/(w+b), type='b', lwd=2)
   }

   return(result.k)
}

readinteger <- function()
{ 
   n <- readline(prompt="Enter an integer: ")
   if(!grepl("^[0-9]+$",n))
   {
      return(readinteger())
   }
   
   return(as.integer(n))
}
