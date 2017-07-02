# script per il classificatore 1D
# sia X il vettore con le caratteristiche unidimensionali e group il vettore
# dei fattori

X <- as.numeric(X$lunghezza)
idx1 <- which(group == levels(group)[1])
idx2 <- which(group == levels(group)[2])
X1 <- X[idx1]
X2 <- X[idx2]
# verifico assunzioni
p.gauss <- c(shapiro.test(X1)$p.value,shapiro.test(X2)$p.value)
p.sigma <- var.test(X1,X2)$p.value
print("##################### CHECK HP ##################################")
cat("shapiro test p.value = ", p.gauss, "\n")
cat("chi.squared test p.value = ", p.sigma, "\n")

ng <- table(group)
n <- sum(ng)
M1 <- mean(X1)
M2 <- mean(X2)
S <- (var(X1) * (ng[1] - 1) + var(X2) * (ng[2] - 1)) / (ng[1] + ng[2] - 2)
SD <- sqrt(S)

####### DEFINISCO PRIOR E MIS.COST:
# nella prima posizione metto c(2|1) il costo che va a moltiplicare la prima distribuzione
# nella seconda posizione metto c(1|2) il costo che va a moltiplicare la seconda distribuzione
mis.cost <- c(1,1)
# nella prima posizione metto la prior che va a moltiplicare la prima distribuzione
# e nella seconda posizione quella che va a moltiplicare la seconda.
prior <- c(0.75,0.25)

# definisco l'expected misclassification cost: resta da decidere come integrare:
# devo integrare da x a +inf quella con la media minore (quindi 1-pnorm) mentre
# devon integrare da -Inf a x quella con la media maggiore
ECM <- function(x){
   prior[2]*mis.cost[2]*(1- pnorm(x, M2, SD)) + 
   (pnorm(x, M1, SD) )* prior[1] * mis.cost[1]
}
R   <- optimize(f=ECM, lower=min(X), upper=max(X))$minimum
AER <- optimize(f=ECM, lower=min(X), upper=max(X))$objective

# check
AER == prior[1]*(pnorm(R, M1, SD)) + prior[2]*(1-pnorm(R, M2, SD))



X <- data.frame(feature = X)
# fitto LDA
fit <- lda(group ~ feature, prior= prior,data = X)
# griglia dove plottare
x <- data.frame(feature = seq(min(X), max(X), by = 0.05))

LDA <- predict(fit, x)$posterior[,2] # classe 2

x11()
par(mfrow=c(2,1))
plot(x[,1], prior[2]*(dnorm(x[,1], M2, SD)), type='l', col='red',
     lty=1, xlab='x', ylab='density * prior', ylim=c(0,0.6))
lines(x[,1], prior[1]*(dnorm(x[,1], M1, SD)), type='l', col='blue',
      lty=1, xlab='x', ylab='density * prior')
abline(v=R, lty=2)
points(X1, rep(0, length(X1)), pch=3, col='blue')
points(X2, rep(0, length(X2)), pch=3, col='red')
legend(7,.5,legend=c('X1','X2'),fill=c('blue','red'),cex=.7)

plot(x[,1], LDA, type='l', col='red', lty=1, xlab='x',
     ylab='estimated posterior')# rosso = iberica
lines(x[,1], 1 - LDA, type='l', col='blue', lty=1, xlab='x',
      ylab='estimated posterior')# blu = atlantica
abline(h = 0.5)
abline(v=R, lty=2)

points(X1, rep(0, length(X1)), pch=3, col='blue')
points(X2, rep(0, length(X2)), pch=3, col='red')

dev.off()

# Compute the APER
prior <- c(0.75,0.25)
G <- 2
misc <- table(classe.vera=group, classe.allocata=predict(fit)$class)
misc

APER <- 0
for(g in 1:G)
   APER <- APER + sum(misc[g,-g])/sum(sum(misc[g,])) * prior[g]  
APER
AER



}