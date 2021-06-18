#### Application 1

require(coda)
require(R2jags)
require(rjags)
library(mcmc)
library(mcmcplots)
require(LaplacesDemon)
require(e1071)
require(gamlss)
require(nortest)

#Time between failures of secondary reactor pumps,
D3<-c(2.160, 0.746, 0.402, 0.954, 0.491, 6.560, 4.992, 0.347,
         0.150, 0.358, 0.101, 1.359, 3.465, 1.060, 0.614, 1.921,
         4.082, 0.199, 0.605, 0.273, 0.070, 0.062, 5.320)

hist(D3,freq=F)
lines(density(D3),lwd=3,col=2)

summary(D3)
sd(D3)
skewness(D3)
kurtosis(D3)
length(D3)
mean(D3)
sd(D3)


#fit1 <-gamlss(D3~1,family=FW1(mu.link = "identity",sigma.link = "identity"),
#              method = CG(100),
#              sigma.start = 0.0018)
#fit1 <-gamlss(D3~1,family=FW1(mu.link = "identity",sigma.link = "identity"),
#              method = mixed(5,100),
#              sigma.start = 0.0018)


fit1 <-gamlss(D3~1,family=FW1(mu.link = "identity",sigma.link = "identity"),n.cyc = 1200,
              sigma.start = 0.0018)

summary(fit1)

prof.dev(fit1, which="mu", min=0.1, max=1)
prof.dev(fit1, which="sigma", min=0.0005, max=0.5)




#Cramer Von Mises:
n<-length(D3)
x<-1:n # para calcular o CVM

(cvm<-((0.5/n)+1)*( (1/(12*n)) +sum((pFW1(sort(D3),
mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients)-(2*x-1)/(2*n))^2)))

#KS
#(ks<-max(x/n-pFW1(sort(D3),
 #            mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients),pFW1(sort(D3),
  #                                                                        mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients)-(x-1)/n))

ad.test(pFW1(D3,mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients))

cvm.test(pFW1(D3,mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients))

ks.test(D3,"pFW1",mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients)



#(ad<-((2.25/n^2)+(0.75/n)+1)*( -n - (1/n)*sum((2*x-1)*log(pWF1(sort(D3),
 # mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients)*(1-pWF1(sort(D3,decreasing = T),
  #                                                                mu=fit1$mu.coefficients, sigma=fit1$sigma.coefficients))))))




fit2 <-gamlss(D3~1,family=WEI3(mu.link = "identity",sigma.link = "identity"),n.cyc = 1200)

summary(fit2)

#Cramer Von Mises:
(cvm<-1/(12*n)+sum(((2*x-1)/(2*n)-pWEI3(sort(D3), mu=fit2$mu.coefficients, sigma=fit2$sigma.coefficients))^2))

#KS
ks.test(D3,"pWEI3",mu=fit2$mu.coefficients, sigma=fit2$sigma.coefficients)

#AD.tet
ad.test(pWEI3(D3,mu=fit2$mu.coefficients, sigma=fit2$sigma.coefficients))


cvm.test(pWEI3(D3,mu=fit2$mu.coefficients, sigma=fit2$sigma.coefficients))



dFW1 <- function(x,a,b){
  g <- exp((a*x - b/x)-exp(a*x - b/x))*(a+b/(x^2)) 
  return(g)
}

pFW1 <- function(x,a,b){
   G <- 1-exp(-exp(a*x - b/x)) 
    return(G)
}



hFW1 <- function(x,a,b){
  G <- 1-exp(-exp(a*x - b/x)) 
  g <- exp((a*x - b/x)-exp(a*x - b/x))*(a+b/(x^2)) 
  h <- g/(1-G)
  return(h)
}

(ekm <- survfit(Surv(D3)~1))
summary(ekm)


plot(ekm,conf.int=F,mark.time = T,lwd=2,col=c(1,2,3,4,5))
curve(1-pFW1(x, a =0.2071,b= 0.258),col="red",add=T,lwd=4,lty=2) # pdf



fit.sa2<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2]))) #minus the loglik  
  lower <- c(0.011,0.011) #may need some changes here
  upper <- c(1000,1000)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=7000,max.time=3))
  return(out[c("value","par","counts")])
}
dFW <- function(x,a,b){
  theta=1
  alpha=1
  g <- exp((a*x - b/x)-exp(a*x - b/x))*(a+b/(x^2)) 
  G <- 1-exp(-exp(a*x - b/x)) 
  f <- (alpha*theta*g*(G^(alpha*theta-1)) *(1-G^theta)^(alpha-1)) /(((G^(alpha*theta)) + (1-G^theta)^(alpha))^2)
  F <- (G^(alpha*theta))/(((G^(alpha*theta))+(1-G^theta)^alpha))
  return(f)
}
fit.sa2(D3,dFW)
#Densidade
dFW <- function(par,x){
  a=par[1]
  b=par[2]
  g <- exp((a*x - b/x)-exp(a*x - b/x))*(a+b/(x^2)) 
  G <- 1-exp(-exp(a*x - b/x)) 
    return(g)
}
#Acumulada
pFW <- function(par,x){
  a=par[1]
  b=par[2]
  g <- exp((a*x - b/x)-exp(a*x - b/x))*(a+b/(x^2)) 
  G <- 1-exp(-exp(a*x - b/x)) 
    return(G)
}

gFW <- goodness.fit(pdf=dFW, cdf=pFW,
                    starts = c(0.3456835,2.3478205), data = D3,method="BFGS", domain=c(0,Inf),mle=NULL);gFW




X11()
hist(D3,freq=F,15,main="",col="gray85",xlab="x",ylab="f(x)",xlim=c(0.001,max(D3)+0.5),ylim=c(0,1.3),cex.axis=1.35,cex.lab=1.35)
curve(dFW1(x, mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients),col="gray35",add=T,lwd=4) # pdf
curve(dWEI3(x, mu=fit2$mu.coefficients ,sigma=fit2$sigma.coefficients),col="gray55",add=T,lwd=4,lty=2) # pdf
curve(dFW1(x, a =0.2071,b= 0.258),col="red",add=T,lwd=4,lty=2) # pdf

legend("topright",lty=c(1,2),lwd=c(3,3),c("FW1","WEI3"),col=c("gray25","gray55"),bty="n",cex=1.35)
box()

x11()
plot(ecdf(D3),main="",col="gray1",xlab="x",ylab="F(x)",xlim=c(min(D3),max(D3)+0.5),cex.axis=1.35,cex.lab=1.35)
curve(pFW1(x, mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients),col="gray25",add=T,lwd=4) # pdf
curve(pWEI3(x, mu=fit2$mu.coefficients ,sigma=fit2$sigma.coefficients),col="gray55",add=T,lwd=4,lty=2) # pdf
curve(pFW(x, a= 0.0207,b= 2.5875),col="red",add=T,lwd=4,lty=2) # pdf

legend("topleft",lty=c(1,1,2),lwd=c(3,3,3),c("Empirical distribution","FW1","WEI3"),col=c("gray25","gray55"),bty="n",cex=1.35)


x11()
curve(hFW1(x, mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients),col="gray25",add=F,lwd=4,xlim=c(0.001,max(D3)+0.5),cex.axis=1.35,cex.lab=1.35,ylab="h(x)",xlab="x") # pdf
curve(dWEI3(x,mu=fit2$mu.coefficients ,sigma=fit2$sigma.coefficients)/(1-pWEI3(x,mu=fit2$mu.coefficients ,sigma=fit2$sigma.coefficients)),add=T,lwd=4,lty=2,col="gray55")
curve(hFW(x, a= 0.0207,b= 2.5875),col="red",add=T,lwd=4,lty=2) # pdf


legend("topright",lty=c(1,2),lwd=c(3,3),c("FW1","WEI3"),col=c("gray25","gray55"),bty="n",cex=1.35)

curve(hFW(x, a= 0.0207,b= 2.5875),col="red",add=T,lwd=4,lty=2,xlim=c(0,60)) # pdf


ekm$surv[23]<-0.001

y <- D3

y <- log(-log(ekm$surv))
X1 <- D3
X2 <- -1/D3

X = cbind(X1,X2)

solve(t(X)%*%X)%*%t(X)%*%y

lm(y~X1+X2-1)


hist(log(-log(ekm$surv)),freq=F)
curve(dFW(x, a= 0.0207,b= 2.5875),col="red",add=T,lwd=4,lty=2,xlim=c(0,60)) # pdf

curve(dFW(x, a= 0.0207,b= 2.5875),col="red",add=F,lwd=4,lty=2,xlim=c(0,60)) # pdf


#require(bshazard)
#fit.bshazard11 <- bshazard(Surv(D2, rep(1,length(D2))) ~ 1)
#x11()
#plot(fit.bshazard11$time, fit.bshazard11$hazard,lty=2,lwd=4,xlab='x', ylab='h(x)',col=1,cex.lab=1.2,cex.axis=1.2,type = "s")

f <- function(x,mu,sigma){
  fw1 <- dWF1(x,mu=mu,sigma=sigma)
  return(fw1)
}

f1 <- function(x,a,b){
  fg <- dgamma(x,a,b)
  return(fg)
}
integrate(f1,lower = 0.0001,upper = Inf,a=90,b=1)

integrate(f,lower = 0.0001,upper = Inf,mu=20,sigma=0.5)



###################################
##An?lise Bayesiana
cat("
data{
 C <- 100
 for (i in 1:n){ zeros[i] <- 0}
}
model{
for (i in 1:n){
zeros[i] ~ dpois(phi[i])
phi[i] <- - l[i] + C
l[i] <-  log((exp((sigma*y[i] - ((((2*sigma*mu-(log(log(2))))^2)
-(log(log(2)))^2)/(4*sigma))/y[i])-exp(sigma*y[i] - 
((((2*sigma*mu-(log(log(2))))^2)-
(log(log(2)))^2)/(4*sigma))/y[i])))*
(sigma+((((2*sigma*mu-(log(log(2))))^2)-(log(log(2)))^2)/
(4*sigma))/(y[i]^2)))
}
 #prior's distribuition
mu ~ dgamma(0.001,0.001)
sigma ~ dgamma(0.001,0.001)
}", file="WF_mu3.txt")


y <- D3
n <- length(y)

data3 <- list(n=n,  y=y)

inits1 <-list("mu"=median(D3), "sigma"=min(D3)) 
#inits2 <- list("mu"=mean(D3), "sigma"=0.1) 
inits2 <- list("mu"=3, "sigma"=1.2) 
#bayes.mod.inits <- list(inits1, inits2, inits3)
bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(943)
Modelo.FW <- jags(data=data3, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                   n.chains=2, n.iter=60000,n.burnin = 20000,
                   n.thin =10,
                   model.file="WF_mu3.txt")

#Outputs
print(Modelo.FW)


#valor de DIC  c=100
DIC <- Modelo.FW$BUGSoutput$DIC
DIC
#valor de DIC corrigido 
DIC.c <- (DIC-2*23*100)
DIC.c

#Deviance m?dio
d.bar <- mean(Modelo.FW$BUGSoutput$sims.list$deviance)
d.bar
#Deviance m?dio corrigido
d.bar.c <- (d.bar-2*23*100)
d.bar.c

#pD
Modelo.FW$BUGSoutput$pD

DIC-d.bar
DIC.c - d.bar.c

#Estatistica EAIC
p <- 2 #n?mero de par?metros do modelo

EAIC.c <- d.bar.c +2*p
EAIC.c 

EBIC.c <- d.bar.c +p*log(n)
EBIC.c 


#m?dia em rela??o a cadeia de mu e sigma
mean(Modelo.FW$BUGSoutput$sims.list$mu)
mean(Modelo.FW$BUGSoutput$sims.list$sigma)


hist(Modelo.FW$BUGSoutput$sims.list$mu,freq=F,xlim=c(0,1))
lines(density(Modelo.FW$BUGSoutput$sims.list$mu))
curve(dgamma(x,0.01,0.001),add=T,lwd=3)


#Intervalos para os par?metros
plot(Modelo.FW)

x11()
traceplot(Modelo.FW)

Modelo.FW.mcmc <- as.mcmc(Modelo.FW)
summary(Modelo.FW.mcmc)

gelman.plot(Modelo.FW.mcmc)

x11()
plot(Modelo.FW.mcmc)

require(superdiag)
mcmcplot(Modelo.FW.mcmc)
superdiag(Modelo.FW.mcmc, burnin = 1000)


traceplot(Modelo.FW.mcmc)
x11()
denplot(Modelo.FW.mcmc)
x11()
autocorr.plot(Modelo.FW.mcmc)

x11()
gelman.plot( Modelo.FW.mcmc )
geweke.diag( Modelo.FW.mcmc )

#################################################
#Weibullllll
cat("
data{
 C <- 100
 for (i in 1:n){ zeros[i] <- 0}
}
model{
for (i in 1:n){
zeros[i] ~ dpois(phi[i])
phi[i] <- - l[i] + C
l[i] <-  log((sigma/(mu/(exp(loggam((1/sigma)+1)))))*((y[i]/(mu/(exp(loggam((1/sigma)+1)))))^(sigma-1))*
exp(-(y[i]/(mu/(exp(loggam((1/sigma)+1)))))^sigma))
}
 #prior's distribuition
mu ~ dgamma(0.001,0.001)
sigma ~ dgamma(0.001,0.001)
}", file="WEI3.txt")

y <- D3
n <- length(y)
data3 <- list(n=n,  y=y)

inits1 <-list("mu"=mean(D3), "sigma"=min(D3)) 
#inits2 <- list("mu"=mean(D3), "sigma"=0.1) 
inits2 <- list("mu"=3, "sigma"=1.2) 
#bayes.mod.inits <- list(inits1, inits2, inits3)
bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(563)
Modelo.WEI <- jags(data=data3, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                  n.chains=2, n.iter=60000,n.burnin = 20000,
                  n.thin =10,
                  model.file="WEI3.txt")

#Outputs
print(Modelo.WEI)


#valor de DIC  c=100
DIC1 <- Modelo.WEI$BUGSoutput$DIC
DIC1
#valor de DIC corrigido menos c=100
DIC.c1 <- (DIC1-2*23*100)
DIC.c1

#Deviance m?dio
d.bar1 <- mean(Modelo.WEI$BUGSoutput$sims.list$deviance)
d.bar1
#Deviance m?dio corrigido
d.bar.c1 <- (d.bar1-2*23*100)
d.bar.c1

#pD
Modelo.WEI$BUGSoutput$pD

DIC1-d.bar1
DIC.c1 - d.bar.c1

#Estatistica EAIC
p <- 2 #n?mero de par?metros do modelo

EAIC.c1 <- d.bar.c1 +2*p
EAIC.c1 

EBIC.c1 <- d.bar.c1 +p*log(n)
EBIC.c1 


#m?dia em rela??o a cadeia de mu e sigma
mean(Modelo.WEI$BUGSoutput$sims.list$mu)
mean(Modelo.WEI$BUGSoutput$sims.list$sigma)


hist(Modelo.WEI$BUGSoutput$sims.list$mu,freq=F,xlim=c(0,3))
lines(density(Modelo.WEI$BUGSoutput$sims.list$mu))
curve(dgamma(x,0.001,0.001),add=T,lwd=3)


#Intervalos para os par?metros
plot(Modelo.WEI)

x11()
traceplot(Modelo.WEI)

Modelo.WEI.mcmc <- as.mcmc(Modelo.WEI)
summary(Modelo.WEI.mcmc)

gelman.plot(Modelo.WEI.mcmc)

x11()
plot(Modelo.WEI.mcmc)

#require(superdiag)
mcmcplot(Modelo.WEI.mcmc)
superdiag(Modelo.WEI.mcmc, burnin = 1000)


traceplot(Modelo.WEI.mcmc)
x11()
denplot(Modelo.WEI.mcmc)
x11()
autocorr.plot(Modelo.WEI.mcmc)

x11()
gelman.plot( Modelo.WEI.mcmc )
geweke.diag( Modelo.WEI.mcmc )
