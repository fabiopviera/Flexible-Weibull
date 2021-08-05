#### Application 2
require(coda)
require(R2jags)
require(rjags)
library(mcmc)
library(mcmcplots)
require(LaplacesDemon)
require(e1071)
require(gamlss)
require(gamlss.cens)
require(survival)

source("R/FW1_GAMLSS.R")


t    <- c(7,34,42,63,64,74,83,84,91,108,112,129,133,133,139,140,140,146,149,154,157,160,160,165,173,176,
          185,218,225,241,248,273,277,279,297,319,405,417,420,440,523,523,583,594,1101,1116,1146,1226,1349,1412,1417)
cens <- c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          0,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,0,1)


#Descriptive statistics by K-M
(ekm <- survfit(Surv(t,cens,type ="right" )~1))
summary(ekm)

plot(ekm,conf.int=T,mark.time = T,col="black",lwd=2)




# FW1 Model
fit1 <-gamlss(Surv(t, cens)~1,family=cens(FW1)(mu.link = "identity",sigma.link = "identity"),n.cyc = 1200,
              sigma.start = 0.0018)

summary(fit1, type="qr")

# WEI3 Model
fit2 <-gamlss(Surv(t, cens)~1,family=cens(WEI3)(sigma.link = "identity"),n.cyc = 1200)
summary(fit2, type="qr")





x11()
plot(ekm,conf.int=F,mark.time = T,col="gray45",lwd=3,cex.axis=1.35,cex.lab=1.35,xlab="Time (in days)",ylab="S(x)")
curve(1-pFW1(x, mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients),col="gray25",add=T,lwd=4) # pdf
curve(1-pWEI3(x, mu=exp(fit2$mu.coefficients) ,sigma=fit2$sigma.coefficients),col="gray55",add=T,lwd=4,lty=2) # pdf
legend("topright",lty=c(1,1,2),lwd=c(3,3,3),c("Empirical survival","FW1","WEI3"),col=c("gray45","gray25","gray55"),bty="n",cex=1.35)


x11()
curve(hFW1(x, mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients),col="gray25",add=F,lwd=4,xlim=c(0.001,1450),cex.axis=1.35,cex.lab=1.35,ylab="h(x)",xlab="Time (in days)") # pdf
curve(dWEI3(x,mu=exp(fit2$mu.coefficients) ,sigma=fit2$sigma.coefficients)/(1-pWEI3(x,mu=exp(fit2$mu.coefficients) ,sigma=fit2$sigma.coefficients)),add=T,lwd=4,lty=2,col="gray55")
legend("topright",lty=c(1,2),lwd=c(3,3),c("FW1","WEI3"),col=c("gray25","gray55"),bty="n",cex=1.35)


####################################################################
# Bayesian FW1

cat("
data{
 C <- 100
 for (i in 1:n){ zeros[i] <- 0}
}
model{
for (i in 1:n){
zeros[i] ~ dpois(phi[i])
phi[i] <- - l[i] + C
l[i] <-  cens[i]*log((exp((sigma*y[i] - ((((2*sigma*mu-(log(log(2))))^2)
-(log(log(2)))^2)/(4*sigma))/y[i])-exp(sigma*y[i] - 
((((2*sigma*mu-(log(log(2))))^2)-
(log(log(2)))^2)/(4*sigma))/y[i])))*
(sigma+((((2*sigma*mu-(log(log(2))))^2)-(log(log(2)))^2)/
(4*sigma))/(y[i]^2)))    +     (1-cens[i])*log((exp(-exp(sigma*y[i] - ((((2*sigma*mu-(log(log(2))))^2)-(log(log(2)))^2)/(4*sigma))/y[i]))))
}
 #prior's distribuition
mu ~ dgamma(0.001,0.001)
sigma ~ dgamma(0.001,0.001)
}", file="WF_mu4.txt")


y <- t
n <- length(y)
cens <- cens

data.D4 <- list(n=n,  y=y, cens=cens)

inits1 <-list("mu"=214, "sigma"=0.002) 
inits2 <- list("mu"=214, "sigma"=0.0005) 

bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(872)
Modelo.FW <- jags(data=data.D4, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                  n.chains=2, n.iter=60000,n.burnin = 20000,
                  n.thin =10,
                  model.file="WF_mu4.txt")

#Outputs
print(Modelo.FW)


#DIC  c=100
DIC <- Modelo.FW$BUGSoutput$DIC
DIC

#DIC corrected
DIC.c <- (DIC-2*n*100)
DIC.c

#Deviance 
d.bar <- mean(Modelo.FW$BUGSoutput$sims.list$deviance)
d.bar

#Deviance corrected
d.bar.c <- (d.bar-2*n*100)
d.bar.c

#pD
Modelo.FW$BUGSoutput$pD

DIC-d.bar
DIC.c - d.bar.c

#Statistics EAIC
p <- 2 #n?mero de par?metros do modelo

EAIC.c <- d.bar.c +2*p
EAIC.c 

EBIC.c <- d.bar.c +p*log(n)
EBIC.c 


#m?dia em rela??o a cadeia de mu e sigma
mean(Modelo.FW$BUGSoutput$sims.list$mu)
mean(Modelo.FW$BUGSoutput$sims.list$sigma)

Modelo.FW$BUGSoutput$sd
Modelo.FW$BUGSoutput$
  

quantile(Modelo.FW$BUGSoutput$sims.list$sigma)
hist(Modelo.FW$BUGSoutput$sims.list$sigma,freq=F)
lines(density(Modelo.FW$BUGSoutput$sims.list$sigma))
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

ll <- log(log(2))
b <- 
cdf <-  1- 

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
l[i] <-  cens[i]*log((sigma/(mu/(exp(loggam((1/sigma)+1)))))*((y[i]/(mu/(exp(loggam((1/sigma)+1)))))^(sigma-1))*
exp(-(y[i]/(mu/(exp(loggam((1/sigma)+1)))))^sigma)) +
(1-cens[i])*log( exp(-(y[i]/(mu/(exp(loggam((1/sigma)+1)))))^sigma) ) 
}
 #prior's distribuition
mu ~ dgamma(0.001,0.001)
sigma ~ dgamma(0.001,0.001)
}", file="WEI4.txt")

y <- t
n <- length(y)
cens <- cens

data.D4 <- list(n=n,  y=y, cens=cens)

inits1 <-list("mu"=214, "sigma"=0.02) 
#inits2 <- list("mu"=mean(D3), "sigma"=0.1) 
inits2 <- list("mu"=214, "sigma"=0.05) 
#bayes.mod.inits <- list(inits1, inits2, inits3)
bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(742)
Modelo.WEI <- jags(data=data.D4, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                   n.chains=2, n.iter=60000,n.burnin = 20000,
                   n.thin =10,
                   model.file="WEI4.txt")

#Outputs
print(Modelo.WEI)


#valor de DIC  c=100
DIC1 <- Modelo.WEI$BUGSoutput$DIC
DIC1
#valor de DIC corrigido menos c=100
DIC.c1 <- (DIC1-2*n*100)
DIC.c1

#Deviance m?dio
d.bar1 <- mean(Modelo.WEI$BUGSoutput$sims.list$deviance)
d.bar1
#Deviance m?dio corrigido
d.bar.c1 <- (d.bar1-2*n*100)
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
