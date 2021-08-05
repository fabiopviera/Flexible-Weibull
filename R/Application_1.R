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


#Descriptive statistics
summary(D3)
sd(D3)
skewness(D3)
kurtosis(D3)
length(D3)
mean(D3)
sd(D3)


# FW1 Model
fit1 <-gamlss(D3~1,family=FW1(mu.link = "identity",sigma.link = "identity"),n.cyc = 1200,
              sigma.start = 0.0018)

summary(fit1)

prof.dev(fit1, which="mu", min=0.1, max=1)
prof.dev(fit1, which="sigma", min=0.0005, max=0.5)



# WEI3 Model
fit2 <-gamlss(D3~1,family=WEI3(mu.link = "identity",sigma.link = "identity"),n.cyc = 1200)

summary(fit2)

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
inits2 <- list("mu"=3, "sigma"=1.2) 

bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(943)
Modelo.FW <- jags(data=data3, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                   n.chains=2, n.iter=60000,n.burnin = 20000,
                   n.thin =10,
                   model.file="WF_mu3.txt")

#Outputs
print(Modelo.FW)


#DIC  c=100
DIC <- Modelo.FW$BUGSoutput$DIC
DIC

#DIC adjusted
DIC.c <- (DIC-2*23*100)
DIC.c

#Deviance 
d.bar <- mean(Modelo.FW$BUGSoutput$sims.list$deviance)
d.bar

#Deviance adjusted
d.bar.c <- (d.bar-2*23*100)
d.bar.c

#pD
Modelo.FW$BUGSoutput$pD

DIC-d.bar
DIC.c - d.bar.c

#Estatistics EAIC
#EAIC
p <- 2 #parameters

EAIC.c <- d.bar.c +2*p
EAIC.c 

EBIC.c <- d.bar.c +p*log(n)
EBIC.c 


#mu and sigma
mean(Modelo.FW$BUGSoutput$sims.list$mu)
mean(Modelo.FW$BUGSoutput$sims.list$sigma)



####################################################################
# Bayesian WEI3

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

inits2 <- list("mu"=3, "sigma"=1.2) 

bayes.mod.params <- c("mu", "sigma")
bayes.mod.inits <- list(inits1, inits2)

set.seed(563)
Modelo.WEI <- jags(data=data3, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                  n.chains=2, n.iter=60000,n.burnin = 20000,
                  n.thin =10,
                  model.file="WEI3.txt")

#Outputs
print(Modelo.WEI)


#DIC  c=100
DIC1 <- Modelo.WEI$BUGSoutput$DIC
DIC1

#DIC adjusted
DIC.c1 <- (DIC1-2*23*100)
DIC.c1

#Deviance 
d.bar1 <- mean(Modelo.WEI$BUGSoutput$sims.list$deviance)
d.bar1

#Deviance adjusted
d.bar.c1 <- (d.bar1-2*23*100)
d.bar.c1

#pD
Modelo.WEI$BUGSoutput$pD

DIC1-d.bar1
DIC.c1 - d.bar.c1

#Estatistics EAIC
p <- 2 #parameters

EAIC.c1 <- d.bar.c1 +2*p
EAIC.c1 

EBIC.c1 <- d.bar.c1 +p*log(n)
EBIC.c1 


#mu and sigma
mean(Modelo.WEI$BUGSoutput$sims.list$mu)
mean(Modelo.WEI$BUGSoutput$sims.list$sigma)

