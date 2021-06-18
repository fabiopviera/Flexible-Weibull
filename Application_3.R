###### Aplicando o modelo Weibull Poisson

require(coda)
require(R2jags)
require(rjags)
#require(runjags)
library(mcmc)
library(mcmcplots)
require(LaplacesDemon)
require(e1071)
require(gamlss)
require(gamlss.cens)
require(MCMCvis)
require(superdiag)
require(survival)
require(muhaz)
require(bshazard)

dados <- read.table("dados_voltagem.txt",header = T)

dados$kV<- as.factor(dados$kV)
str(dados)

plot(dados$tempo,dados$kV)
ss1 <- dados$tempo[dados$kV=="52.5"]
dd1 <- dados$cens[dados$kV=="52.5"]

ss2 <- dados$tempo[dados$kV=="55"]
dd2 <- dados$cens[dados$kV=="55"]

ss3 <- dados$tempo[dados$kV=="57.5"]
dd3 <- dados$cens[dados$kV=="57.5"]

sd(ss1[dd1==1])
sd(ss2[dd2==1])
sd(ss3[dd3==1])
max(dados$tempo)
mean(dados$tempo[dados$kV=="55"])


fitW.Null <- gamlss(Surv(tempo,cens)~1,
                sigma.formula =~1,family=cens(WEI3)
                ,n.cyc=200,data=dados)
summary(fitW.Null)
fitW.Media <- gamlss(Surv(tempo,cens)~kV,
                sigma.formula =~1,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Media)
fitW.Shape <- gamlss(Surv(tempo,cens)~1,
                sigma.formula =~kV,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Shape)
fitW.Full <- gamlss(Surv(tempo,cens)~kV,
                sigma.formula =~kV,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Full)

plot(fitW.Full)


fitFW.Null <- gamlss(Surv(tempo,cens)~1,sigma.formula =~ 1,family=cens(FW1),
                n.cyc=200,data=dados,sigma.start = 0.0003)
summary(fitFW.Null,type = "qr")
fitFW.Median <- gamlss(Surv(tempo,cens)~kV,sigma.formula =~ 1,family=cens(FW1),
                n.cyc=200,data=dados,sigma.start = 0.0003)
summary(fitFW.Median)
fitFW.Shape <- gamlss(Surv(tempo,cens)~1,sigma.formula =~ kV,family=cens(FW1),
                n.cyc=200,data=dados,sigma.start = fitFW.Null$sigma.fv)
summary(fitFW.Shape)
fitFW.Full <- gamlss(Surv(tempo,cens)~kV,sigma.formula =~ kV,family=cens(FW1),
                n.cyc=200,data=dados,sigma.start = 0.0003,mu.start = fitFW.Median$mu.fv)
summary(fitFW.Full,type = "qr")

summary(fitFW.Full,type = "qr", 
        robust=T,  hessian.fun = "R")


plot(fitFW.Full)


confint(fitFW.Full)

deviance(fit10);deviance(fit11)
AIC(fit10,fit11)
BIC(fit10,fit11)

Res.qFW <- fitFW.Full$residuals

x11()
qqnorm(Res.qFW,pch=19,ylim = c(-3,3),xlim=c(-3,3),lwd=4,
       ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.2)
qqline(Res.qFW, datax = F, distribution = qnorm,
       col="gray65",lwd=4)

(ekm <- survfit(Surv(tempo,cens,type ="right" )~kV,data=dados))
summary(ekm)


plot(ekm,conf.int=F,mark.time = T,lwd=2,col=c(1,2,3,4,5))

x11()
plot(ekm,mark.time = T,lwd=c(4,4,4),lty=c(1,2,1),col=c("gray45","gray45","gray45"),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="S(x)")

curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=T,lwd=4,col="darkorchid",lty=1)
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="chartreuse3")
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="firebrick2")

curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="darkorchid",lty=2)
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="chartreuse3")
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="firebrick2")
legend("topright",c("Empirical survival","FW1 - Level: 52.5","FW1 - Level: 55.0","FW1 - Level: 57.5","WEI3 - Level: 52.5","WEI3 - Level: 55.0","WEI3 - Level: 57.5"),lty=c(1,1,1,1,2,2,2),lwd=c(3,3,3,3,3,3,3),col=c("gray45","darkorchid","chartreuse3","firebrick2","darkorchid","chartreuse3","firebrick2"),bty='n',cex=1.2)

x11()
plot(ekm,mark.time = T,lwd=c(4,4,4),lty=c(1,2,1),col=c("gray45","gray45","gray45"),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="S(x)")

curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=1)
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="gray1")
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="gray75")

curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=2)
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="gray1")
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="gray75")
legend(4500,1,c("Empirical survival","FW1 - Level: 52.5kV","FW1 - Level: 55.0kV","FW1 - Level: 57.5kV","WEI3 - Level: 52.5kV","WEI3 - Level: 55.0kV","WEI3 - Level: 57.5kV"),lty=c(1,1,1,1,2,2,2),lwd=c(3,3,3,3,3,3,3),col=c("gray45","gray35","gray1","gray75","gray35","gray1","gray75"),bty='n',cex=1.2)


x11()
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=F,lwd=4,col="gray35",lty=1,xlim=c(0,500),ylim=c(0.4,1))
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="gray1")
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="gray75")

curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=2)
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="gray1")
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="gray75")
legend(4500,1,c("Empirical survival","FW1 - Level: 52.5kV","FW1 - Level: 55.0kV","FW1 - Level: 57.5kV","WEI3 - Level: 52.5kV","WEI3 - Level: 55.0kV","WEI3 - Level: 57.5kV"),lty=c(1,1,1,1,2,2,2),lwd=c(3,3,3,3,3,3,3),col=c("gray45","gray35","gray1","gray75","gray35","gray1","gray75"),bty='n',cex=1.2)



par(new=TRUE, oma=c(3,1,1,2))
## create a layout to plot the subplot in the right bottom corner
layout(matrix(1:4,2))
## use xlim and ylim to zoom the subplot
plot(y ~ x, data = lin,xlim=c(0,2), ylim=c(0,2))
abline(linm)

x11()
plot(ekm,mark.time = T,lwd=c(4,4,4),lty=c(1,2,1),col=c("gray45","gray45","gray45"),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="S(x)",
     xlim=c(0,500),ylim=c(0.4,1))

curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=1)
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="gray1")
curve(1-pFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="gray75")

curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=2)
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="gray1")
curve(1-pWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="gray75")
legend("topright",c("Empirical survival","FW1 - Level: 52.5kV","FW1 - Level: 55.0kV","FW1 - Level: 57.5kV","WEI3 - Level: 52.5kV","WEI3 - Level: 55.0kV","WEI3 - Level: 57.5kV"),lty=c(1,1,1,1,2,2,2),lwd=c(3,3,3,3,3,3,3),col=c("gray45","gray35","gray1","gray75","gray35","gray1","gray75"),bty='n',cex=1.2)




x11()

curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=F,lwd=4,col="darkorchid",lty=1,xlim = c(0.001,6500),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="h(x)",ylim=c(0,0.002))
curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="chartreuse3",xlim = c(0.001,3500))
curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="firebrick2",xlim = c(0.001,1500))

curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="darkorchid",lty=2)
curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="chartreuse3",xlim = c(0.001,3500))
curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="firebrick2",xlim = c(0.001,1500))
legend("topright",c("FW_mu - Level: 52.5","FW_mu - Level: 55.0","FW_mu - Level: 57.5","Weibull - Level: 52.5","Weibull - Level: 55.0","Weibull - Level: 57.5"),lty=c(1,1,1,2,2,2),lwd=c(3,3,3,3,3,3),col=c("darkorchid","chartreuse3","firebrick2","darkorchid","chartreuse3","firebrick2"),bty='n',cex=1.2)


x11()
curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=F,lwd=4,col="gray35",lty=1,xlim = c(0.001,6500),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="h(x)",ylim=c(0,0.002))
curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="gray1",xlim = c(0.001,4500))
curve(hFW1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="gray75",xlim = c(0.001,1500))

curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]),                         sigma=exp(fitW.Full$sigma.coefficients[1])),add=T,lwd=4,col="gray35",lty=2)
curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2])),add=T,lwd=4,lty=2,col="gray1",xlim = c(0.001,4500))
curve(hWEI3(x,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[3])),add=T,lwd=4,lty=2,col="gray75",xlim = c(0.001,1500))
legend("topright",c("FW1 - Level: 52.5kV","FW1 - Level: 55.0kV","FW1 - Level: 57.5kV","WEI3 - Level: 52.5kV","WEI3 - Level: 55.0kV","WEI3 - Level: 57.5kV"),lty=c(1,1,1,2,2,2),lwd=c(3,3,3,3,3,3),col=c("gray35","gray1","gray75","gray35","gray1","gray75"),bty='n',cex=1.2)



hWEI3 <-function(x,mu,sigma) dWEI3(x,mu=mu,sigma=sigma)/(1-pWEI3(x,mu=mu,sigma=sigma))

#########################################################

s1 <- 1- pFW1(dados$tempo,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1]))
predict_dist <- data.frame(
  time = c(0,sort(dados$tempo)),
  pred_d1 = c(1,s1))


t1 = p_(ggsurv(ekm, plot.cens = F, CI = F, size.est = 2))
p_(t1 +
  geom_line(data = predict_dist,
            mapping = aes(x = time, y = pred_d1, colour = 'D1'), 
            size = 1, col = 'red')) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, .2)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  annotate("text",
           label  = c("Kaplan-Meier", "PW"), 
           x = c(70, 77),
           y = c(1, 0.92),
           col = c('black', 'red'), size   = 6) + 
  xlab('Time') + ylab('Survival Function') +
  theme(legend.position = 'none') +
  theme_gray((base_size=15))
t1


cbind(EM.FW=  c(
        exp(fitFW.Full$mu.coefficients[1]),
        exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),
        exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3])),
      EM.W=  c(
        exp(fitW.Full$mu.coefficients[1]),
        exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]),
        exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[3]))
)



x11()
plot(density(fitFW.Full$residuals),lwd=5,xlab="Quantile Residuals",ylab="Density",main="",col="1",cex.lab=1.2,cex.axis=1.2)
rug(jitter(fitFW.Full$residuals))

x11()
plot(density(fitW.Full$residuals),lwd=5,xlab="Quantile Residuals",ylab="Density",main="",col="1",cex.lab=1.2,cex.axis=1.2)
rug(jitter(fitW.Full$residuals))

x11()
wp1(fitFW.Full,cor1 = "gray65",ylim.all = 2)
wp(fitFW.Full,ylim.all = 2)


x11()
wp1(fitW.Full,cor1 = "gray65",ylim.all = 2)



plot(ekm,mark.time = T,xlim=c(0,1000),lwd=c(4,4,4),lty=c(1,2,1),col=c("gray45","gray45","gray45"),cex.lab=1.2,cex.lab=1.2,xlab ="Time (in minutes)",ylab="S(x)")

curve(1-pWF1(x,mu=exp(fitFW.Full$mu.coefficients[1]),                         sigma=exp(fitFW.Full$sigma.coefficients[1])),add=T,lwd=4,col="darkorchid",lty=1)
curve(1-pWF1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2])),add=T,lwd=4,col="chartreuse3")
curve(1-pWF1(x,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[3]),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[3])),add=T,lwd=4,col="firebrick2")
abline(v=250,lwd=3)
abline(v=300,lwd=3)
abline(v=350,lwd=3)


Res.q1 <- fitFW.Full$residuals
Res.qo <- sort(Res.q1)
index <- 1:length(dados$tempo)
plot(index,Res.q1, ylim=c(-4,4))

iter<-0  
j<-1
B<-200
n <- length(dados$tempo)
set.seed(78)
mrq <- matrix(0,  ncol = B, nrow = n)
while(j<B+1){
  simula <- rWF1(n, fitted(fitFW.Full),fitted(fitFW.Full,what = "sigma"))
  m1s <- try(gamlss(Surv(simula,cens) ~kV,sigma.formula=~kV,family=cens(WF1),n.cyc=200,data=dados,
                    sigma.start =0.00003))
  if((class(m1s) != "try-error")==T){
    
    Res.qs <- m1s$residuals
    mrq[,j] <- Res.qs
    j=j+1
  }
  cat("iteration = ", iter <- iter + 1, j,"\n")
}

dim(mrq)
Res.q.Ords <- apply(mrq,2,sort)

i <- 1:n
Z <- qnorm((i - 3/8) / (n + 1/4))
#Z <- qnorm((i - 1/8) / (2*n + 1/2))


rqi.m  <-  apply(Res.q.Ords,1, mean)
rqi.min <- apply(Res.q.Ords,1, min)
rqi.max <- apply(Res.q.Ords,1, max)


(res.out <- sum(c(sum(Res.qo>rqi.max),sum(Res.qo<rqi.min))))
(per.out <- round(res.out/n*100,2))

x11()
qqnorm(Res.qo,pch=19,col="gray35",ylim = c(-3.5,3.5),xlim=c(-3.5,3.5),
       ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.35)
lines(Z,rqi.max,col="gray20",lwd=3)
lines(Z,rqi.m, lty = 2,col="gray20",lwd=3)
lines(Z,rqi.min,col="gray20",lwd=3)
legend("topleft", c(paste("Total points:",n), paste("Points out of envelope:",res.out,"(",per.out,"%)")), bty="n", cex=1.25)


Res.q <- fitW.Full$residuals
Res.q <- sort(Res.q)
index <- 1:length(dados$tempo)
plot(index,Res.q, ylim=c(-4,4))

iter<-0  
j<-1
B<-200
n <- length(dados$tempo)
set.seed(78)
mrq <- matrix(0,  ncol = B, nrow = n)
while(j<B+1){
  simula <- rWEI3(n, fitted(fitW.Full),fitted(fitW.Full,what = "sigma"))
  m1s <- try(gamlss(Surv(simula,cens) ~kV,sigma.formula=~kV,family=cens(WEI3),n.cyc=200,data=dados))
  if((class(m1s) != "try-error")==T){
    
    Res.qs <- m1s$residuals
    mrq[,j] <- Res.qs
    j=j+1
  }
  cat("iteration = ", iter <- iter + 1, j,"\n")
}

dim(mrq)
Res.q.Ords <- apply(mrq,2,sort)

i <- 1:n
Z <- qnorm((i - 3/8) / (n + 1/4))
#Z <- qnorm((i - 1/8) / (2*n + 1/2))


rqi.m  <-  apply(Res.q.Ords,1, mean)
rqi.min <- apply(Res.q.Ords,1, min)
rqi.max <- apply(Res.q.Ords,1, max)


(res.out <- sum(c(sum(Res.q>rqi.max),sum(Res.q<rqi.min))))
(per.out <- round(res.out/n*100,2))

x11()
qqnorm(Res.q,pch=19,col="gray35",ylim = c(-3.5,3.5),xlim=c(-3.5,3.5),
       ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.35)
lines(Z,rqi.max,col="gray20",lwd=3)
lines(Z,rqi.m, lty = 2,col="gray20",lwd=3)
lines(Z,rqi.min,col="gray20",lwd=3)
legend("topleft", c(paste("Total points:",n), paste("Points out of envelope:",res.out,"(",per.out,"%)")), bty="n", cex=1.25)



Res.qFW <- fitFW.Full$residuals
Res.qW <- fitW.Full$residuals

x11()
qqnorm(Res.qFW,pch=19,ylim = c(-3,3),xlim=c(-3,3),lwd=4,
       ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.2)
qqline(Res.qFW, datax = F, distribution = qnorm,
       col="gray65",lwd=4)

# grafico 1
x11()
qqnorm(Res.qW,pch=19,ylim = c(-3,3),xlim=c(-3,3),lwd=4,
       ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.2)
qqline(Res.qW, datax = F, distribution = qnorm,
       col="gray65",lwd=4)

wp(fitFW.Full)
wp(fitW.Full)

wp(resid=rd)
wp(fitW.Full)


#################################################
d1 <- c(rep(0,20),rep(1,20),rep(0,20))
d2 <- c(rep(0,20),rep(0,20),rep(1,20))

rm <- 1 + log(1-pWF1(dados$tempo,mu=exp(fitFW.Full$mu.coefficients[1]+fitFW.Full$mu.coefficients[2]*d1+fitFW.Full$mu.coefficients[3]*d2),sigma=exp(fitFW.Full$sigma.coefficients[1]+fitFW.Full$sigma.coefficients[2]*d1+fitFW.Full$sigma.coefficients[3]*d2)))
rd  <- sign(rm)*(-2 *(rm+1*log(1-rm)))^0.5

Resw <- 1 + log(1-pWEI3(dados$tempo,mu=exp(fitW.Full$mu.coefficients[1]+fitW.Full$mu.coefficients[2]*d1+fitW.Full$mu.coefficients[3]*d2),sigma=exp(fitW.Full$sigma.coefficients[1]+fitW.Full$sigma.coefficients[2]*d1+fitW.Full$sigma.coefficients[3]*d2)))
Resdw <- sign(Resw)*(-2 *(Resw+1*log(1-Resw)))^0.5


x11()
qqnorm(rd,pch=19,ylim = c(-3,3),xlim=c(-3,3),lwd=4,
       ylab="Deviance residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.2)
points()
qqline(rd, datax = F, distribution = qnorm,
       col="gray65",lwd=4)




# grafico 1
x11()
qqnorm(Resdw,pch=19,ylim = c(-3,3),xlim=c(-3,3),lwd=4,
       ylab="Deviance residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.2)
qqline(Resdw, datax = F, distribution = qnorm,
       col="gray65",lwd=4)
#legend("topleft",c("n=50","30% censored"),col=c(NA,NA),bty = "n",cex=1.2,seg.len = 1.5)



#############################################################
###################################
##An?lise Bayesiana

cat("
data{
 C <- 10000
 for (i in 1:n){ zeros[i] <- 0}
}
model{
for (i in 1:n){
zeros[i] ~ dpois(phi[i])
phi[i] <- - l[i] + C
l[i] <-  cens[i]*log((exp((sigma[i]*y[i] - ((((2*sigma[i]*mu[i]-(log(log(2))))^2)
-(log(log(2)))^2)/(4*sigma[i]))/y[i])-exp(sigma[i]*y[i] - 
((((2*sigma[i]*mu[i]-(log(log(2))))^2)-
(log(log(2)))^2)/(4*sigma[i]))/y[i])))*
(sigma[i]+((((2*sigma[i]*mu[i]-(log(log(2))))^2)-(log(log(2)))^2)/
(4*sigma[i]))/(y[i]^2)))    +     (1-cens[i])*log((exp(-exp(sigma[i]*y[i] - ((((2*sigma[i]*mu[i]-(log(log(2))))^2)-(log(log(2)))^2)/(4*sigma[i]))/y[i]))))
log(mu[i]) <- X[i,1]*beta10 + beta11*X[i,2]+ beta12*X[i,3]
log(sigma[i]) <- X[i,1]*beta20 + beta21*X[i,2]+ beta22*X[i,3]
}
beta10 ~ dnorm(0, 0.0001)
beta11 ~ dnorm(0, 0.0001)
beta12 ~ dnorm(0, 0.0001)
beta20 ~ dnorm(0, 0.0001)
beta21 ~ dnorm(0, 0.0001)
beta22 ~ dnorm(0, 0.0001)

}", file="WF_mu_regre1.txt")

#curve(dnorm(x,0,0.0001),xlim=c(-0.01,0.01))
y <- dados$tempo
cens <- dados$cens

d0 <- rep(1,length(y))
d1 <- c(rep(0,20),rep(1,20),rep(0,20))
d2 <- c(rep(0,20),rep(0,20),rep(1,20))
X<- cbind(d0,d1,d2)
#J<- ncol(X)
n <- length(y)


data.regre <- list(n=n,y=y,X=X, cens=cens)

#inits1 <- list("beta10"=3,"beta11"=0.5,"beta12"=0.3, "beta20"=-8)
#inits1 <- list("beta10"=3, "beta20"=-8,"beta21"=0.5,"beta22"=2)
inits1 <- list("beta10"=3,"beta11"=0.5,"beta12"=0.3, "beta20"=-8,"beta21"=0.8,"beta22"=2)
#inits2 <- list("mu"=mean(D3), "sigma"=0.1) 
inits2 <- list("beta10"=5,"beta11"=0.7,"beta12"=1.3, "beta20"=-7,"beta21"=.7,"beta22"=0.3) 
#inits2 <- list("beta10"=5,"beta11"=0.7,"beta12"=1.3, "beta20"=-15) 
#inits2 <- list("beta10"=5, "beta20"=-15,"beta21"=0.7,"beta22"=2.3) 

#bayes.mod.params <- c("beta10", "beta11","beta12","beta20")
#bayes.mod.params <- c("beta10","beta20", "beta21","beta22")
bayes.mod.params <- c("beta10","beta11","beta12","beta20", "beta21","beta22")
bayes.mod.inits <- list(inits1,inits2)

set.seed(987)
Modelo.FW <- jags(data=data.regre, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                  n.chains=2, n.iter=60000,n.burnin = 20000,
                  n.thin =10,
                  model.file="WF_mu_regre1.txt")

#Outputs
print(Modelo.FW)


#valor de DIC  c=100
c=10000
DIC <- Modelo.FW$BUGSoutput$DIC
DIC
#valor de DIC corrigido 
DIC.c <- (DIC-2*n*c)
DIC.c

#Deviance m?dio
d.bar <- mean(Modelo.FW$BUGSoutput$sims.list$deviance)
d.bar
#Deviance m?dio corrigido
d.bar.c <- (d.bar-2*n*c)
d.bar.c

#pD
Modelo.FW$BUGSoutput$pD

DIC-d.bar
DIC.c - d.bar.c

#Estatistica EAIC


p <- length(Modelo.FW$BUGSoutput$mean)- 1 #n?mero de par?metros do modelo

EAIC.c <- d.bar.c +2*p
EAIC.c 

EBIC.c <- d.bar.c +p*log(n)
EBIC.c 



  
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

x11()
MCMCplot(Modelo.FW, 
         object2 = Modelo.WEI,
         params = c('beta10','beta11','beta12','beta20','beta21','beta22'),
         offset=0.15,
         col="blue",col2="orange")


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
 C <- 10000
 for (i in 1:n){ zeros[i] <- 0}
}
model{
for (i in 1:n){
zeros[i] ~ dpois(phi[i])
phi[i] <- - l[i] + C
l[i] <-  cens[i]*log((sigma[i]/(mu[i]/(exp(loggam((1/sigma[i])+1)))))*((y[i]/(mu[i]/(exp(loggam((1/sigma[i])+1)))))^(sigma[i]-1))*
exp(-(y[i]/(mu[i]/(exp(loggam((1/sigma[i])+1)))))^sigma[i])) +
(1-cens[i])*log( exp(-(y[i]/(mu[i]/(exp(loggam((1/sigma[i])+1)))))^sigma[i]) ) 
log(mu[i]) <- X[i,1]*beta10 + beta11*X[i,2]+ beta12*X[i,3]
log(sigma[i]) <- X[i,1]*beta20  + beta21*X[i,2]+ beta22*X[i,3]
}
beta10 ~ dnorm(0, 0.0001)
beta11 ~ dnorm(0, 0.0001)
beta12 ~ dnorm(0, 0.0001)
beta20 ~ dnorm(0, 0.0001)
beta21 ~ dnorm(0, 0.0001)
beta22 ~ dnorm(0, 0.0001)
}", file="WEI_regre2.txt")

#curve(dnorm(x,0,0.0001),xlim=c(-0.01,0.01))
y <- dados$tempo
cens <- dados$cens

d0 <- rep(1,length(y))
d1 <- c(rep(0,20),rep(1,20),rep(0,20))
d2 <- c(rep(0,20),rep(0,20),rep(1,20))
X<- cbind(d0,d1,d2)
#J<- ncol(X)
n <- length(y)


data.regre2 <- list(n=n,y=y,X=X, cens=cens)

#inits1 <- list("beta10"=5,"beta11"=1.5,"beta12"=0.5, "beta20"=-0.5)
#inits1 <- list("beta10"=3, "beta20"=-3,"beta21"=-0.5,"beta22"=1)
inits1 <- list("beta10"=3,"beta11"=1.5,"beta12"=0.5, "beta20"=-3,"beta21"=-0.5,"beta22"=1)
#inits2 <- list("mu"=mean(D3), "sigma"=0.1) 
inits2 <- list("beta10"=5,"beta11"=0.7,"beta12"=0.5, "beta20"=-3,"beta21"=-0.7,"beta22"=1.3) 
#inits2 <- list("beta10"=3,"beta11"=0.7,"beta12"=0.5, "beta20"=-0.7) 
#inits2 <- list("beta10"=5, "beta20"=-3,"beta21"=-0.7,"beta22"=1.3) 

#bayes.mod.params <- c("beta10", "beta11","beta12","beta20")
#bayes.mod.params <- c("beta10","beta20", "beta21","beta22")
bayes.mod.params <- c("beta10","beta11","beta12","beta20", "beta21","beta22")
bayes.mod.inits <- list(inits1, inits2)


set.seed(982)
Modelo.WEI <- jags(data=data.regre2, inits=bayes.mod.inits,parameters.to.save = bayes.mod.params,
                   n.chains=2, n.iter=60000,n.burnin = 20000,
                   n.thin =10,
                   model.file="WEI_regre2.txt")

#Outputs
print(Modelo.WEI)


#valor de DIC  c=10000
c=10000
DIC1 <- Modelo.WEI$BUGSoutput$DIC
DIC1
#valor de DIC corrigido menos c=100
DIC.c1 <- (DIC1-2*n*c)
DIC.c1

#Deviance m?dio
d.bar1 <- mean(Modelo.WEI$BUGSoutput$sims.list$deviance)
d.bar1
#Deviance m?dio corrigido
d.bar.c1 <- (d.bar1-2*n*c)
d.bar.c1

#pD
Modelo.WEI$BUGSoutput$pD

DIC1-d.bar1
DIC.c1 - d.bar.c1

#Estatistica EAIC
p <- length(Modelo.WEI$BUGSoutput$mean)- 1 #n?mero de par?metros do modelo

EAIC.c1 <- d.bar.c1 +2*p
EAIC.c1 

EBIC.c1 <- d.bar.c1 +p*log(n)
EBIC.c1 


#m?dia em rela??o a cadeia de mu e sigma
mean(Modelo.WEI$BUGSoutput$sims.list$beta10)
mean(Modelo.WEI$BUGSoutput$sims.list$beta20)


hist(Modelo.WEI$BUGSoutput$sims.list$beta10,freq=F)
lines(density(Modelo.WEI$BUGSoutput$sims.list$beta10))
curve(dnorm(x,0,0.0001),add=T,lwd=3)


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

require(MCMCvis)

x11()
MCMCplot(Modelo.WEI, 
         params = c('beta10','beta11','beta12','beta20','beta21','beta22'))


traceplot(Modelo.WEI.mcmc)
x11()
denplot(Modelo.WEI.mcmc)
x11()
autocorr.plot(Modelo.WEI.mcmc)

x11()
gelman.plot( Modelo.WEI.mcmc )
geweke.diag( Modelo.WEI.mcmc )
