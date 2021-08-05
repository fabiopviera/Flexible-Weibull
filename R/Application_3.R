#### Application 3

require(gamlss)
require(gamlss.cens)
require(survival)

source("R/FW1_GAMLSS.R")


data.1 <- read.table("dados_voltagem.txt",header = T)

data.1$kV<- as.factor(data.1$kV)



# FW1 Model
fitFW.Null <- gamlss(Surv(tempo,cens)~1,sigma.formula =~ 1,family=cens(FW1),
                n.cyc=200,data=data.1,sigma.start = 0.0003)
summary(fitFW.Null,type = "qr")

fitFW.Median <- gamlss(Surv(tempo,cens)~kV,sigma.formula =~ 1,family=cens(FW1),
                n.cyc=200,data=data.1,sigma.start = 0.0003)
summary(fitFW.Median, type="qr")

fitFW.Shape <- gamlss(Surv(tempo,cens)~1,sigma.formula =~ kV,family=cens(FW1),
                n.cyc=200,data=data.1,sigma.start = fitFW.Null$sigma.fv)
summary(fitFW.Shape, type="qr")

fitFW.Full <- gamlss(Surv(tempo,cens)~kV,sigma.formula =~ kV,family=cens(FW1),
                n.cyc=200,data=data.1,sigma.start = 0.0003,mu.start = fitFW.Median$mu.fv)
summary(fitFW.Full,type = "qr")


# WEI3 model
fitW.Null <- gamlss(Surv(tempo,cens)~1,
                    sigma.formula =~1,family=cens(WEI3)
                    ,n.cyc=200,data=data.1)
summary(fitW.Null, type="qr")

fitW.Media <- gamlss(Surv(tempo,cens)~kV,
                     sigma.formula =~1,family=cens(WEI3),n.cyc=200,data=data.1)
summary(fitW.Media, type="qr")

fitW.Shape <- gamlss(Surv(tempo,cens)~1,
                     sigma.formula =~kV,family=cens(WEI3),n.cyc=200,data=data.1)
summary(fitW.Shape, type="qr")

fitW.Full <- gamlss(Surv(tempo,cens)~kV,
                    sigma.formula =~kV,family=cens(WEI3),n.cyc=200,data=data.1)
summary(fitW.Full, type="qr")



