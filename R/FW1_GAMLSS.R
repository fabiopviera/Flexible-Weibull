library(gamlss)
library(gamlss.dist)
FW1 <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link", "flexible Weibull 1", substitute(mu.link),c("1/mu^2", "log ", "identity"))
  dstats <- checklink("sigma.link", "flexible Weibull 1", substitute(sigma.link),c("inverse", "log", "identity"))
  structure(list(family = c("FW1 ", "flexible Weibull 1"),
                 parameters = list(mu=TRUE, sigma=TRUE), 
                 nopar = 2, 
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),  
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 
                 dldm = function(y,mu,sigma){ #----------------------------------------------------- ok
                   nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log = T),"mu", delta = 1e-04)
                   dldm = as.vector(attr(nd1, "gradient")) 
                   dldm
                 },
                 d2ldm2 = function(y,mu,sigma){#----------------------------------------------------- ok
                   nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log = T),"mu", delta = 1e-04)
                   dldm = as.vector(attr(nd1, "gradient")) 
                   d2ldm2 = -dldm * dldm
                 },     
                 dldd = function(y,mu,sigma){#----------------------------------------------------- ok  
                   nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log =T),"sigma", delta = 1e-04)
                   dldd = as.vector(attr(nd1, "gradient"))
                   dldd
                 } ,
                 d2ldd2 = function(y,mu,sigma){#----------------------------------------------------- ok
                   nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log =T),"sigma", delta = 1e-04)
                   dldd = as.vector(attr(nd1, "gradient"))
                   d2ldd2 = -dldd*dldd
                   d2ldd2 
                   },
                 d2ldmdd = function(y,mu,sigma){#----------------------------------------------------- ok
                     nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log= TRUE), "mu", delta = 1e-04)
                     dldm = as.vector(attr(nd1, "gradient"))
                     nd1 = gamlss:::numeric.deriv(dFW1(y, mu, sigma,log=TRUE), "sigma", delta = 1e-04)
                     dldd = as.vector(attr(nd1, "gradient"))           
                     d2ldmdd = -dldm * dldd
                     d2ldmdd               
                   },
                   
                 #----------------------------------------------------- ok  
                 G.dev.incr  = function(y,mu,sigma,...) -2*dFW1(y,mu=mu,sigma=sigma,log=TRUE), 
                 rqres = expression(rqres(pfun="pFW1", type="Continuous", y=y, mu=mu, sigma=sigma)),
                
                 mu.initial = expression( mu <- rep(median(y),length(y))),
                 sigma.initial = expression( sigma <- rep(min(y)/(10*max(y)),length(y))), 
                 
                 mu.valid = function(mu) all(mu > 0), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y>0)
  ),
  class = c("gamlss.family","family"))
}
#--------------------------------------------------------------
dFW1 <- function(x, mu=1, sigma=0.5, log = FALSE)
{
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x < 0))  stop(paste("x must be positive", "\n", "")) 
  
  ll <- log(log(2))
  b <- (((2*sigma*mu-ll)^2)-ll^2)/(4*sigma)
  f <-  (exp((sigma*x - b/x)-exp(sigma*x - b/x)))*(sigma+b/(x^2))
  
  if(log==FALSE) fx  <- f else fx <- log(f) 
  fx
}    
#--------------------------------------------------------------  
pFW1 <- function(q, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE)
{  
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
  
  ll <- log(log(2))
  b <- (((2*sigma*mu-ll)^2)-ll^2)/(4*sigma)
  cdf <-  1- (exp(-exp(sigma*q - b/q)))
  
  
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
  cdf
}
#--------------------------------------------------------------  
hFW1 <- function(x, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE)
{  
  h <- dFW1(x,mu,sigma)/(1-pFW1(x,mu,sigma))
  
  h
}


#--------------------------------------------------------------
qFW1 <- function(p, mu=1, sigma=0.5,  lower.tail = TRUE, log.p = FALSE )
{ 
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
  
  ll <- log(log(2))
  b <- (((2*sigma*mu-ll)^2)-ll^2)/(4*sigma)
  y <- (log(-log(1-p))   +  sqrt( ((log(-log(1-p)))^2) + 4*sigma*b))/(2*sigma)
  y
}
#--------------------------------------------------------------
rFW1 <- function(n, mu=1, sigma=0.5, nu=1)
{
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  n <- ceiling(n)
  u <- runif(n)
  r <- qFW1(u,mu=mu,sigma=sigma)
  r
}
#--------------------------------------------------------------
