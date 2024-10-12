urw<-expression(log(
  sigma*y^(sigma-1)*(1-mu)^sigma/(mu^sigma*(1-y)^(sigma+1))*
    log((1-tau)^(-1))*(1-tau)^(y^sigma*(1-mu)^sigma/((mu^sigma)*(1-y)^sigma))
)
)
m1urw<-D(urw,"mu")
m2urw<-D(m1urw,"mu")
s1urw<-D(urw,"sigma")
s2urw<-D(s1urw,"sigma")
ms2urw<-D(m1urw,"sigma")


URW<-function (mu.link = "logit", sigma.link = "identity") 
{
  tau<-.5
  mstats <- checklink("mu.link", "URW", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "URW", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("URW", "URWeibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
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
                 
                 dldm = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldm <- eval(m1urw)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1urw)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldd <- eval(s1urw)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau<-.5
                   dldd <- eval(s1urw)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1urw)
                   dldd <- eval(s1urw)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dURW(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pURW",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(     mu <- rep(median(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(.3, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}


# density function
dURW<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  fy1 <-   sigma*y^(sigma-1)*(1-mu)^sigma/(mu^sigma*(1-y)^(sigma+1))*
    log((1-tau)^(-1))*(1-tau)^(y^sigma*(1-mu)^sigma/((mu^sigma)*(1-y)^sigma))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function
pURW<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  cdf1<- 1-(1-tau)^(q^sigma*(1-mu)^sigma/(mu^sigma*(1-q)^sigma))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#------------------------------------------------------------------------------------------ #ok
# quantile function
qURW<-function(u,mu,sigma)
{
  tau<-.5
  q<- (mu^sigma*log(1-u)/((1-mu)^sigma*log(1-tau)))^(1/sigma)/
    (1+(mu^sigma*log(1-u)/((1-mu)^sigma*log(1-tau)))^(1/sigma))
  q
}

# inversion method for randon generation
rURW<-function(n,mu,sigma,a=0,b=1)
{
  u<- runif(n)
  tau<-.5
  y<- qURW(u,mu =mu, sigma =sigma)
  y
}
