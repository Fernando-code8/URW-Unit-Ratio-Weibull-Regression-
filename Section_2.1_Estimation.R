################################################################################
# PAPER: The COVID-19 mortality rate in Latin America: a cross-country analysis
# SUBSECTION: 2.1. Parameter estimation
# GOAL: Calculation of the likelihood, vector score and Fisher's information.
# AUTHOR: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#         Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
################################################################################

####################################
##      log-likelihood   URW     ##
####################################
sigma<- 1
mu<- 0.6
tau<-0.5
  
URWd<-function(y){sigma*y^(sigma-1)*(1-mu)^sigma/(mu^sigma*(1-y)^(sigma+1))*
    log((1-tau)^(-1))*(1-tau)^(y^sigma*(1-mu)^sigma/((mu^sigma)*(1-y)^sigma))}
integrate(URWd,0,1)

y=c(.1,.2,.3,.4)
prod(URWd(y))

llURW<- function(y){
  n<-length(y)
  exp(sum(log(sigma/y))+sum(log(log(1-tau)/(y-1)))
      +sum(log(((y*(1-mu))/(mu*(1-y)))^sigma))
      +log(1-tau)*sum(((y*(1-mu))/(mu*(1-y)))^sigma))

  }

llURW(y)

################################################################################
#################       DERIVATIVES  URW                       #################
################################################################################
set.seed(2020)
n=100
u = runif(n)
y = (mu^sigma*log(1-u)/((1-mu)^sigma*log(1-tau)))^(1/sigma)/
  (1+(mu^sigma*log(1-u)/((1-mu)^sigma*log(1-tau)))^(1/sigma))

### LOGLIKEHOOD
loglik<-expression((log(sigma/y))+(log(log(1-tau)/(y-1)))
                   +(log(((y*(1-mu))/(mu*(1-y)))^sigma))
                   +log(1-tau)*(((y*(1-mu))/(mu*(1-y)))^sigma))

#############################################
############## FIRST DERIVATION #############
#############################################

### mu 

## R function
dmu<-D(loglik,"mu")
sum(eval(dmu))  

## Derivaton 
k=(y*(mu-1))/(mu*(y-1))
sum((sigma/(mu*(mu-1)))*(1+log(1-tau)*(k)^sigma))

### sigma

## R function
dsigma<-D(loglik,"sigma")
sum(eval(dsigma))  

## Derivation
sum(1/sigma+log(k)*(1+log(1-tau)*(k)^sigma))

#############################################
########### SECOND DERIVATION ###############
#############################################

### mu^2 

## R function
ddmu<-D(dmu,"mu")
sum(eval(ddmu)) 

## Derivation
sum((sigma/(mu^2*(mu-1)^2))*(1-2*mu+((k)^sigma)*log(1-tau)*(1-2*mu+sigma)
                               ))
### sigma^2

## R function
ddsigma<-D(dsigma,"sigma")
sum(eval(ddsigma)) 

## Derivation
sum(((k)^sigma)*(log(k)^2)*log(1-tau)
    -(1/(sigma^2)))

### mu sigma

## R function
ddmusigma<-D(dmu,"sigma")
sum(eval(ddmusigma)) 

## Derivation
sum((1/(mu*(mu-1)))*(1+((((k)^sigma)*log(1-tau)))*
      (1+sigma*log(k))))

