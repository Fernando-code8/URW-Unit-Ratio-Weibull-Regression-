################################################################################
# PAPER: The COVID-19 mortality rate in Latin America: a cross-country analysis
# SECTION: 3. NUMERICAL EVIDENCE
# GOAL: Conducting a simulation study for the URW regression model with variable 
# shape parameter.
# AUTHORS: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#          Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
################################################################################

rm(list=ls()) 
source("URW_reg.R")
library("moments")
library("gamlss")
library("beepr")

vn<-c(30,70,150,300)
alpha10<-.10
alpha5<-.05
alpha1<-.01
R=10000
set.seed(2020) 
tempo <- Sys.time()
for ( n in vn){
  print(n)
  x1<-runif(n)
  x2<-runif(n)
  X<-cbind(x1)
  X2<-cbind(rep(1,length(x1)),x2)
  X1<-cbind(rep(1,length(x1)),X)
  
  # beta<-c(-3.75,0.25)
  beta <- c(-5,1.75)
  eta <- X1%*%beta
  mu <- make.link("logit")$linkinv(eta)
  # beta2<-c(1.5,1.5)
  beta2<-c(2,1.25)
  eta2 <- X2%*%beta2
  sigma <- make.link("identity")$linkinv(eta2)
  tau=.75
  i <- 1
  c_conv <- 0
  est <- c()
  k<-length(c(beta,beta2))
  LI<-LS<-EP<-vres<-CAres<-Kres<-mres<-pv<-est
  vetor=rep(0, k)
  normres5=normres1=normres10=c_conv
  while(i <= R){
    print(i)
    y <-rURW(n,mu,sigma)
    fit <- try(gamlss(y~X, sigma.formula=~X2[,-1],
                      family = URW(
                        mu.link = "logit",sigma.link = "identity"
                      ), trace= F),silent = T )
    
    if(class(fit)[1] != "try-error"){
      est<-rbind(est,c(fit$mu.coefficients,fit$sigma.coefficients))
      mres[i]<-mean(fit$residuals)
      vres[i]<-var(fit$residuals)
      CAres[i]<-skewness(fit$residuals)
      Kres[i]<-kurtosis(fit$residuals)
      
      ## saving stderror
      a<-summary(fit)
      EP<-rbind(EP,a[,2])
      
      # Shapiro-Wilk normality test for residuals
      teste<-shapiro.test(fit$residuals)
      pv<- teste$p.value
      normres10=(pv<alpha10)+normres10
      normres5=(pv<alpha5)+normres5
      normres1=(pv<alpha1)+normres1
      
      ## Confidence interval
      Li=c(fit$mu.coefficients,fit$sigma.coefficients)-qnorm(1-alpha5/2)*a[,2]
      Ls=c(fit$mu.coefficients,fit$sigma.coefficients)-qnorm(alpha5/2)*a[,2]
      theta=c(beta,beta2)
      vetor<-(Ls<=theta)+(Li>=theta)+vetor
      LI<-rbind(LI,Li)
      LS<-rbind(LS,Ls)
      
      i<-i+1
      
    }else{
      c_conv <- c_conv + 1 
      print("No Convergence")
    }
  }
  
  parameters <- c(beta,beta2)
  media<-colMeans(est)
  vies<- media-parameters
  RB<-100*vies/parameters
  CA<-apply(est,2,skewness)
  K<-apply(est,2,kurtosis)
  dp<-apply(est,2,sd)
  CV<-100*dp/media
  EQM<- dp^2+ vies^2
  Taxa<-1-(vetor/R) # coverage rate
  
  # Analise descritiva dos Residuos ##colocar assimetria e curtose
  media_res<-mean(mres)
  var_res<-mean(vres)
  CAres<-mean(CAres)
  Kres<-mean(Kres)
  teste_res10<-normres10/R
  teste_res5<-normres5/R
  teste_res1<-normres1/R
  
  resultado<- rbind(
    parameters, RB, EQM,CA,K,Taxa
  )
  
  resultado_res<-rbind(media_res, var_res, CAres, Kres,teste_res10,teste_res5,teste_res1)
  
  row.names(resultado)<-c("Parameters","RB%", "MSE","CS","K", "CR")
  row.names(resultado_res)<-c("Mean","Variance", "CS", "K",
                              "NRR(10%)", "NRRR(5%)", "NRR(1%)")
  colnames(resultado_res)<-c("Residuals")
  print(c("Tamanho amostral=",n))
  print(round(resultado,4))
  print(round(t(resultado_res),4))
  name<-paste0("Scenario2_n",n,"_R",R,".RData")
  save.image(name)
}
