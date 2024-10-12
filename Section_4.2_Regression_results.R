################################################################################
# PAPER:The COVID-19 mortality rate in Latin America: a cross-country analysis
# SUBSECTION: Regression results and discussion
# GOAL: Fiting the  KW, UW and URW regression models
# AUTHOR: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#         Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
################################################################################

rm(list = ls())
library(tidyverse)
library(readxl)
library(gamlss)
library(nortest)

###############################  AUXILIARY FUNCTIONS  ##########################
source("KW_reg.R")
source("UW_reg.R")
source("URW_reg.R")
source("unit_regressions_fit.R")

### results matrix
good_fit_meas_mat <- matrix(NA, nrow = 3,ncol = 6)
colnames(good_fit_meas_mat) <- c("AIC","BIC","R-square","p-val(AD)", "MAE","MAE/MR")
rownames(good_fit_meas_mat)<- c("Kw","URW","UW")

###############################  DATA PREPARATION  #############################
data.set <- read_excel("dados_final.xlsx") %>%
  mutate(GDP=GDP/1000000)
# Response and sample size
z <- data.set$MR
n <- length(z)

#Selected covariates
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    GDP,
    HDI,
    CHE,
    DP,
    DGGHE
  )

data.set.select2 <- data.set%>%   
  dplyr::select(
    MR,
    GDP,
    UP
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  
Y <- model.matrix(MR~.,data = data.set.select2)  
tau<-.5

###############################  FITTED REGRESSIONS  ###########################

################################################################################
################################# Kumaraswamy ##################################
################################################################################
mod_Kuma<-gamlss(z~X,sigma.formula = ~Y, family="Kuma")
# fitted coefficients
Kuma.coef<-summary(mod_Kuma)
# quantile residuals
quant_res_Kuma <- mod_Kuma$residuals

# MAE/mean(MR)
(res_KW=round(LOOCV.unit_reg(z,X,Y, regression = "Kuma"),4))
mae_KW = as.numeric(res_KW)
mp_KW<-(comp_KW = round(mae_KW/mean(z),4))

# goodness-of-fit measures
good_fit_meas_mat[1,] <- c( AIC(mod_Kuma,k=2),BIC(mod_Kuma),
                            Rsq(mod_Kuma),
                            as.numeric(ad.test(quant_res_Kuma)[2]),mae_KW,mp_KW)


plots_quantres(quant_res_Kuma,mod_Kuma$mu.fv,n)

################################################################################
###################################### URW #####################################
################################################################################

mod_URW<-gamlss(z~X,sigma.formula = ~Y, family="URW")
# fitted coefficients
URW.coef<-summary(mod_URW)
# quantile residuals
quant_res_URW <- mod_URW$residuals 

# MAE/mean(MR)
(res_URW=round(LOOCV.unit_reg(z,X,Y, regression = "URW"),4))
mae_URW = as.numeric(res_URW)
mp_URW<-(comp_URW = round(mae_URW/mean(z),4))

# goodness-of-fit measures
good_fit_meas_mat[2,] <- c( AIC(mod_URW,k=2),BIC(mod_URW),
                            Rsq(mod_URW),
                            as.numeric(ad.test(quant_res_URW)[2]),mae_URW,mp_URW)

plots_quantres(quant_res_URW,mod_URW$mu.fv,n)


################################################################################
###################################### UW ######################################
################################################################################

mod_UW<-gamlss(z~X,sigma.formula = ~Y, family="UW")
mod_UW<-UnitReg.fit(z, X, Y, n = NA, regression = "UW")
# fitted coefficients
UW.coef<-summary(mod_UW)
# quantile residuals
quant_res_UW <- mod_UW$residuals 

# MAE/mean(MR)
(res_UW=round(LOOCV.unit_reg(z,X,Y, regression = "UW"),4))
mae_UW = as.numeric(res_UW)
mp_UW<-(comp_UW = round(mae_UW/mean(z),4))

# goodness-of-fit measures
good_fit_meas_mat[3,] <- c( AIC(mod_UW,k=2),BIC(mod_UW),
                            Rsq(mod_UW),
                            as.numeric(ad.test(quant_res_UW)[2]),mae_UW,mp_UW)

plots_quantres(quant_res_UW,mod_UW$mu.fv,n)


### Results tables
round(cbind(Kuma.coef[,c(1,4)],URW.coef[,c(1,4)],UW.coef[,c(1,4)]),4)
round(good_fit_meas_mat,4)

