########################################################################################
# PAPER: The COVID-19 mortality rate in Latin America: a cross-country analysis
# SUBSECTION: 4.2. Regression results and discussion
# GOAL: Providing functions to fit a KW, UW and URW regression model
# AUTHOR: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#         Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
########################################################################################

#Packages
library(nortest)
library(gamlss)
library(extRemes)

####################################################################################
#                        IT FITS A UNIT REGRESSION MODEL
####################################################################################
tau = 0.5

#Main function
UnitReg.fit <- function(z, X, Y, n = NA, regression = "URW"){ 

  source("URW_reg.R")
  source("UW_reg.R")
  source("KW_reg.R")
  
  if (any(regression == c("URW", "Kuma", "UW"))){
    if (regression == "URW"){
      source("URW_reg.R")
      mod <- gamlss(z~X[,-1], sigma.formula=~Y[,-1],
                    family = URW(
                      mu.link = "logit",sigma.link = "identity"
                    ), trace= F)
    } 
    if (regression == "Kuma"){
      source("KW_reg.R")
      mod <- gamlss(z~X[,-1], sigma.formula=~Y[,-1],
                    family = Kuma(
                      mu.link = "logit",sigma.link = "identity"
                    ), trace= F)
    } 
    if (regression == "UW"){
      source("UW_reg.R")
      mod <- gamlss(z~X[,-1], sigma.formula=~Y[,-1],
                    family = UW(
                      mu.link = "logit",sigma.link = "identity"
                    ), trace= F,c.crit=.0001, n.cyc=1000)
    } 
    
    #Quantile residuals
    mu_fv <- as.numeric(mod$mu.fv)
    sigma_fv <- as.numeric(mod$sigma.fv)
    quant_res_mod <- mod$residuals 
    
    #goodness-of-fit measures
    good_fit_meas_mat <- matrix(NA, nrow = 1,ncol = 3)
    colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(AD)")
    rownames(good_fit_meas_mat) <- c(regression)
    good_fit_meas_mat[1,] <- c( AIC(mod,k=2),
                                Rsq(mod),
                                as.numeric(ad.test(quant_res_mod)[2]))
    
    #summary BETA
    K <- ncol(X)+length(mod$sigma.coefficients)
    mod.coef<-summary(mod)
    mles_mod<- mod.coef[,1]
    se_mod<- mod.coef[,2]
    t_value <- mod.coef[,3]
    p_value_mod <- mod.coef[,4]
    
    # if(length(mod$sigma.coefficients)==1){sigma<-c("Sigma")}else{sigma<-names(data.frame(Y))}
    
    #SEs
    summ_MOD = matrix(NA,nrow = K, ncol = 4)
    row.names(summ_MOD) <- c(names(data.frame(X)),names(data.frame(Y)))
    colnames(summ_MOD) <- c("Estimate", "Std. Error", "t value","Pr(>|t|)")
    summ_MOD[,1] <- mles_mod
    summ_MOD[,2] <- se_mod
    summ_MOD[,3] <- t_value
    summ_MOD[,4] <- p_value_mod
    #round(summ_BETA, digits = 4)
    
    log.lik.mod = gen.likelihood(mod)
    ell_mod = -log.lik.mod()
    
    results = list(summary = summ_MOD, 
                   mod.coeff = mles_mod[1:(K-(length(mod$sigma.coefficients)))],
                   sigma.coeff = mles_mod[((K+1)-(length(mod$sigma.coefficients))):K],
                   mu.fv = mu_fv, 
                   sigma.fv = sigma_fv, 
                   residuals = quant_res_mod,
                   diagnosticMEASURES = good_fit_meas_mat,
                   max_ell = ell_mod,
                   eta_hat = X%*%mles_mod[1:(K-(length(mod$sigma.coefficients)))]
    )
    
    return(results)
  }
  else {
    stop(paste(regression, "Regression not available, available regressions are \"URW\",","\"Kuma\", and \"UW\""))
  }
}

####################################################################################
#                                 LOOCV FUNCTION
####################################################################################
LOOCV.unit_reg <- function(z_comp,X_comp,Y_comp, regression = "Kuma"){
  
  if (any(regression == c("URW", "Kuma", "UW"))){
    data_comp = data.frame(cbind(z_comp,X_comp[,-1],Y_comp[,-1]))
    n_comp.CV <- length(data_comp$z_comp)
    X_comp.CV <- model.matrix(z_comp~.,data = data_comp[,1:(ncol(data_comp)-ncol(Y_comp[,-1]))]) 
    z_comp.CV <- data_comp$z_comp
    Y_comp.CV <- model.matrix(z_comp~.,data = data_comp[,(ncol(X_comp.CV)+1):length(data_comp)])
    n_train.CV <- n_comp.CV-1
    # LOOCV LOOP
    mae_vec.CV = vector() 
    bug<-0
    for(c in 1:n_comp.CV){
      
      # Salve the cth observation
      z_test.CV <- z_comp.CV[c]
      X_test.CV <- X_comp.CV[c,]
      Y_test.CV<- Y_comp.CV[c,]
      
      # Extract the cth observation
      z_train.CV <- z_comp.CV[-c]
      X_train.CV <- X_comp.CV[-c,]
      Y_train.CV<- Y_comp.CV[-c,]
      
     
      if (regression == "URW"){
        #Fit the model
        mod_train <- try(gamlss(z_train.CV~X_train.CV[,-1],sigma.formula=~Y_train.CV[,-1], 
                                     family = URW(
                                       mu.link = "logit",sigma.link = "identity"
                                     ), trace= F),T)}
      if (regression == "Kuma"){
        #Fit the model
        mod_train <- try(gamlss(z_train.CV~X_train.CV[,-1],sigma.formula=~Y_train.CV[,-1], 
                                     family = Kuma(
                                       mu.link = "logit",sigma.link = "identity"
                                     ), trace= F),T)}
      if (regression == "UW"){
        #Fit the model
        mod_train <- try(gamlss(z_train.CV~X_train.CV[,-1],sigma.formula=~Y_train.CV[,-1], 
                                     family = UW(
                                       mu.link = "logit",sigma.link = "identity"
                                     ), trace= F,c.crit=.0001, n.cyc=1000),T)}
      if(class(mod_train)[1] != "try-error"){
        ## prediction
        source("Link_function.R")
        vec_hat.CV = as.numeric(mod_train$mu.coefficients)
        z_predict.CV = lfunc(c(vec_hat.CV),X_test.CV)

        # lfunc<-make.link("logit") # define a fun??o de liga??o logit
        # z_predict.CV = lfunc$linkfun(vec_hat.CV) # aplica a fun??o de liga??o logit em y
      
        # adequacy measures
        mae.CV = mean(abs(z_test.CV-z_predict.CV))
        
        # salve the adeq. measures
        mae_vec.CV[c] = mae.CV
      }else{
        bug<-bug+1
      }
      }
      
    a<-cbind(rep(1:length(mae_vec.CV)),mae_vec.CV)
    a<-as.data.frame(a)
    a <- na.omit(a)
    a$mae_vec.CV<-as.numeric(a$mae_vec.CV)
    #Results matrix
    res_CV = matrix(NA, nrow = 1,ncol = 1)
    colnames(res_CV) <- c("MAE")
    res_CV[1,1] <- c(mean(a$mae_vec.CV))
    
    return(res_CV)
  }
    else {
      stop(paste(regression, "Regression not available, available regressions are \"URW\",","\"Kuma\", and \"UW\""))
    }
}

################################################################################
#                                RESIDUALS PLOT FUNCTION                       #
################################################################################

plots_quantres <- function(quant_residuals, fitted_values, sample_size){
  ##PLOT 2: qrs versus index
  plot(1:n,quant_residuals, main = "",      #We expected no trend.
       xlab = "Index", pch=20,ylim=c(-2.5,2.5),
       ylab = "Quantile residuals")  
  abline(h=-2, lty=2)
  abline(h=0)
  abline(h=2,lty=2)
  
  ##PLOT 1: Worm plot (wp)/
  wp(resid = quant_residuals, col = 1, ylim.all =2)
  title(main = "")
  
  ##PLOT 3: NORMAL QQ-PLOT
  qqnorm(quant_residuals,pch=1,frame=T, main = "",
         make.plot = T, lwd=1,ylim=c(-2.5,2.5))
}