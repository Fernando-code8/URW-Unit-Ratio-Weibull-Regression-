###################################################################################
# PAPER: The COVID-19 mortality rate in Latin America: a cross-country analysis
# GOAL: Providing link function for the regression models.
# AUTHOR: Fernando Arturo Pe\~na-Ram\'irez, Renata Rojas Guerra and 
#         Fernando Jos\'e Monteiro de Ara\'ujo
# LAST UPDATE: October 09, 2024
###################################################################################

# Logit link function
lfunc <- function(beta_vec,X){
  qi_logit_link <- 1/(1+exp(-(X%*%beta_vec)))
  return(qi_logit_link)
}