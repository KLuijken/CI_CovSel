#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper functions to process model output
#------------------------------------------------------------------------------#

# Obtain warnings (convergence, separation)
obtain_warnings <- function(premodel){
  warning_mod <- ifelse(is.null(premodel$warning),NA,paste(premodel$warning))
  
  return(warning_mod)
}

# Store model coefficients and number of variables
obtain_coefficients <- function(model,
                                data,
                                datagen_scenario = datagen_scenario){
  coef_est    <- coef(model)[c("(Intercept)","A",
                                 paste0("L",seq(1:datagen_scenario[['nL']])))]
  coef_est[is.na(coef_est)] <- NA
    
  coef_se     <-  rep(NA, times= datagen_scenario[['nL']]+2) # plus exposure and intercept
  names(coef_se) <- c("(Intercept)","A",
                      paste0("L",seq(1:datagen_scenario[['nL']])))
  coef_se[c(names(model$coefficients))] <- sqrt(diag(vcov(model, complete=T)))
  names(coef_se) <- c("se(Intercept)","se_A",
                      paste0("se_L",seq(1:datagen_scenario[['nL']])))
  
  
  r_squared_M <- lm(model$linear.predictors ~ data[,"Y"])
  r_squared   <- r_squared_M$coefficients[2]
    
  alloutput   <- c(coef_est,coef_se,r_squared)
  
  return(alloutput)
}
