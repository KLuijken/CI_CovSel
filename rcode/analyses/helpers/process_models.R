#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper functions to process model output
#------------------------------------------------------------------------------#

# Store try_catch values and warnings in objects
pre_model <- function(inputmodel,
                      datagen_scenario = datagen_scenario,
                      data = data){
  preM      <- inputmodel
  M         <- preM$value
  Mpred     <- as.matrix(data[,names(coef(M)[-1])]) %*% coef(M)[-1]
  preM_int  <- tryCatch.W.E(glm(data$Y ~ offset(Mpred), family = binomial))   # Re-estimate the intercept
  M_int     <- preM_int$value
  
  return(list(preM = preM, M = M, preM_int = preM_int, M_int = M_int))
}


# Obtain warnings (convergence, separation)
obtain_warnings <- function(premodel, preintmodel){
  warning_mod <- ifelse(is.null(premodel$warning),NA,paste(premodel$warning))
  warning_int <- ifelse(is.null(preintmodel$warning),NA,paste(preintmodel$warning))
  
  return(list(warning_mod = warning_mod, warning_int = warning_int))
}

# Store model coefficients and number of variables
obtain_coefficients <- function(warning, model, intmodel,
                                datagen_scenario = datagen_scenario){
  ifelse(is.na(warning$warning_mod),{
    int         <- ifelse(is.na(warning$warning_int),intmodel$coefficients[1],model$coefficients[1])
    se_int      <- coef(summary(intmodel))[1,2]
    coef_est    <- coef(model)[c("(Intercept)","A",
                                 paste0("L",seq(1:datagen_scenario[['nL']])))]
    coef_est[is.na(coef_est)] <- NA
    coef_est    <- coef_est[-1] 
    
    coef_se     <-  rep(NA, times= datagen_scenario[['nL']]+2)
    names(coef_se) <- c("(Intercept)","A",
                        paste0("L",seq(1:datagen_scenario[['nL']])))
    coef_se[c(names(model$coefficients))] <- sqrt(diag(vcov(model, complete=T)))
    coef_se     <- coef_se[-1]
    
    alloutput   <- c(int,coef_est,se_int,coef_se)
  },
  {
    alloutput   <- rep(NA, times= (2*(datagen_scenario[['nL']]+2)))
  })
  
  
  return(alloutput)
}
