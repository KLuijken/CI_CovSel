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


# Obtain warnings (convergence, separation) [note: not used yet]
obtain_warnings <- function(premodel, preintmodel, scen, rep){
  warning     <- ifelse(is.null(premodel$warning),"NULL",{paste0(premodel$warning,"scen=",scen,"nsim=",rep)}) # Store warning, if any, and location in simulation repetitions
  warning_int <- ifelse(is.null(preintmodel$warning),"NULL",{paste0(preintmodel$warning,"scen=",scen,"nsim=",rep)}) # Store warning, if any

  return(list(warning = warning, warning_int = warning_int))
}
  
# Store model coefficients and number of variables
obtain_coefficients <- function(model, intmodel,
                                datagen_scenario = datagen_scenario){
  int         <- intmodel$coefficients[1]
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
  
  return(alloutput)
}



