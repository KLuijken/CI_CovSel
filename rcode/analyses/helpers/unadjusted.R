#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function perform unadjusted analysis 
#------------------------------------------------------------------------------#

analyse_unadjusted <- function(data, method){
  # conditional estimator
  model <- logistf(Y~A,
          data = data,
          firth = isTRUE(method == "FLIC"),
          pl = F)
  COR <- coef(model)["A"]
  
  # marginal estimator
  pY1 <- mean(data$Y[data$A==1])
  pY0 <- mean(data$Y[data$A==0])
  
  MRR <- pY1/pY0
  MOR <- (pY1 * (1-pY0)) / (pY0 * (1-pY1))

  return(list(COR = COR,
              MRR = MRR,
              MOR = MOR))
}


