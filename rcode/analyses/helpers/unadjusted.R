#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function perform unadjusted analysis 
#------------------------------------------------------------------------------#

analyse_unadjusted <- function(data){
  pY1 <- mean(data$Y[data$A==1])
  pY0 <- mean(data$Y[data$A==0])
  
  MRR <- pY1/pY0
  MOR <- (pY1 * (1-pY0)) / (pY0 * (1-pY1))

  return(list(MRR=MRR,
              MOR=MOR))
}


