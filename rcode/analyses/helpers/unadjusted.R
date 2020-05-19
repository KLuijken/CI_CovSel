#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function perform unadjusted analysis [note: not added yet]
#------------------------------------------------------------------------------#

analyse_unadjusted <- function(data){
  pY1 <- mean(data$Y[data$A==1])
  pY0 <- mean(data$Y[data$A==0])
  
  MRR <- mean(pY1)/mean(pY0)
  MOR <- (mean(pY1) * (1-mean(pY0))) / (mean(pY0) * (1-mean(pY1)))

  return(list(MRR=MRR,
              MOR=MOR))
}


