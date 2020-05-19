#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function to estimate Marginal Risk Ratio
#------------------------------------------------------------------------------#

estimate_marginals <- function(data, int, modelcoefs){
  newdat <- data.frame(data[,c(names(modelcoefs))])
  colnames(newdat) <- c(names(modelcoefs))
  
  ifelse("A" %in% c(names(modelcoefs)),{             # In order to catch boostrap samples that have too strong multicollinearity
    # Note to self: add counter to keep track in how many samples this occurs
    newdat[,"A"] <- rep(0,times=nrow(newdat))
    PredA0 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
    newdat[,"A"] <- rep(1,times=nrow(newdat))
    PredA1 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
    
    MRR <- mean(PredA1)/mean(PredA0)
    MOR <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))},
    
    MRR <- MOR <- NA)
    
  return(list(MRR = MRR, MOR = MOR))
}
