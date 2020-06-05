#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function to obtain Marginal Risk Ratio
#------------------------------------------------------------------------------#

# Estimate Marginal Risk Ratio from logistic regression models ----
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

# Obtain true Marginal Risk Ratio by numerical integration ----
#------------------------------------------------------------------------------#

integrate_MRR <- function(nL, bYA, bL_const, bAL1, bAL2, bAL3, bYL1, bYL2, bYL3, Yint, rhoL){
  var.L1.star <- (nL/2) + (nL/2)*((nL/2)-1)*rhoL
  var.L2.star <- (nL/6) + (nL/6)*((nL/6)-1)*rhoL
  var.L3.star <- (nL/6) + (nL/6)*((nL/6)-1)*rhoL
  var.L4.star <- (nL/6) + (nL/6)*((nL/6)-1)*rhoL
  cov.L       <- (nL/2)*(nL/6)*(nL/6)*(nL/6)*rhoL
  sigma       <- matrix(cov.L,ncol=4,nrow=4)
  diag(sigma) <- c(var.L1.star,var.L2.star,var.L3.star,var.L3.star)
  dx          <- 0.2
  values      <- seq(-4,4,by=dx)
  x           <- expand.grid(values,values,values,values)
  dL          <- dmvnorm(x,mean=c(rep(0,times=4)), sigma=sigma)
  
  pY0     		<- plogis((Yint + bL_const*x[,1] + 
                         bYL1*x[,2] + bYL2*x[,3] + 
                         bYL3*x[,4]))*(dL*dx) 
  pY1     		<- plogis((Yint + bL_const*x[,1] + 
                         bYL1*x[,2] + bYL2*x[,3] + 
                         bYL3*x[,4] + 
                         bYA))*(dL*dx) 
  
  MRR <- mean(pY1)/mean(pY0)
  MOR <- (mean(pY1) * (1- mean(pY0)))/((1-mean(pY1)) * mean(pY0))
  
  
  return(list(MRR=MRR, MOR=MOR))
}
