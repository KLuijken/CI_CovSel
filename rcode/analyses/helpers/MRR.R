#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function to obtain Marginal Risk Ratio and Marginal Odds Ratio
#------------------------------------------------------------------------------#

# Estimate Marginal Ratios from logistic regression models ----
#------------------------------------------------------------------------------#

estimate_marginals <- function(data, modelcoefs){
  newdat <- data.frame(cbind(1,
                               data[,c(names(modelcoefs[-1]))]))
  colnames(newdat) <- c(names(modelcoefs))
    
  newdat[,"A"] <- 0
  PredA0 <- plogis(as.matrix(newdat) %*% matrix(modelcoefs))
  newdat[,"A"] <- 1
  PredA1 <- plogis(as.matrix(newdat) %*% matrix(modelcoefs))
    
  MRR <- mean(PredA1)/mean(PredA0)
  MOR <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))

  return(list(MRR = MRR, MOR = MOR))
}

# Obtain true Marginal Ratios by large sample approximation ----
#------------------------------------------------------------------------------#


largesample_marginals <- function(nL, bYA, bL_const, bYL1, bYL2, bYL3, bYL4, Yint, rhoL){
  N <- 1e6
  
  sigL         <- matrix(rhoL, 
                         nrow = nL,
                         ncol = nL)
  diag(sigL)   <- 1
  L <- mvrnorm(N, mu=rep(0, times = nL), Sigma = sigL)
  
  # Store covariate effects
  betasYL <- c(rep(bL_const, times=(nL/2)), # half of covariates are fixed confounders
               rep(bYL1,times=(nL/8)),      # other half varies across scenarios
               rep(bYL2,times=(nL/8)),      # (predictor/noise/instrument/confounder)
               rep(bYL3,times=(nL/8)),
               rep(bYL4,times=(nL/8)))
  
  # Generate outcome:
  Y0 <- plogis(Yint +
                 L %*% betasYL)
  
  Y1 <- plogis(Yint + 
                bYA + 
                L %*% betasYL)
  
  MRR <- mean(Y1)/mean(Y0)
  MOR <- (mean(Y1) * (1- mean(Y0)))/((1-mean(Y1)) * mean(Y0))
  
  
  return(list(MRR=MRR, MOR=MOR))
}