#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function to obtain Marginal Risk Ratio and Marginal Odds Ratio
#------------------------------------------------------------------------------#

# Estimate Marginal Ratios from logistic regression models ----
#------------------------------------------------------------------------------#

estimate_marginals <- function(warning, data, int, modelcoefs){
  ifelse(all(is.na(warning)),{
    newdat <- data.frame(data[,c(names(modelcoefs))])
    colnames(newdat) <- c(names(modelcoefs))
    
    newdat[,"A"] <- 0
    PredA0 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
    newdat[,"A"] <- 1
    PredA1 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
    
    MRR <- mean(PredA1)/mean(PredA0)
    MOR <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))
  },
  {
    MRR <- MOR <- NA
  })
    
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
  
  # Make code flexible for number of confounders.
  # Generate outcome:
  Y0 <- plogis(Yint + 
                # half of the covariates are fixed confounders:
                L[,1:(nL/2)] %*% rep(bL_const, times=(nL/2)) +
                # other half of the covariates is a mix of noise/instruments/predictors/confounders
                L[,((nL/2)+1):((nL/2)+(nL/8))] %*% rep(bYL1,times=(nL/8)) +
                L[,((nL/2)+(nL/8)+1):((nL/2)+(2*nL/8))] %*% rep(bYL2,times=(nL/8)) +
                L[,((nL/2)+(2*nL/8)+1):((nL/2)+(3*nL/8))] %*% rep(bYL3,times=(nL/8)) +
                L[,((nL/2)+(3*nL/8)+1):nL] %*% rep(bYL4,times=(nL/8)))
  
  Y1 <- plogis(Yint + 
                bYA + 
                # half of the covariates are fixed confounders:
                L[,1:(nL/2)] %*% rep(bL_const, times=(nL/2)) +
                # other half of the covariates is a mix of noise/instruments/predictors/confounders
                L[,((nL/2)+1):((nL/2)+(nL/8))] %*% rep(bYL1,times=(nL/8)) +
                L[,((nL/2)+(nL/8)+1):((nL/2)+(2*nL/8))] %*% rep(bYL2,times=(nL/8)) +
                L[,((nL/2)+(2*nL/8)+1):((nL/2)+(3*nL/8))] %*% rep(bYL3,times=(nL/8)) +
                L[,((nL/2)+(3*nL/8)+1):nL] %*% rep(bYL4,times=(nL/8)))
  
  MRR <- mean(Y1)/mean(Y0)
  MOR <- (mean(Y1) * (1- mean(Y0)))/((1-mean(Y1)) * mean(Y0))
  
  
  return(list(MRR=MRR, MOR=MOR))
}