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

# Obtain true Marginal Ratios by numerical integration ----
#------------------------------------------------------------------------------#

# integrate_marginals <- function(nL, bYA, bL_const, bYL1, bYL2, bYL3, bYL4, Yint, rhoL){
#   var.L1.star <- (nL/2) + (nL/2)*((nL/2)-1)*rhoL
#   var.L2.star <- (nL/8) + (nL/8)*((nL/8)-1)*rhoL
#   var.L3.star <- (nL/8) + (nL/8)*((nL/8)-1)*rhoL
#   var.L4.star <- (nL/8) + (nL/8)*((nL/8)-1)*rhoL
#   var.L5.star <- (nL/8) + (nL/8)*((nL/8)-1)*rhoL
#   cov.L       <- (nL/2)*(nL/8)*(nL/8)*(nL/8)*(nL/8)*rhoL
#   sigma       <- matrix(cov.L,ncol=5,nrow=5)
#   diag(sigma) <- c(var.L1.star,var.L2.star,var.L3.star,var.L4.star,var.L5.star)
#   dx          <- 0.4
#   values      <- seq(-4,4,by=dx)
#   x           <- as.matrix(expand.grid(values,values,values,values,values))
#   dL          <- mvnfast::dmvn(x,mu=c(rep(0,times=5)), sigma=sigma)
#   
#   pY0     		<- sum(plogis(Yint + bL_const*x[,1] + 
#                             bYL1*x[,2] + bYL2*x[,3] + 
#                             bYL3*x[,4] + bYL4*x[,5])*dL)/sum(dL)  
#   pY1     		<- sum(plogis(Yint + bL_const*x[,1] + 
#                             bYL1*x[,2] + bYL2*x[,3] + 
#                             bYL3*x[,4] + bYL4*x[,5] + bYA)*dL)/sum(dL) 
#   
#   MRR <- pY1/pY0
#   MOR <- (pY1 * (1- pY0))/((1-pY1) * pY0)
#   
#   
#   return(list(MRR=MRR, MOR=MOR))
# }

largesample_marginals <- function(nL, bYA, bL_const, bAL1, bAL2, bAL3, bAL4, bYL1, bYL2, bYL3, bYL4, Yint, rhoL){
  N <- 1e7
  
  sigL         <- matrix(rhoL, 
                         nrow = nL,
                         ncol = nL)
  diag(sigL)   <- 1
  L <- mvrnorm(N, mu=rep(0, times = nL), Sigma = sigL)
  
  # Make code flexible for number of confounders.
  # Generate exposure:
  A <- rbinom(N, 1,
              # half of the covariates are fixed confounders:
              plogis(L[,1:(nL/2)] %*% rep(bL_const,times=(nL/2)) + 
                       # other half of the covariates is a mix of noise/instruments/predictors/confounders
                       L[,((nL/2)+1):((nL/2)+(nL/8))] %*% rep(bAL1,times=(nL/8)) +
                       L[,((nL/2)+(nL/8)+1):((nL/2)+(2*nL/8))] %*% rep(bAL2,times=(nL/8)) + 
                       L[,((nL/2)+(2*nL/8)+1):((nL/2)+(3*nL/8))] %*% rep(bAL3,times=(nL/8)) + 
                       L[,((nL/2)+(3*nL/8)+1):nL] %*% rep(bAL4,times=(nL/8))))
  # Generate outcome:
  Y0 <- rbinom(N, 1,
              plogis(Yint + 
                       # half of the covariates are fixed confounders:
                       L[,1:(nL/2)] %*% rep(bL_const, times=(nL/2)) +
                       # other half of the covariates is a mix of noise/instruments/predictors/confounders
                       L[,((nL/2)+1):((nL/2)+(nL/8))] %*% rep(bYL1,times=(nL/8)) +
                       L[,((nL/2)+(nL/8)+1):((nL/2)+(2*nL/8))] %*% rep(bYL2,times=(nL/8)) +
                       L[,((nL/2)+(2*nL/8)+1):((nL/2)+(3*nL/8))] %*% rep(bYL3,times=(nL/8)) +
                       L[,((nL/2)+(3*nL/8)+1):nL] %*% rep(bYL4,times=(nL/8))))
  
  Y1 <- rbinom(N, 1,
              plogis(Yint + 
                       bYA * A + 
                       # half of the covariates are fixed confounders:
                       L[,1:(nL/2)] %*% rep(bL_const, times=(nL/2)) +
                       # other half of the covariates is a mix of noise/instruments/predictors/confounders
                       L[,((nL/2)+1):((nL/2)+(nL/8))] %*% rep(bYL1,times=(nL/8)) +
                       L[,((nL/2)+(nL/8)+1):((nL/2)+(2*nL/8))] %*% rep(bYL2,times=(nL/8)) +
                       L[,((nL/2)+(2*nL/8)+1):((nL/2)+(3*nL/8))] %*% rep(bYL3,times=(nL/8)) +
                       L[,((nL/2)+(3*nL/8)+1):nL] %*% rep(bYL4,times=(nL/8))))
  
  MRR <- mean(Y1)/mean(Y0)
  MOR <- (mean(Y1) * (1- mean(Y0)))/((1-mean(Y1)) * mean(Y0))
  
  
  return(list(MRR=MRR, MOR=MOR))
}