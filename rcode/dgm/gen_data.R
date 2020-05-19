#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Data generating mechanisms
#------------------------------------------------------------------------------#

# Some notes on how the data are generated:
# Each dataset contains a binary exposure A, a binary outcome Y and a set of 
# covariates, L, where the number of covariates, nL, can be specified by the user in
# sim_scen.R. Note that the number of covariates should be a multiplicative of 3,
# since there are three levels of association between a (set of) covariates and
# the exposure and outcome.


# Load libraries  ----
#------------------------------------------------------------------------------#

# Generate data  ----
#------------------------------------------------------------------------------#

# use_datagen_senarios
gen_data <- function(nobs,
                     nL,
                     bAL1,
                     bAL2,
                     bAL3,
                     bYL1,
                     bYL2,
                     bYL3,
                     sd_UY,
                     rhoL,
                     seed){
  # Generate data
  set.seed(seed)
  
  # Generate Y, A and L
  UY           <- rnorm(nobs,sd_UY)
  
  sigL         <- matrix(rhoL, 
                         nrow = nL,
                         ncol = nL)
  diag(sigL)   <- 1
  L <- mvrnorm(nobs, mu=rep(0, times = nL), Sigma = sigL)
  
  # Make code flexible for number of confounders.
  # Generate exposure:
  A <- rbinom(nobs, 1,
                         plogis(L[,1:(nL/3)] * bAL1 +
                                L[,(nL/3+1):(2*nL/3)] * bAL2 + 
                                L[,(2*nL/3+1):nL] * bAL3))
  # Generate outcome:
  Y <- rbinom(nobs, 1,
                         plogis(L[,1:(nL/3)] * bYL1 +
                                L[,(nL/3+1):(2*nL/3)] * bYL2 + 
                                L[,(2*nL/3+1):nL] * bYL3 + UY))
  
  out_df <- data.frame(Y,A,L)
  colnames(out_df) <- c("Y","A",paste0("L",1:nL))
  
  return(out_df)
}


