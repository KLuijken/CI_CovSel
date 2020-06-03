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
gen_data <- function(nevents,
                     nL,
                     bYA,
                     bL_const,
                     bAL1,
                     bAL2,
                     bAL3,
                     bYL1,
                     bYL2,
                     bYL3,
                     eventrate,
                     Yint,
                     sd_UY,
                     rhoL,
                     seed){
  # Generate data
  set.seed(seed)
  
  # Compute nobs from set number of events and varying eventrate
  nobs         <- nevents * 1/eventrate
  
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
                         # half of the covariates are fixed confounders:
                         plogis(L[,1:(nL/2)] %*% rep(bL_const,times=(nL/2)) + 
                                # other half of the covariates is a mix of noise/instruments/predictors/confounders
                                L[,((nL/2)+1):((nL/2)+(nL/6))] %*% rep(bAL1,times=(nL/6)) +
                                L[,((nL/2)+(nL/6)+1):((nL/2)+(2*nL/6))] %*% rep(bAL2,times=(nL/6)) + 
                                L[,((nL/2)+(2*nL/6)+1):nL] %*% rep(bAL3,times=(nL/6))))
  # Generate outcome:
  Y <- rbinom(nobs, 1,
                         plogis(Yint + 
                                bYA * A + 
                                # half of the covariates are fixed confounders:
                                L[,1:(nL/2)] %*% rep(bL_const, times=(nL/2)) +
                                # other half of the covariates is a mix of noise/instruments/predictors/confounders
                                L[,((nL/2)+1):((nL/2)+(nL/6))] %*% rep(bYL1,times=(nL/6)) +
                                L[,((nL/2)+(nL/6)+1):((nL/2)+(2*nL/6))] %*% rep(bYL2,times=(nL/6)) +
                                L[,((nL/2)+(2*nL/6)+1):nL] %*% rep(bYL3,times=(nL/6)) +
                                UY))
  
  out_df <- data.frame(Y,A,L)
  colnames(out_df) <- c("Y","A",paste0("L",1:nL))
  
  return(out_df)
}


