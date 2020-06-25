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
                     bAL4,
                     bYL1,
                     bYL2,
                     bYL3,
                     bYL4,
                     eventrate,
                     Yint,
                     sd_UY,
                     rhoL,
                     seed){
  # Generate data
  set.seed(seed)
  
  # Compute nobs from set number of events and varying eventrate
  nobs         <- round(nevents * 1/eventrate, digits = 0)
  
  # Generate Y, A and L
  UY           <- rnorm(nobs,sd_UY)
  
  sigL         <- matrix(rhoL, 
                         nrow = nL,
                         ncol = nL)
  diag(sigL)   <- 1
  L <- mvrnorm(nobs, mu=rep(0, times = nL), Sigma = sigL)
  
  # Store covariate effects
  betasAL <- c(rep(bL_const, times=(nL/2)), # half of covariates are fixed confounders
               rep(bAL1,times=(nL/8)),      # other half varies across scenarios
               rep(bAL2,times=(nL/8)),      # (predictor/noise/instrument/confounder)
               rep(bAL3,times=(nL/8)),
               rep(bAL4,times=(nL/8)))
  
  betasYL <- c(rep(bL_const, times=(nL/2)), # half of covariates are fixed confounders
               rep(bYL1,times=(nL/8)),      # other half varies across scenarios
               rep(bYL2,times=(nL/8)),      # (predictor/noise/instrument/confounder)
               rep(bYL3,times=(nL/8)),
               rep(bYL4,times=(nL/8)))
  
  # Generate exposure:
  A <- rbinom(nobs, 1, plogis(L %*% betasAL))
  
  # Generate outcome:
  Y <- rbinom(nobs, 1, plogis(Yint + 
                                A * bYA + 
                                L %*% betasYL +
                                UY))
  
  out_df <- data.frame(Y,A,L)
  colnames(out_df) <- c("Y","A",paste0("L",1:nL))
  
  return(out_df)
}


