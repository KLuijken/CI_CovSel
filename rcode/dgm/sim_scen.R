#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Simulation scenarios
#------------------------------------------------------------------------------#


# Load libraries  ----
#------------------------------------------------------------------------------#



# Initialize simulation scenarios ----
# Function output: data.frame with 8748 scenarios for the dgm, varying the 
# strength of association between covariates and exposure/outcome, residual variance
# of exposure and outcome, correlation between covariates, sample size and
# number of covariates.
#------------------------------------------------------------------------------#


datagen_scenarios <- function(){
  nobs  <- c(120,300)
  nL    <- 9
  bAL1  <- bAL2  <- bAL3  <- c(0,0.5,1)
  bYL1  <- bYL2  <- bYL3  <- c(0,0.5,1)
  Yint  <- 0
  sd_UY <- c(0.01,1)
  rhoL  <- c(0,0.3,0.7)
  
  # data.frame with simulation scenarios
  datagen_scenarios <- expand.grid(nobs,
                                   nL,
                                   bAL1,
                                   bAL2,
                                   bAL3,
                                   bYL1,
                                   bYL2,
                                   bYL3,
                                   Yint,
                                   sd_UY,
                                   rhoL)
  # Add scenarios number to keep track of results
  datagen_scenarios$scen_num <- c(1:nrow(datagen_scenarios))
  colnames(datagen_scenarios) <- c("nobs",
                                   "nL",
                                   "bAL1",
                                   "bAL2",
                                   "bAL3",
                                   "bYL1",
                                   "bYL2",
                                   "bYL3",
                                   "Yint",
                                   "sd_UY",
                                   "rhoL",
                                   "scen_num")
  
  return(datagen_scenarios)
}


# Initialize simulation analyses ----
# Specify penalisation in estimation no/yes (ML/FLIC) and p-value cutoff criterion
# in backward elimination.
# Function output: data.frame with 2 methods of analysis applied in the simulation
#------------------------------------------------------------------------------#

analysis_scenarios <- function(){
  method <- c("ML","FLIC")
  pcutoff <- 0.157
  
  analyse_scenarios <- expand.grid(method, pcutoff)
  colnames(analyse_scenarios) <- c("method","pcutoff")
  
  return(analyse_scenarios)
}
 