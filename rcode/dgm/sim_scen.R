#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Simulation scenarios
#------------------------------------------------------------------------------#


# Helper function ----
#------------------------------------------------------------------------------#
compute_rowtotal <- function(x){sum(x, combn(x, m = 2, FUN = prod))}

# Initialize simulation scenarios ----
# Function output: data.frame with 8748 scenarios for the dgm, varying the 
# strength of association between covariates and exposure/outcome, residual variance
# of exposure and outcome, correlation between covariates, sample size and
# number of covariates.
#------------------------------------------------------------------------------#


datagen_scenarios <- function(){
  nobs  <- c(120,300)
  nL    <- 9
  bAY   <- c(0,0.5)
  bAL1  <- bAL2  <- bAL3  <- c(0,0.5,1)
  bYL1  <- bYL2  <- bYL3  <- c(0,0.5,1)
  Yint  <- 0
  sd_UY <- c(0.01,1)
  rhoL  <- c(0,0.3,0.7)
  
  # data.frame with simulation scenarios
  datagen_scenarios <- expand.grid(nobs=nobs,
                                   nL=nL,
                                   bAY=bAY,
                                   bAL1=bAL1,
                                   bAL2=bAL2,
                                   bAL3=bAL3,
                                   bYL1=bYL1,
                                   bYL2=bYL2,
                                   bYL3=bYL3,
                                   Yint=Yint,
                                   sd_UY=sd_UY,
                                   rhoL=rhoL)
  
  # Remove redundant scenarios using a counter
  datagen_scenarios$counter <- apply(datagen_scenarios,
                                     MARGIN=1,
                                     FUN = function(x) compute_rowtotal(x))
  datagen_scenarios         <- datagen_scenarios[!duplicated(
                                     datagen_scenarios[['counter']]),]
  datagen_scenarios$counter <- NULL
  
  # Add scenarios number to keep track of results
  datagen_scenarios$scen_num <- c(1:nrow(datagen_scenarios))

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
 