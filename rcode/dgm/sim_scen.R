#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Simulation scenarios
#------------------------------------------------------------------------------#


# Helper function ----
#------------------------------------------------------------------------------#
# Helper functions to remove scenarios that only differ with respect to ordering
# of L effects and similar otherwise
# 
# Desired effect combinations and notation thereof:
#	1:(bAL==0 & bYL==0) , 2:(bAL==0 & bYL==0.2) , 3:(bAL==0 & bYL==0.4) + 
#	4:(bAL==0.2 & bYL==0) , 5:(bAL==0.2 & bYL==0.2) , 6:(bAL==0.2 & bYL==0.4) + 
#	7:(bAL==0.4 & bYL==0) , 8:(bAL==0.4 & bYL==0.2) , 9:(bAL==0.4 & bYL==0.4)} 

# Take exponent in computing rowtotals to avoid accidental doubles
add_counter     <- function(x) sum(exp(x)) 
fun_patterns_AL <- function(z) 1*(z==1 | z==2 | z==3) + 2*(z==4 | z==5 |z==6) + 3*(z==7 | z==8 |z==9)
fun_patterns_YL <- function(z) 1*(z==1 | z==4 | z==7) + 2*(z==2 | z==5 |z==8) + 3*(z==3 | z==6 |z==9)


# Initialize simulation scenarios ----
# Function output: data.frame with 3960 scenarios for the dgm, varying the 
# strength of association between covariates and exposure/outcome, residual variance
# of the outcome, correlation between covariates, effective sample size and
# event rate.
#------------------------------------------------------------------------------#


datagen_scenarios <- function(){
  nevents   <- c(50, 200)
  nL        <- 24
  bYA       <- c(0, log(1.5))
  bL_const  <- log(1.05)
  bAL       <- c(0, log(1.05), log(1.2))
  bYL       <- c(0, log(1.05), log(1.2))
  eventrate <- c(0.2, 0.03)
  Yint      <- c(-2.1, -4.05)
  rhoL      <- 0.3
  
  # data.frame with simulation scenarios
  variable_grid <- expand.grid(nevents=nevents,
                               nL=nL,
                               bYA=bYA,
                               bL_const=bL_const,
                               eventrate=eventrate,
                               rhoL=rhoL)
  
  # 9 possible covariates
  cov1 <- cov2  <- cov3  <- cov4  <- 1:9
  # Identify duplicate combinations (differ only with respect to order)
  covariate_grid <- expand.grid(cov1 = cov1,
                                cov2 = cov2,
                                cov3 = cov3,
                                cov4 = cov4)
  covariate_grid$counter    <- apply(covariate_grid,1,add_counter) 
  covariate_grid$duplicates <- 0 
  for (i in 2:nrow(covariate_grid)){
    covariate_grid$duplicates[i] <- 1*(covariate_grid$counter[i] %in% 	covariate_grid$counter[1:(i-1)])		
  }
  # Remove duplicates
  covariate_select <- subset(covariate_grid,subset=covariate_grid$duplicates==0)
  
  datagen_scenarios <- data.frame(nevents = rep(variable_grid$nevents,each=nrow(covariate_select)),
                                  nL = rep(variable_grid$nL,each=nrow(covariate_select)),
                                  bYA = rep(variable_grid$bYA,each=nrow(covariate_select)),
                                  bL_const = rep(variable_grid$bL_const,each=nrow(covariate_select)),
                                  eventrate = rep(variable_grid$eventrate,each=nrow(covariate_select)),
                                  rhoL = rep(variable_grid$rhoL,each=nrow(covariate_select)),
                                  
                                  cov1 = rep(covariate_select$cov1,nrow(variable_grid)),
                                  cov2 = rep(covariate_select$cov2,nrow(variable_grid)),
                                  cov3 = rep(covariate_select$cov3,nrow(variable_grid)),
                                  cov4 = rep(covariate_select$cov4,nrow(variable_grid)))	
  
  datagen_scenarios$bAL1 <- bAL[fun_patterns_AL(datagen_scenarios$cov1)]
  datagen_scenarios$bYL1 <- bYL[fun_patterns_YL(datagen_scenarios$cov1)]
  datagen_scenarios$bAL2 <- bAL[fun_patterns_AL(datagen_scenarios$cov2)]
  datagen_scenarios$bYL2 <- bYL[fun_patterns_YL(datagen_scenarios$cov2)]
  datagen_scenarios$bAL3 <- bAL[fun_patterns_AL(datagen_scenarios$cov3)]
  datagen_scenarios$bYL3 <- bYL[fun_patterns_YL(datagen_scenarios$cov3)]
  datagen_scenarios$bAL4 <- bAL[fun_patterns_AL(datagen_scenarios$cov4)]
  datagen_scenarios$bYL4 <- bYL[fun_patterns_YL(datagen_scenarios$cov4)]
  
  datagen_scenarios$cov1 <- 
    datagen_scenarios$cov2 <- 
    datagen_scenarios$cov3 <- 
    datagen_scenarios$cov4 <- NULL
  
  # Read in intercepts (estimated based on large sample approximation of event rate
  # in ./rcode/dgm/estimate_intercepts.R)
  datagen_scenarios$Yint <- c(readRDS("./rcode/dgm/intercepts.rds"))
  
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
