#------------------------------------------------------------------------------#
# Large sample simulation for intercept estimation
# This script is for documentation
# Results are incorporated in ./rcode/dgm/sim_scen.R
#------------------------------------------------------------------------------#

# Take exponent in computing rowtotals to avoid accidental doubles
add_counter     <- function(x) sum(exp(x)) 
fun_patterns_AL <- function(z) 1*(z==1 | z==2 | z==3) + 2*(z==4 | z==5 |z==6) + 3*(z==7 | z==8 |z==9)
fun_patterns_YL <- function(z) 1*(z==1 | z==4 | z==7) + 2*(z==2 | z==5 |z==8) + 3*(z==3 | z==6 |z==9)


# Initialize simulation scenarios ----
# Function output: data.frame with 7920 scenarios for the dgm, varying the 
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
  
  # Add scenarios number to keep track of results
  datagen_scenarios$scen_num <- c(1:nrow(datagen_scenarios))
  
  return(datagen_scenarios)
}


gen_data <- function(nobs,
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
                     rhoL){
  # Generate Y, A and L
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
  
  # Function to optimize for intercept
  Y <- rbinom(nobs, 1, plogis(Yint + 
                                 A * bYA + 
                                 L %*% betasYL))
  
  prev_diff <- eventrate - mean(Y)
  return(prev_diff)
}

scenarios <- datagen_scenarios()
Yint <- matrix(NA, nrow = nrow(scenarios), ncol = 1)

# Bisection method
for(i in 1:nrow(scenarios)){
  a <- 0
  b <- -6
  c <-(a+b)/2
  current_value <- gen_data(nobs = 100000,
                            nL = scenarios[i,"nL"],
                            bYA = scenarios[i,"bYA"],
                            bL_const = scenarios[i,"bL_const"],
                            bAL1 = scenarios[i,"bAL1"],
                            bAL2 = scenarios[i,"bAL2"],
                            bAL3 = scenarios[i,"bAL3"],
                            bAL4 = scenarios[i,"bAL4"],
                            bYL1 = scenarios[i,"bYL1"],
                            bYL2 = scenarios[i,"bYL2"],
                            bYL3 = scenarios[i,"bYL3"],
                            bYL4 = scenarios[i,"bYL4"],
                            eventrate = scenarios[i,"eventrate"],
                            rhoL = scenarios[i,"rhoL"],
                            Yint = c)
  while(current_value != 0 && abs(b - a) >.002){
    if(current_value < 0)
    {
      a <- c
    }
    else{
      b <- c
    }
    c <- (a + b)/ 2
    
    current_value <- gen_data(nobs = 100000,
                              nL = scenarios[i,"nL"],
                              bYA = scenarios[i,"bYA"],
                              bL_const = scenarios[i,"bL_const"],
                              bAL1 = scenarios[i,"bAL1"],
                              bAL2 = scenarios[i,"bAL2"],
                              bAL3 = scenarios[i,"bAL3"],
                              bAL4 = scenarios[i,"bAL4"],
                              bYL1 = scenarios[i,"bYL1"],
                              bYL2 = scenarios[i,"bYL2"],
                              bYL3 = scenarios[i,"bYL3"],
                              bYL4 = scenarios[i,"bYL4"],
                              eventrate = scenarios[i,"eventrate"],
                              rhoL = scenarios[i,"rhoL"],
                              Yint = c)
    
    c <- (a + b)/ 2;c
  }
  
  Yint[i,] <- c
  
}

saveRDS(Yint,"./rcode/dgm/intercepts.rds")
