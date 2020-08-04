#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Proof of concept: when to use selection of potential confounders
#------------------------------------------------------------------------------#

# Running this proof-of-concept simulation takes around 6.5 hours


# Description central principle:
# Dgm: binary exposure (A), single continuous covariate (L), binary outcome (Y).
# Two strategies: always include or always omit covariable, L.
# Assumption: bias_full = 0.
# Then, MSE_omit < MSE_full, when: bias_omit^2 < var_full - var_omit



# Load libraries  ----
#------------------------------------------------------------------------------#
library(logistf)

# Simulation  ----
#------------------------------------------------------------------------------#
# Parameters 
nsim      <- 10000
nobs      <- c(60,120,300)
Yint      <- c(0,-1.65)   # Control eventrate by setting intercept 0 or -1.65 on log-scale
bLA       <- seq(0,0.5,by=0.1)
bLY       <- seq(0,0.5,by=0.1)
truth     <- 1

# define scenarios
scen         <- expand.grid(nobs = nobs,
                              Yint = Yint, 
                              bLA = bLA,
                              bLY = bLY)

# Output matrices 
Out_OR_Omit  <- matrix(NA, nrow=nrow(scen), ncol = nsim)
Out_OR_Full  <- matrix(NA, nrow=nrow(scen), ncol = nsim)
Out_MRR_Omit <- matrix(NA, nrow=nrow(scen), ncol = nsim)
Out_MRR_Full <- matrix(NA, nrow=nrow(scen), ncol = nsim)

set.seed(20200429) 
seeds <- matrix(sample(1:1000000000, size = nsim*nrow(scen), replace = T),
                nrow=nrow(scen), ncol = nsim)

# Run basic sims
for(j in 1:nrow(scen)){
  for(i in 1:nsim){
    
    bLA <- scen$bLA[j]
    bLY <- scen$bLY[j]
    nobs<- scen$nobs[j]
    Yint<- scen$Yint[j]
    
    set.seed(seeds[j,i])
    
    # Data generating mechanism: generate under null
    L  <- rnorm(nobs)
    A  <- rbinom(nobs,1,plogis(bLA * L))
    Y  <- rbinom(nobs,1,plogis(Yint + bLY * L))
    df <- data.frame(Y,A,L)
    
    # Estimated MRR
    
    ## Covariate L Omitted
    M_Omit     <- logistf(Y~A, pl = F)                                       # Estimate logistic regression model using 
    # Firth's correction
    Out_OR_Omit[j,i]  <- M_Omit$coefficients["A"]
    M_pred_Omit       <- df$A * M_Omit$coefficients[-1]                             # Re-estimate the intercept
    M_Omit_int        <- glm(df$Y ~ offset(M_pred_Omit), family = binomial)         # Re-estimate the intercept
    
    
    p0_Omit           <- plogis(M_Omit_int$coefficients[1] + 0 * M_Omit$coefficients[-1])
    p1_Omit           <- plogis(M_Omit_int$coefficients[1] + 1 * M_Omit$coefficients[-1])
    Out_MRR_Omit[j,i] <- sum(p1_Omit)/sum(p0_Omit)
    
    
    ## Covariate L Forced
    M_Full            <- logistf(Y~A+L, pl = F)                                     # Estimate logistic regression model using
    # Firth's correction
    Out_OR_Full[j,i]  <- M_Full$coefficients["A"]
    M_pred_Full       <- as.matrix(df[,names(df)!="Y"]) %*% M_Full$coefficients[-1] # Re-estimate the intercept
    M_Full_int        <- glm(df$Y ~ offset(M_pred_Full), family = binomial)         # Re-estimate the intercept
    
    p0_Full           <- plogis(M_Full_int$coefficients[1] + cbind(0,L) %*% matrix(M_Full$coefficients[-1],ncol=1))
    p1_Full           <- plogis(M_Full_int$coefficients[1] + cbind(1,L) %*% matrix(M_Full$coefficients[-1],ncol=1))
    Out_MRR_Full[j,i] <- sum(p1_Full)/sum(p0_Full)
    
  }
}

# Save raw output
saveRDS(list(Out_OR_Omit = Out_OR_Omit, Out_OR_Full = Out_OR_Full,
             Out_MRR_Omit = Out_MRR_Omit, Out_MRR_Full = Out_MRR_Full), 
        file = "./data/poc/CICovSel_POC_raw.rds")

# Output
Bias2_Omit_OR <- (rowMeans(Out_MRR_Omit) - log(truth))^2
Var_Omit_OR    <- rowMeans(
  sweep(Out_MRR_Omit, 
        MARGIN = 1, 
        STATS = rowMeans(Out_MRR_Omit))
  ^2)

Bias2_Full_OR <- (rowMeans(Out_MRR_Full) - log(truth))^2
Var_Full_OR    <- rowMeans(
  sweep(Out_MRR_Full, 
        MARGIN = 1, 
        STATS = rowMeans(Out_MRR_Full))
  ^2)

Bias2_Omit_MRR <- (rowMeans(log(Out_MRR_Omit)) - log(truth))^2
Var_Omit_MRR    <- rowMeans(
  sweep(log(Out_MRR_Omit), 
        MARGIN = 1, 
        STATS = rowMeans(log(Out_MRR_Omit)))
  ^2)

Bias2_Full_MRR <- (rowMeans(log(Out_MRR_Full)) - log(truth))^2
Var_Full_MRR    <- rowMeans(
  sweep(log(Out_MRR_Full), 
        MARGIN = 1, 
        STATS = rowMeans(log(Out_MRR_Full)))
  ^2)

# Comparison
Comparison <- scen
Comparison$Bias2_omit_MRR <- Bias2_Omit_MRR
Comparison$Var_MRR        <- Var_Full_MRR - Var_Omit_MRR
Comparison$Bias2_full_MRR <- Bias2_Full_MRR
Comparison$Bias2_omit_OR <- Bias2_Omit_OR
Comparison$Var_OR        <- Var_Full_OR - Var_Omit_OR
Comparison$Bias2_full_OR <- Bias2_Full_OR

saveRDS(Comparison, 
        file = "./data/poc/CICovSel_POC.rds")
