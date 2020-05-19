#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Proof of concept: when to use selection of potential confounders
#------------------------------------------------------------------------------#


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
nsim      <- 1
nobs      <- c(60,120,300)
Yint      <- c(0,-1.65)   # Control eventrate by setting intercept 0 or -1.6 on log-scale
bLA       <- seq(0,0.5,by=0.1)
bLY       <- seq(0,0.5,by=0.1)
MRR.truth <- 1

scen <- expand.grid(nobs, Yint, bLA, bLY)
colnames(scen) <-c("nobs","Yint","bLA","bLY")

# Output matrices 
Out.MRR.Omit <- matrix(NA, nrow=nrow(scen),ncol = nsim)
Out.MRR.Full <- matrix(NA, nrow=nrow(scen),ncol = nsim)

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
    M.Omit     <- logistf(Y~A, pl = F)                                       # Estimate logistic regression model using 
    # Firth's correction
    M.pred.Omit<- df$A * M.Omit$coefficients[-1]                             # Re-estimate the intercept
    M.Omit.int <- glm(df$Y ~ offset(M.pred.Omit), family = binomial)         # Re-estimate the intercept
    
    
    p0.Omit           <- plogis(M.Omit.int$coefficients[1] + 0 * M.Omit$coefficients[-1])
    p1.Omit           <- plogis(M.Omit.int$coefficients[1] + 1 * M.Omit$coefficients[-1])
    Out.MRR.Omit[j,i] <- sum(p1.Omit)/sum(p0.Omit)
    
    
    ## Covariate L Forced
    M.Full     <- logistf(Y~A+L, pl = F)                                     # Estimate logistic regression model using
    # Firth's correction
    M.pred.Full<- as.matrix(df[,names(df)!="Y"]) %*% M.Full$coefficients[-1] # Re-estimate the intercept
    M.Full.int <- glm(df$Y ~ offset(M.pred.Full), family = binomial)         # Re-estimate the intercept
    
    p0.Full           <- plogis(M.Full.int$coefficients[1] + cbind(0,L) %*% matrix(M.Full$coefficients[-1],ncol=1))
    p1.Full           <- plogis(M.Full.int$coefficients[1] + cbind(1,L) %*% matrix(M.Full$coefficients[-1],ncol=1))
    Out.MRR.Full[j,i] <- sum(p1.Full)/sum(p0.Full)
    
  }
}

# Output
Bias2.Omit <- (rowMeans(log(Out.MRR.Omit)) - log(MRR.truth))^2
Var.Omit    <- rowMeans(
  sweep(log(Out.MRR.Omit), 
        MARGIN = 1, 
        STATS = rowMeans(log(Out.MRR.Omit)))
  ^2)

Bias2.Full <- (rowMeans(log(Out.MRR.Full)) - log(MRR.truth))^2
Var.Full    <- rowMeans(
  sweep(log(Out.MRR.Full), 
        MARGIN = 1, 
        STATS = rowMeans(log(Out.MRR.Full)))
  ^2)

# Comparison
Comparison <- scen
Comparison$Bias2_omit <- Bias2.Omit
Comparison$Var        <- Var.Full - Var.Omit
Comparison$Bias2_full <- Bias2.Full

saveRDS(Comparison, 
        file = "./rcode/poc/CICovSel_POC.rds")
