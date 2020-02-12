##############################################
# CovSel simulations
#
# December 2019
# Author: K Luijken
##############################################

# This directory is used to save and load documents
filepath <- "//vf-d2-home/d2home$/kluijken/MyDocs/Projects/CovSelection/Simulations/CI_CovSel/"
# Note that simulation output is stored in files filepath/Simoutput/S1...S3
# Make sure filepath ends with "/


# Source packages and functions
source("packages.R")
source("subfunctions.R")
source("RUN.R")

#-----------------------------------------------------------------------#
# Simulation parameters
#-----------------------------------------------------------------------#



# Determine which models should be estimated
Allmodels    <- c("Full", "Selectedp0.157", "Selectedp0.320")
Allmarginals <- c("Truth","Unadjusted","Full", "Selectedp0.157", "Selectedp0.320")


# Set dgm parameters
intA    <- 0                    # Intercept for generating exposure
setbL1Y <- 1                    # Static effect of covariate(s) L1 on exposure A (unconditional, later on standardized over number of covariates in L1)
setbL1A <- 1                    # Static effect of covariate(s) L1 on outcome Y (unconditional, later on standardized over number of covariates in L1)

# Simulation parameters
nsim   <- 25                    # Number of times the simulations are repeated
B      <- 100                   # Number of bootstrap samples
rhoL   <- seq(0,0.6, by=0.3)    # Mutual correlation of covariates, L
nobs   <- seq(100, 300, by=100) # Number of observations in one sample
nL     <- c(2,6,10)             # Number of covariates, L
bAY    <- c(0,1)                # Effect of exposure A on outcome Y (unconditional)
S1bL2A <- seq(0,0.5,by=0.1)     # Scenario 1: Effect of covariate(s) L2 on exposure A (unconditional)
S1bL2Y <- seq(0,0.5,by=0.1)     # Scenario 1: Effect of covariate(s) L2 on outcome Y (unconditional)
S2bL2A <- seq(0.5,1,by=0.1)     # Scenario 2:  Effect of covariate(s) L2 on exposure A (unconditional)
S2bL2Y <- seq(0,0.1,by=0.1)     # Scenario 2: Effect of covariate(s) L2 on outcome Y (unconditional)
S3bL2A <- seq(0,0.1,by=0.1)     # Scenario 3:  Effect of covariate(s) L2 on exposure A (unconditional)
S3bL2Y <- seq(0.5,1,by=0.1)     # Scenario 3: Effect of covariate(s) L2 on outcome Y (unconditional)


#-----------------------------------------------------------------------#
# Scenario 1
#-----------------------------------------------------------------------#

# S1, simulation matrix
S1 <- expand.grid(S1bL2Y, rhoL, nobs, nL, bAY)
S1 <- cbind(S1[,1],S1)
colnames(S1) <- c("bL2A","bL2Y", "rhoL", "nobs", "nL", "bAY")
S1[,"bL2A"] <- S1[,"bL2A"]/(sqrt(S1[,"nL"]/2))
S1[,"bL2Y"] <- S1[,"bL2Y"]/(sqrt(S1[,"nL"]/2))
S1 <- data.matrix(S1)

#-----------------------------------------------------------------------#
# Scenario 2
#-----------------------------------------------------------------------#

# S2, simulation matrix
S2 <- expand.grid(S2bL2A,S2bL2Y, rhoL, nobs, nL, bAY)
colnames(S2) <- c("bL2A","bL2Y", "rhoL", "nobs", "nL", "bAY")
S2[,"bL2A"] <- S2[,"bL2A"]/(sqrt(S2[,"nL"]/2))
S2[,"bL2Y"] <- S2[,"bL2Y"]/(sqrt(S2[,"nL"]/2))
S2 <- data.matrix(S2)


#-----------------------------------------------------------------------#
# Scenario 3
#-----------------------------------------------------------------------#

# S3, simulation matrix
S3 <- expand.grid(S3bL2A,S3bL2Y, rhoL, nobs, nL, bAY)
colnames(S3) <- c("bL2A","bL2Y", "rhoL", "nobs", "nL", "bAY")
S3[,"bL2A"] <- S3[,"bL2A"]/(sqrt(S3[,"nL"]/2))
S3[,"bL2Y"] <- S3[,"bL2Y"]/(sqrt(S3[,"nL"]/2))
S3 <- data.matrix(S3)

#-----------------------------------------------------------------------#
# Save setup
#-----------------------------------------------------------------------#

Setup <- list(
    filepath = filepath,
    Allmodels = Allmodels,
    Allmarginals = Allmarginals,
    intA = intA,
    setbL1Y = setbL1Y,
    setbL1A = setbL1A,
    nsim   = nsim,
    B      = B,
    rhoL   = rhoL,
    nobs   = nobs,
    nL     = nL,
    bAY    = bAY,
    S1bL2A = S1bL2A,
    S1bL2Y = S1bL2Y,
    S2bL2A = S2bL2A,
    S2bL2Y = S2bL2Y,
    S3bL2A = S3bL2A,
    S3bL2Y = S3bL2Y,
    S1     = S1,
    S2     = S2,
    S3     = S3)

saveRDS(Setup,paste0(filepath,"PlotsTables/Setup.rds"))

#-----------------------------------------------------------------------#
# Run simulation, output stored in RDS-files
#-----------------------------------------------------------------------#

RUN(scenario="S1",
    nsim=nsim,
    B=B,
    Allmodels=Allmodels,
    Allmarginals=Allmarginals)

RUN(scenario="S2",
    nsim=nsim,
    B=B,
    Allmodels=Allmodels,
    Allmarginals=Allmarginals)

RUN(scenario="S3",
    nsim=nsim,
    B=B,
    Allmodels=Allmodels,
    Allmarginals=Allmarginals)


