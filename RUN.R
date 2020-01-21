####################################################
# Covariate selection project
#
# December 2019 
# Author: K Luijken
# File content: run simulations
####################################################



RUN <- function(scenario, nsim,B, Allmodels, Allmarginals){
  # Arguments:
  #    scenario     = "S1", "S2" or "S3" (character) referring to weak confounder, near-instrument or predictor
  #    nsim         = number of times the simulation is repeated
  #    B            = number or bootstrap samples
  #    Allmodels    = models used, c("Full", "Selectedp0.157", "Selectedp0.320")
  #    Allmarginals = marginals computed or estimated, c("Truth","Unadjusted","Full", "Selectedp0.157", "Selectedp0.320")
  
  S <- get(scenario)
  
  for(j in 1:nrow(S)){
    progress(j, max.value = nrow(S))
    
    # Store output exposure/outcome frequency 
    Setup         <- matrix(NA,
                            nrow=2,ncol=nsim,
                            dimnames = list(c("Pexposure",
                                              "Poutcome"),
                                            NULL))
    
    # Store model and marginal outcome
    Model         <- array(NA,
                           dim = c(nrow=36,ncol=nsim,height=length(Allmodels)),
                           dimnames = list(c("Int", 
                                             "bAY",
                                             paste0("bL1",seq(1:(max(nL)/2))), 
                                             paste0("bL2",seq(1:(max(nL)/2))), 
                                             "se(Int)",
                                             "se(bAY)",
                                             paste0("se(bL1)",seq(1:(max(nL)/2))), 
                                             paste0("se(bL2)",seq(1:(max(nL)/2))),
                                             "RMSD bAY",
                                             "RMSD bAY ratio",
                                             paste0("incl. freq. bL1",seq(1:(max(nL)/2))), 
                                             paste0("incl. freq. bL2",seq(1:(max(nL)/2)))),
                                           NULL,
                                           Allmodels))
    Marginal      <- array(NA,
                           dim=c(nrow=10, ncol= nsim,heigth=length(Allmarginals)),
                           dimnames = list(c("MRR",
                                             "log(MRR)",
                                             "error",
                                             "log(error)",
                                             "sq.error", 
                                             "log(sq.error)",
                                             "MOR", 
                                             "cilow", 
                                             "ciup",
                                             "coverage"),
                                           NULL,
                                           Allmarginals))
    
    # Storage warnings in lists
    Warning      <- list(Full = list (),
                    Selectedp0.157 = list(),
                    Selectedp0.320 = list())
    
    Warning.int  <- list(Full = list (),
                    Selectedp0.157 = list(),
                    Selectedp0.320 = list())
    
    # Storage Bootstrap warnings
    BSWarning    <- list(Full = rep(list(NA),times=nsim*B), 
                      Selectedp0.157 = rep(list(NA),times=nsim*B),
                      Selectedp0.320 = rep(list(NA),times=nsim*B))
    
    BSWarning.int<- list(Full = rep(list(NA),times=nsim*B), 
                        Selectedp0.157 = rep(list(NA),times=nsim*B),
                        Selectedp0.320 = rep(list(NA),times=nsim*B))
    
    
    # Set the parameters for bL1Y, bL1A
    bL1Y <- setbL1Y/(sqrt(S[j,"nL"]/2))
    bL1A <- setbL1A/(sqrt(S[j,"nL"]/2))
    
    #-----------------------------------------------------------------------#
    # Truth
    #-----------------------------------------------------------------------#
    
    # Compute the true marginal risk ratio in the dgm using numerical integration
    Truthoutput <- Truth(S=S, j=j, bL1Y = bL1Y, bL1A=bL1A)
    
    Marginal[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"), ,
             "Truth"] <- as.numeric(Marginalcomp(pY1 = Truthoutput$pY1, # User-defined function Marginalcomp
                                                      pY0 = Truthoutput$pY0,
                                                      truth = sum(Truthoutput$pY1)/sum(Truthoutput$pY0)))
    
    
    for(i in 1:nsim){
      # generate A, Y and set of L
      set.seed(readRDS(paste0(filepath,"Siminput/",scenario,"/seeds",scenario,"_",j,".rds"))[i])
      sigma                  <- matrix(S[j,"rhoL"], nrow=S[j,"nL"],ncol=S[j,"nL"])
      diag(sigma)            <- 1
      L                      <- mvrnorm(S[j,"nobs"],c(rep(0,times=S[j,"nL"])), sigma)
      L1                     <- L[,1:(S[j,"nL"]/2)]
      L2                     <- L[,((S[j,"nL"]/2)+1):S[j,"nL"]]
      A                      <- rbinom(S[j,"nobs"],1,plogis(intA + as.matrix(L1)%*%rep(bL1A,times=(S[j,"nL"]/2)) + as.matrix(L2)%*%rep(S[j,"bL2A"],times=(S[j,"nL"]/2))))
      Setup["Pexposure",i]   <- sum(A)/S[j,"nobs"]
      Y                      <- rbinom(S[j,"nobs"],1,plogis(intY + S[j,"bAY"]*A + as.matrix(L1)%*%rep(bL1Y,times=(S[j,"nL"]/2)) + as.matrix(L2)%*%rep(S[j,"bL2A"],times=(S[j,"nL"]/2))))
      Setup["Poutcome",i]    <- sum(Y)/S[j,"nobs"]
      
      # Store data such that it is available within the function
      assign("datafull", 
             data.frame(matrix(c(Y,A,L1,L2),
                               nrow=S[j,"nobs"], 
                               ncol=2+S[j,"nL"], 
                               dimnames = list(NULL,c("Y","A",paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2)))))), 
             envir = parent.frame())
      
      #-----------------------------------------------------------------------#
      # Unadjusted
      #-----------------------------------------------------------------------#
      
      # Compute unadjusted marginal risk ratio (MRR, log(MRR), squared error, log(SqE))
      Marginal[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
               i,
               "Unadjusted"] <- as.numeric(Marginalcomp(pY1=with(datafull,mean(Y[A==1])), # User-defined function Marginalcomp
                                                             pY0=with(datafull,mean(Y[A==0])),
                                                             truth = Marginal["MRR",i,"Truth"]))
      
      #-----------------------------------------------------------------------#
      # Full model
      #-----------------------------------------------------------------------#
      # Estimate a model using all L (full model), using Firth type Logistic regression with Intercept Correction (FLIC)
      Fulloutput <- Full(S=S,
                         j=j,
                         datafull=datafull,
                         i=i)
      
      # Store full model output
      Warning$Full[[i]]           <- Fulloutput$warning
      Warning.int$Full[[i]]       <- Fulloutput$warning.int
      
      Model[1:24,i,"Full"]        <- Modeloutput(intmodel = Fulloutput$M.int, # User-defined function Modeloutput
                                                 model = Fulloutput$M)
      Marginal[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
              i,
              "Full"]             <- as.numeric(Marginalest(data = datafull,  # User-defined function Marginalest
                                                 int = Fulloutput$M.int$coefficients["(Intercept)"],
                                                 modelcoefs = Fulloutput$M$coefficients[-1],
                                                 truth = Marginal["MRR",i,"Truth"]))
      
      
      #-----------------------------------------------------------------------#
      # Backward elimination model, p = 0.157
      #-----------------------------------------------------------------------#
      
      # Estimate a model using backward elimination (selected model), using FLIC
      Selectedoutput0157 <- Selected(Mod=Fulloutput$M,
                                     pcutoff = 0.157,
                                     S=S,
                                     j=j,
                                     datafull=datafull,
                                     i=i)
      
      
      # Store selected model output
      Warning$Selectedp0.157[[i]]    <- Selectedoutput0157$warning
      Warning.int$Selectedp0.157[[i]]<- Selectedoutput0157$warning.int
      
      Model[1:24,i,"Selectedp0.157"] <- Modeloutput(intmodel = Selectedoutput0157$Msel.int,  # User-defined function Modeloutput
                                                                  model = Selectedoutput0157$Msel)
    
      Marginal[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
               i,
               "Selectedp0.157"] <- as.numeric(Marginalest(data = datafull,              # User-defined function Marginalest
                                               int = Selectedoutput0157$Msel.int$coefficients["(Intercept)"],
                                               modelcoefs = Selectedoutput0157$Msel$coefficients[-1],
                                               truth = Marginal["MRR",i,"Truth"]))
      #-----------------------------------------------------------------------#
      # Backward elimination model, p = 0.320
      #-----------------------------------------------------------------------#
      
      # Estimate a model using backward elimination (selected model), using FLIC
      Selectedoutput0320 <- Selected(Mod = Fulloutput$M,
                                     pcutoff = 0.320,
                                     S=S,
                                     j=j,
                                     datafull=datafull,
                                     i=i)

      # Store selected model output
      Warning$Selectedp0.320[[i]]     <- Selectedoutput0320$warning
      Warning.int$Selectedp0.320[[i]] <- Selectedoutput0320$warning.int
      
      Model[1:24,i,"Selectedp0.320"]  <- Modeloutput(intmodel = Selectedoutput0320$Msel.int,  # User-defined function Modeloutput
                                                model = Selectedoutput0320$Msel)
      
      Marginal[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
               i,
               "Selectedp0.320"] <- as.numeric(Marginalest(data = datafull,              # User-defined function Marginalest
                                                           int = Selectedoutput0320$Msel.int$coefficients["(Intercept)"],
                                                           modelcoefs = Selectedoutput0320$Msel$coefficients[-1],
                                                           truth = Marginal["MRR",i,"Truth"]))
      
      # Bootstrap the MRRs to obtain confidence intervals and coverage
      # Create Bootstrap storage
      bsoutput  <- array(NA,
                             dim=c(nrow=9, ncol=B,heigth=length(Allmodels)),
                             dimnames = list(c("MRR",
                                               "log(MRR)",
                                               "error",
                                               "log(error)",
                                               "sq.error", 
                                               "log(sq.error)",
                                               "MOR", 
                                               "bAYboot",
                                               "sq.error(bAYboot)"),
                                             NULL,
                                             Allmodels))
      bsest    <- array(NA,
                         dim=c(nrow=12, ncol=B,heigth=length(Allmodels)),
                         dimnames = list(c("Int", 
                                           "bAY",
                                           paste0("bL1",seq(1:(max(nL)/2))), 
                                           paste0("bL2",seq(1:(max(nL)/2)))),
                                         NULL,
                                         Allmodels))
      
      # # Resample
      for(k in 1:B){
         set.seed(readRDS(paste0(filepath,"Siminput/",scenario,"/bsseeds",scenario,"_",j,".rds"))[k,i])
         bs <- sample(1:S[j,"nobs"],size=S[j,"nobs"],replace=T)
         assign("databs", data.frame(matrix(c(Y,A,L1,L2),nrow=S[j,"nobs"], ncol=2+S[j,"nL"], dimnames = list(NULL,c("Y","A",paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2))))))[bs,], envir = parent.frame())
         
         # Full model
         bsFulloutput <- Fullbs(S=S,
                            j=j,
                            databs=databs,
                            i=i,
                            k=k)
         
         # Store full model output
         BSWarning$Full[[(i-1)*B+k]]     <- bsFulloutput$warning
         BSWarning.int$Full[[(i-1)*B+k]] <- bsFulloutput$warning.int
         
         selestfull <- coef(bsFulloutput$M)[c("(Intercept)", 
                                              "A",
                                              paste0("L1",seq(1:(max(nL)/2))), 
                                              paste0("L2",seq(1:(max(nL)/2))))]
         bsest[,k,"Full"] <- selestfull
         bsoutput["bAYboot",k,"Full"] <- bsFulloutput$M$coefficients["A"]
         bsoutput["sq.error(bAYboot)",k,"Full"] <- (bsFulloutput$M$coefficients["A"] - Model["bAY",i,"Full"])^2
         bsoutput[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
                  k,
                  "Full"]             <- as.numeric(Marginalest(data = databs,  # User-defined function Marginalest
                                                                int = bsFulloutput$M.int$coefficients["(Intercept)"],
                                                                modelcoefs = bsFulloutput$M$coefficients[-1],
                                                                truth = Marginal["MRR",i,"Truth"]))
         # Selected model, p = 0.157
         bsSelectedoutput0157 <- Selectedbs(Mod=bsFulloutput$M,
                                        pcutoff = 0.157,
                                        S=S,
                                        j=j,
                                        databs = databs,
                                        i=i,
                                        k=k)
         
         
         # Store selected model output
         BSWarning$Selectedp0.157[[(i-1)*B+k]]     <- bsSelectedoutput0157$warning
         BSWarning.int$Selectedp0.157[[(i-1)*B+k]] <- bsSelectedoutput0157$warning.int
         
         selestsel0157 <- coef(bsSelectedoutput0157$Msel)[c("(Intercept)", 
                                                            "A",
                                                            paste0("L1",seq(1:(max(nL)/2))), 
                                                            paste0("L2",seq(1:(max(nL)/2))))]
         bsest[,k,"Selectedp0.157"] <- selestsel0157
         bsoutput["bAYboot",k,"Selectedp0.157"] <- bsSelectedoutput0157$Msel$coefficients["A"]
         bsoutput["sq.error(bAYboot)",k,"Selectedp0.157"] <- (bsSelectedoutput0157$Msel$coefficients["A"] - Model["bAY",i,"Selectedp0.157"])^2
         bsoutput[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
                  k,
                  "Selectedp0.157"] <- as.numeric(Marginalest(data = databs,              # User-defined function Marginalest
                                                              int = bsSelectedoutput0157$Msel.int$coefficients["(Intercept)"],
                                                              modelcoefs = bsSelectedoutput0157$Msel$coefficients[-1],
                                                              truth = Marginal["MRR",i,"Truth"]))
         
         # Selected model, p = 0.320
         bsSelectedoutput0320 <- Selectedbs(Mod=bsFulloutput$M,
                                          pcutoff = 0.320,
                                          S=S,
                                          j=j,
                                          databs=databs,
                                          i=i,
                                          k=k)
         
         
         # Store selected model output
         BSWarning$Selectedp0.320[[(i-1)*B+k]]     <- bsSelectedoutput0320$warning
         BSWarning.int$Selectedp0.320[[(i-1)*B+k]] <- bsSelectedoutput0320$warning.int
         
         selestsel0320 <- coef(bsSelectedoutput0320$Msel)[c("(Intercept)", 
                                                            "A",
                                                            paste0("L1",seq(1:(max(nL)/2))), 
                                                            paste0("L2",seq(1:(max(nL)/2))))]
         bsest[,k,"Selectedp0.320"] <- selestsel0320
         bsoutput["bAYboot",k,"Selectedp0.320"] <- bsSelectedoutput0320$Msel$coefficients["A"]
         bsoutput["sq.error(bAYboot)",k,"Selectedp0.320"] <- (bsSelectedoutput0320$Msel$coefficients["A"] - Model["bAY",i,"Selectedp0.320"])^2
         bsoutput[c("MRR","log(MRR)", "error", "log(error)", "sq.error", "log(sq.error)","MOR"),
                  k,
                  "Selectedp0.320"] <- as.numeric(Marginalest(data = databs,              # User-defined function Marginalest
                                                              int = bsSelectedoutput0320$Msel.int$coefficients["(Intercept)"],
                                                              modelcoefs = bsSelectedoutput0320$Msel$coefficients[-1],
                                                              truth = Marginal["MRR",i,"Truth"]))
      }

      # Confidence intervals marginal risk ratios full, selected, double selected
      Marginal[c("cilow", "ciup"),i,"Full"]               <- quantile(bsoutput["MRR", ,"Full"],c(0.025,0.975),names=F, na.rm=T)
      Marginal["coverage",i,"Full"]                       <- isTRUE((Marginal["cilow",i,"Full"]<=Marginal["MRR",i,"Truth"])&
                                                                           (Marginal["MRR",i,"Truth"]<=Marginal["ciup",i,"Full"]))
      Marginal[c("cilow", "ciup"),i,"Selectedp0.157"]       <- quantile(bsoutput["MRR", ,"Selectedp0.157"],c(0.025,0.975),names=F, na.rm=T)
      Marginal["coverage",i,"Selectedp0.157"]               <- isTRUE((Marginal["cilow",i,"Selectedp0.157"]<=Marginal["MRR",i,"Truth"])&
                                                                           (Marginal["MRR",i,"Truth"]<=Marginal["ciup",i,"Selectedp0.157"]))
      Marginal[c("cilow", "ciup"),i,"Selectedp0.320"]       <- quantile(bsoutput["MRR", ,"Selectedp0.320"],c(0.025,0.975),names=F, na.rm=T)
      Marginal["coverage",i,"Selectedp0.320"]               <- isTRUE((Marginal["cilow",i,"Selectedp0.320"]<=Marginal["MRR",i,"Truth"])&
                                                                      (Marginal["MRR",i,"Truth"]<=Marginal["ciup",i,"Selectedp0.320"]))
      
      Model["RMSD bAY",i,"Full"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Full"]))
      Model["RMSD bAY ratio",i,"Full"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Full"])) / Model["se(bAY)",i,"Full"]
      Model[c(paste0("incl. freq. bL1",seq(1:(S[j,"nL"]/2))), 
            paste0("incl. freq. bL2",seq(1:(S[j,"nL"]/2)))),i,"Full"] <- apply(!is.na(bsest[c(paste0("bL1",seq(1:(S[j,"nL"]/2))), 
                                                                                              paste0("bL2",seq(1:(S[j,"nL"]/2)))),,"Full"])*1, MARGIN = 1, FUN = function(x) sum(x)/B)
      
      Model["RMSD bAY",i,"Selectedp0.157"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Selectedp0.157"]))
      Model["RMSD bAY ratio",i,"Selectedp0.157"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Selectedp0.157"])) / Model["se(bAY)",i,"Full"]
      Model[c(paste0("incl. freq. bL1",seq(1:(S[j,"nL"]/2))), 
              paste0("incl. freq. bL2",seq(1:(S[j,"nL"]/2)))),i,"Selectedp0.157"] <- apply(!is.na(bsest[c(paste0("bL1",seq(1:(S[j,"nL"]/2))), 
                                                                                                paste0("bL2",seq(1:(S[j,"nL"]/2)))),,"Selectedp0.157"])*1, MARGIN = 1, FUN = function(x) sum(x)/B)
      
      Model["RMSD bAY",i,"Selectedp0.320"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Selectedp0.320"]))
      Model["RMSD bAY ratio",i,"Selectedp0.320"] <- sqrt(mean(bsoutput["sq.error(bAYboot)",,"Selectedp0.320"])) / Model["se(bAY)",i,"Full"]
      Model[c(paste0("incl. freq. bL1",seq(1:(S[j,"nL"]/2))), 
              paste0("incl. freq. bL2",seq(1:(S[j,"nL"]/2)))),i,"Selectedp0.320"] <- apply(!is.na(bsest[c(paste0("bL1",seq(1:(S[j,"nL"]/2))), 
                                                                                                paste0("bL2",seq(1:(S[j,"nL"]/2)))),,"Selectedp0.320"])*1, MARGIN = 1, FUN = function(x) sum(x)/B)
      
      
    }
  output <- list(Setup = Setup, Model = Model, Marginal = Marginal, Warning = Warning, BSWarning = BSWarning, Warning.int=Warning.int,BSWarning.int=BSWarning.int)
  saveRDS(output, file = paste0(filepath,"Simoutput/",scenario,"/",j,".rds"))
  }
  
}
