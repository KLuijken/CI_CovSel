####################################################
# Covariate selection project
#
# December 2019 
# Author: K Luijken
# File content: subfunctions performing;
#               Marginal Risk Ratio computation;
#               Backward elimination with FLIC
#               Model output storage
####################################################

#-----------------------------------------------------------------------#
# Perform numeric integration to obtain 'true MRR'
#-----------------------------------------------------------------------#

Truth <- function(S,j,bL1Y, bL1A){
  # Arguments
  #    S = scenariomatrix
  #    j = current row of scenariomatrix within function
  #    bL1Y = effect confounder L on Y, adjusted to # covariates
  #    bL1A = effect confounder L on A, adjusted to # covariates
  var.L1.star <- (S[j,"nL"]/2) + (S[j,"nL"]/2)*((S[j,"nL"]/2)-1)*S[j,"rhoL"]
  var.L2.star <- (S[j,"nL"]/2) + (S[j,"nL"]/2)*((S[j,"nL"]/2)-1)*S[j,"rhoL"]
  cov.L1L2    <- (S[j,"nL"]/2)*(S[j,"nL"]/2)*S[j,"rhoL"]
  sigma       <- matrix(cov.L1L2,ncol=2,nrow=2)
  diag(sigma) <- c(var.L1.star,var.L2.star)
  dx          <- 0.1
  values      <- seq(-4,4,by=dx)
  x           <- expand.grid(values,values)
  dL          <- dmvnorm(x,mean=c(rep(0,times=2)), sigma=sigma)
  
  pY0     		<- plogis((bL1Y*x[,1] + S[j,"bL2Y"]*x[,2]))*(dL*dx) 
  pY1     		<- plogis((bL1Y*x[,1] + S[j,"bL2Y"]*x[,2] + S[j,"bAY"]))*(dL*dx) 
  
  return(list(pY0=pY0, pY1=pY1))
}

#-----------------------------------------------------------------------#
# Compute MRR (Truth and Unadjusted)
#-----------------------------------------------------------------------#

Marginalcomp <- function(pY1,pY0,truth){
  MRR <- sum(pY1)/sum(pY0)
  log.MRR <- log(MRR)
  
  error <- MRR - truth
  log.error <- log.MRR - log(truth)
  
  sq.error <- (MRR - truth)^2
  log.sq.error <- (log.MRR - log(truth))^2
  
  MOR <- sum(pY1)/(1-sum(pY1)) / sum(pY0)/(1-sum(pY0))
  
  return(list(MRR=MRR, log.MRR=log.MRR, 
              error = error, log.error=log.error, 
              sq.error=sq.error, log.sq.error=log.sq.error,
              MOR=MOR))
}



#-----------------------------------------------------------------------#
# Perform backward elimination in Firth-corrected logistic regression models
#-----------------------------------------------------------------------#

BackwardFR <- function (object, scope, steps = 1000, slstay = 0.05, trace = TRUE, 
                        printwork = FALSE, ...) 
{
  istep <- 0
  working <- object
  if (trace) {
    cat("Step ", istep, ": starting model\n")
    if (printwork) {
      print(working)
      cat("\n\n")
    }
  }
  if (missing(scope)) 
    scope <- attr(terms(working), "term.labels")
  while (istep < steps & working$df >= 2) {
    istep <- istep + 1
    mat <- drop1(working)
    inscope <- match(scope, rownames(mat))  ###GH changed ordering of this and the next three lines
    inscope <- inscope[!is.na(inscope)]
    if (all(mat[inscope, 3] < slstay))    ### GH added inscope to subset the p-values. This was crucial!!
      break
    removal <- rownames(mat)[mat[, 3] == max(mat[inscope, 
                                                 3])]
    newform = as.formula(paste("~.-", removal))
    if (working$df == 2 | working$df == mat[mat[, 3] == max(mat[, 3]), 2]) 
      working <- update(working, formula = newform, pl = FALSE)
    else working <- update(working, formula = newform)
    if (trace) {
      cat("drop1:\n")   ### GH: added it to see what happens, can be delted 
      print(mat)           ### GH: added it to see what happens, can be delted 
      cat("\n\n")
      cat("Step ", istep, ": removed ", removal, 
          " (P=", max(mat[, 3]), ")\n")
      if (printwork) {
        print(working)
        cat("\n\n")
      }
    }
  }
  if (trace) 
    cat("\n")
  return(working)
}


#-----------------------------------------------------------------------#
# Fit full logistic regression model applying FLIC
#-----------------------------------------------------------------------#

Full <- function(S, j, datafull,i){
  # Arguments
  #    S = scenariomatrix
  #    j = current row of scenariomatrix within function
  #    data = data generated in current simulation run/Bootstrap sample
  preM        <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2))),collapse = "+")),
                                    data = datafull,
                                    firth = T,
                                    trace =F,
                                    dataout = T,
                                    pl = F))
  warning     <- ifelse(is.null(preM$warning),"NULL",{paste0(preM$warning,"scen=",j,"nsim=",i)}) # Store warning, if any, and location in simulation repetitions
  M           <- preM$value
  M.pred      <- as.matrix(datafull[,names(datafull)!="Y"]) %*% M$coefficients[-1]
  pre.M.int   <- tryCatch.W.E(glm(datafull$Y ~ offset(M.pred), family = binomial))   # Re-estimate the intercept
  warning.int <- ifelse(is.null(pre.M.int$warning),"NULL",{paste0(pre.M.int$warning,"scen=",j,"nsim=",i)}) # Store warning, if any
  M.int       <- pre.M.int$value
  
  return(list(warning=warning, M=M, M.int=M.int, warning.int=warning.int))
}


#-----------------------------------------------------------------------#
# Fit selected logistic regression model applying FLIC
#-----------------------------------------------------------------------#

Selected <- function(Mod,pcutoff,S, j, datafull,i){
  # Arguments
  #    Mod = Full model
  #    pcutoff = criterion for backward elimination
  #    S = scenariomatrix
  #    j = current row of scenariomatrix within function
  #    datafull = data generated in current simulation run
  preMsel       <- tryCatch.W.E(BackwardFR(Mod, 
                                       sls = pcutoff, 
                                       trace=F,
                                       scope = c(paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2))),
                                       printwork = F,
                                       datafull = datafull,
                                       pl = F))
  warning       <- ifelse(is.null(preMsel$warning),"NULL",{paste0(preMsel$warning,"scen=",j,"nsim=",i)}) # Store warning
  Msel          <- preMsel$value
  Msel.pred     <- as.matrix(datafull[,names(coef(Msel))[-1]]) %*% coef(Msel)[-1]
  pre.Msel.int  <- tryCatch.W.E(glm(datafull$Y ~ offset(Msel.pred), family = binomial)) # Re-estimate intercept
  warning.int   <- ifelse(is.null(pre.Msel.int$warning),"NULL",{paste0(pre.Msel.int$warning,"scen=",j,"nsim=",i)}) # Store warning
  Msel.int      <- pre.Msel.int$value
  
  
  return(list(warning=warning, Msel=Msel, Msel.int=Msel.int, warning.int=warning.int))
}




#-----------------------------------------------------------------------#
# Obtain model coefficients
#-----------------------------------------------------------------------#

Modeloutput <- function(intmodel,model){
  Int     <- intmodel$coefficients[1]
  se.Int  <- coef(summary(intmodel))[1,2]
  coef_est <-  coef(model)[c("(Intercept)","A",
                             paste0("L1",seq(1:(max(nL)/2))),
                             paste0("L2",seq(1:(max(nL)/2))))]
  coef_est[is.na(coef_est)] <- NA
  coef_est                 <- coef_est[-1] 
  
  coef_se <-  rep(NA, times= max(nL)+2)
  names(coef_se) <- c("(Intercept)","A",
                      paste0("L1",seq(1:(max(nL)/2))),
                      paste0("L2",seq(1:(max(nL)/2))))
  coef_se[c(names(model$coefficients))] <- sqrt(diag(vcov(model, complete=T)))
  coef_se                 <- coef_se[-1]
  
  alloutput <- c(Int,coef_est,se.Int,coef_se)
  
  return(alloutput)
}

#-----------------------------------------------------------------------#
# Estimate MRR (Full model and Selected models)
#-----------------------------------------------------------------------#

Marginalest <- function(data,int,modelcoefs,truth){
  newdat <- data.frame(data[,c(names(modelcoefs))])
  colnames(newdat) <- c(names(modelcoefs))
  ifelse("A" %in% c(names(modelcoefs)),{             # In order to catch boostrap samples that have too strong multicollinearity
                                                     # Note to self: add counter to keep track in how many samples this occurs
  newdat[,"A"] <- rep(0,times=nrow(newdat))
  PredA0 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
  newdat[,"A"] <- rep(1,times=nrow(newdat))
  PredA1 <- plogis(rep(1,times=nrow(newdat)) * int + as.matrix(newdat) %*% matrix(modelcoefs))
  
  MRR <- mean(PredA1)/mean(PredA0)
  log.MRR <- log(MRR)
  
  error <- MRR - truth
  log.error <- log.MRR - log(truth)
  
  sq.error <- (MRR - truth)^2
  log.sq.error <- (log.MRR - log(truth))^2
  
  MOR <- sum(PredA1)/(1-sum(PredA1)) / sum(PredA0)/(1-sum(PredA0))},
  
  MRR <- log.MRR <- error <- log.error <- sq.error <- log.sq.error <- MOR <- NA)
  
  return(list(MRR=MRR, log.MRR=log.MRR, 
              error = error, log.error=log.error, 
              sq.error=sq.error, log.sq.error=log.sq.error,
              MOR=MOR))
}


#-----------------------------------------------------------------------#
# Bootstrap
#-----------------------------------------------------------------------#
# NB: the reason to create separate bootstrap functions is that the dataframes
# are not passed on correctly if these aren't separate arguments --> note to self: think about way to improve this

Fullbs <- function(S, j, databs, i, k){
  # Arguments
  #    S = scenariomatrix
  #    j = current row of scenariomatrix within function
  #    databs = data generated in current Bootstrap sample
  preM        <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2))),collapse = "+")),
                                    data = databs,
                                    firth = T,
                                    trace =F,
                                    dataout = T,
                                    pl = F))
  warning     <- ifelse(is.null(preM$warning),"NULL",{paste0(preM$warning,"scen=",j,"nsim=",i,"Boot=",k)}) # Store warning, if any
  M           <- preM$value
  M.pred      <- as.matrix(databs[,names(databs)!="Y"]) %*% M$coefficients[-1]
  pre.M.int   <- tryCatch.W.E(glm(databs$Y ~ offset(M.pred), family = binomial))   # Re-estimate the intercept
  warning.int <- ifelse(is.null(pre.M.int$warning),"NULL",{paste0(pre.M.int$warning,"scen=",j,"nsim=",i,"Boot=",k)}) # Store warning
  M.int       <- pre.M.int$value
  
  return(list(warning=warning, M=M, M.int=M.int, warning.int=warning.int))
}

Selectedbs <- function(Mod,pcutoff,S, j, databs, i, k){
  # Arguments
  #    Mod = Full model
  #    pcutoff = criterion for backward elimination
  #    S = scenariomatrix
  #    j = current row of scenariomatrix within function
  #    databs = data generated in current Bootstrap sample
  preMsel       <- tryCatch.W.E(BackwardFR(Mod, 
                                       sls = pcutoff, 
                                       trace=F,
                                       scope = c(paste0("L1",1:(S[j,"nL"]/2)),paste0("L2",1:(S[j,"nL"]/2))),
                                       printwork = F,
                                       databs = databs,
                                       pl = F))
  warning       <- ifelse(is.null(preMsel$warning),"NULL",{paste0(preMsel$warning,"scen=",j,"nsim=",i,"Boot=",k)}) # Store warning
  Msel          <- preMsel$value
  Msel.pred     <- as.matrix(databs[,names(coef(Msel))[-1]]) %*% coef(Msel)[-1]
  pre.Msel.int  <- tryCatch.W.E(glm(databs$Y ~ offset(Msel.pred), family = binomial)) # Re-estimate intercept
  warning.int   <- ifelse(is.null(pre.Msel.int$warning),"NULL",{paste0(pre.Msel.int$warning,"scen=",j,"nsim=",i,"Boot=",k)}) # Store warning
  Msel.int      <- pre.Msel.int$value
  
  return(list(warning=warning, Msel=Msel, Msel.int=Msel.int, warning.int=warning.int))
}
