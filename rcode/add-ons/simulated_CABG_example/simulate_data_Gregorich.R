#------------------------------------------------------------------------------#
# Causal inference covariate selection
# Mariella Gregorich 
# December 2017
#
# Generate data simular to CABG example
#------------------------------------------------------------------------------#

# Load libraries ----
#------------------------------------------------------------------------------#
library(bindata)
library(mvtnorm)

# Generate data ----
#------------------------------------------------------------------------------#

### Determine Correlations
# x_tmp <- model.matrix(CT~Age+Gender+Smoker+Diabetes.Control+CreaCl+Dialysis+Hypertension+
#                         Peripheral.Vascular.Disease+Cerebrovascular.Disease+Cerebrovascular.Accident+
#                         Myocardial.Infarction+Congestive.Heart.Failure+Angina.Type+Afib.flutter+
#                         Number.of.Diseased.Coronary.Vessels+Left.Main.Disease+Ejection.Fraction+
#                         Status+Dyslipidemia+Lipid.Lowering+Previous.Coronary.Artery.Bypass+Previous.Valve+
#                         Year.CABG, data=dataCABG)
# X <- cor(x_tmp)
# corr_matrix <- X[2:nrow(X), 2:ncol(X)]


generate_data <- function(N, betaTr.zero = F, avsu = F){
  
  ### Categorize the variables
  var_binary <- c("Gender", "Smoker", "Dialysis", "Hypertension", "Peripheral.Vascular.Disease", "Cerebrovascular.Disease",
                  "Cerebrovascular.Accident", "Myocardial.Infarction", "Congestive.Heart.Failure","Afib.flutter",
                  "Left.Main.Disease", "Previous.Coronary.Artery.Bypass",
                  "Previous.Valve")
  var_cont <- c("Age", "CreaCl","Ejection.Fraction", "Year.CABG" )
  var_cat <- c("Diabetes.Control", "Angina.Type", "Number.of.Diseased.Coronary.Vessels", "Status")
  # prev_binary <- sapply(var_binary, function(x)  sum(dataCABG[,x]==1)/nrow(dataCABG))
  prev_binary <- c("Gender"=0.19329214, "Smoker"=0.45542807, "Dialysis"=0.02118270, "Hypertension"=0.90158870, 
                   "Peripheral.Vascular.Disease"=0.20432480, "Cerebrovascular.Disease"=0.21006178, 
                   "Cerebrovascular.Accident"=0.07899382, "Myocardial.Infarction"=0.51368049, 
                   "Congestive.Heart.Failure"=0.11782877, "Afib.flutter"=0.09399823, "Left.Main.Disease"=0.37819947, 
                   "Previous.Coronary.Artery.Bypass"=0.01235658, "Previous.Valve"=0.00353045)

  ## Set 1: CHF (0.1) - -0.5 - Ejection Fraction - Myocardial Infarction
  m.sigma <-  matrix(c(1,0,-0.3, 0,1,-0.5, -0.3,-0.5,1),3,3)
  tmp <- rmvnorm(N, mean = rep(0,3), sigma = m.sigma, method = "chol")
  myoinf <- (tmp[,1] > qnorm(1-prev_binary["Myocardial.Infarction"]))*1
  CHF <- (tmp[,2] > qnorm(1-prev_binary["Congestive.Heart.Failure"]))*1
  eject.frac <- as.integer(12.81585*tmp[,3] + 51.94156)

  ## Set 2: Smoker - Age - CreaCl - Gender
  #m.sigma <- round(cor(model.matrix(~Smoker+Age+CreaCl+Gender, data=dataCABG))[2:5,2:5],2)
  m.sigma <- matrix(c(1.00, -0.32, 0.23, -0.14, 
                      -0.32, 1.00, -0.55, 0.11, 
                      0.23, -0.55, 1.00, -0.25,
                      -0.14, 0.11, -0.25, 1.00),ncol = 4, nrow=4)
  tmp <- rmvnorm(N ,mean = rep(0,4), sigma = m.sigma, method = "chol")
  smoker <- (tmp[,1] > qnorm(1-prev_binary["Smoker"]))*1
  age <- pmin(90, as.integer(10.05103*tmp[,2] + 66.75115))
  m <- 82.09253          # mean(dataCABG$CreaCl)
  s <- 32.90679
  loc <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  crea <- as.integer(exp(tmp[,3]*shape + loc))
  IQR <- quantile(crea, c(0.25,0.75))
  IQRdist <- IQR[2]-IQR[1]
  crea[crea > IQR[2]+5*IQRdist] <- IQR[2]+5*IQRdist
  gender <- (tmp[,4] > qnorm(0.8))*1
  
  subgroup <- data.frame("CreaCl"=crea, "Dialysis"=rnorm(N))
  subgroup[subgroup$CreaCl <= 60,]$Dialysis <- subgroup[subgroup$CreaCl <= 60, ]$Dialysis > qnorm(0.92)
  subgroup[subgroup$CreaCl > 60,]$Dialysis <- subgroup[subgroup$CreaCl > 60, ]$Dialysis > qnorm(0.997)
  dialysis <- subgroup$Dialysis
  
  ## Set 3: Cerebro. Accident (0.07) - 0.6 - Cerebro. Disease (0.2)
  m.sigma <- matrix(c(1,0.6,0.6,1),2,2)
  cereA.cereD <- rmvbin(N, margprob = c(0.1,0.2), bincorr = m.sigma) 
  colnames(cereA.cereD) <- c("Cerebrovascular.Accident", "Cerebrovascular.Disease")
  
  ## Set 4: Status - (Emergent - 0.3 - Unstable) - Angina.Type (subgrouping)
  tmp <- rnorm(N)
  df <- data.frame("Status"=rep(NA,N), "Angina.Type"=rnorm(N))
  prev_status <- c("Elective"= 0.64298323, "Emergent"=0.07149162, "Urgent"=0.28552515) 
  
  df[,"Status"] <- (tmp > qnorm(0.64)) + (tmp > qnorm(0.71)) 
  
  # table_subgroup <- table(dataCABG$Status, dataCABG$Angina.Type)
  # prev_angit <- table_subgroup / apply(table_subgroup,1, sum)
  df[df$Status==0,]$Angina.Type <- (df[df$Status==0,]$Angina.Typ > qnorm(0.24)) + (df[df$Status==0,]$Angina.Typ > qnorm(0.88))
  df[df$Status==1,]$Angina.Type <- (df[df$Status==1,]$Angina.Typ > qnorm(0.12)) + (df[df$Status==1,]$Angina.Typ > qnorm(0.20))
  df[df$Status==2,]$Angina.Type <- (df[df$Status==2,]$Angina.Typ > qnorm(0.16)) + (df[df$Status==2,]$Angina.Typ > qnorm(0.59))
  
  
  ## Set 5: Hypertension - 0.25 - Dyslipidemia (0.8) - 0.3 - Lipid.Lowering (0.7) - 0.1
  m.sigma <- matrix(c(1,0.25,0.1, 0.25, 1,0.3, 0.1,0.3,1),3,3)
  hyp.dys.liplo <- rmvbin(N, margprob = c(0.9,0.8,0.7), bincorr = m.sigma) 
  colnames(hyp.dys.liplo) <- c("Hypertension", "Dyslipidemia", "Lipid.Lowering")
  
  ## Set 6: Previous Valve (0.0035) - 0.3 - Previous CABG (0.01)
  tmp <- data.frame("Previous.Valve"=rnorm(N), "Previous.Coronary.Artery.Bypass"=rnorm(N))
  tmp$Previous.Valve  <-  (tmp[,1] > qnorm(0.996))*1
  tmp[tmp$Previous.Valve==1,]$Previous.Coronary.Artery.Bypass <- (tmp[tmp$Previous.Valve==1,]$Previous.Coronary.Artery.Bypass > qnorm(0.5))*1
  tmp[tmp$Previous.Valve==0,]$Previous.Coronary.Artery.Bypass <- (tmp[tmp$Previous.Valve==0,]$Previous.Coronary.Artery.Bypass > qnorm(0.99))*1
  prev.valve.prev.cabg <- tmp
  
  
  ### Generate binary variables
  prev_binary_uncor <-prev_binary[c("Peripheral.Vascular.Disease", "Afib.flutter", "Left.Main.Disease")]
  X_binvar <- sapply(prev_binary_uncor, function(x) rbinom(N, 1, x))
  
  ### Generate Year of CABG
  year <- as.integer(7*pnorm(rnorm(N))+2009)
  
  ### Generate categorical variables: Diabetes, NodCV
  # prev_diab <- table(dataCABG$Diabetes.Control)/nrow(dataCABG)
  tmp <- rnorm(N)
  diabetes <- (tmp > qnorm(0.57)) + (tmp > qnorm(0.82)) + (tmp > qnorm(0.95)) + (tmp > qnorm(0.986)) 
  
  # prev_NodCV <- table(dataCABG$Number.of.Diseased.Coronary.Vessels)/nrow(dataCABG)
  tmp <- rnorm(N)
  NodCV <- (tmp > qnorm(0.03)) + (tmp > qnorm(0.15)) 
  
  
  ### Assemble simulated covariates
  sim_data <- data.frame("Age"=age, "Gender"=gender, "Smoker"=smoker, "Diabetes.Control"=diabetes, 
                         "CreaCl"=crea, "Dialysis"=dialysis ,"Hypertension"=hyp.dys.liplo[,"Hypertension"],
                         "Peripheral.Vascular.Disease"=X_binvar[,c("Peripheral.Vascular.Disease")],
                         cereA.cereD, "Myocardial.Infarction"=myoinf, 
                         "Congestive.Heart.Failure"=CHF, "Angina.Type"=df$Angina.Type,"Afib.flutter"=X_binvar[,c("Afib.flutter")],
                         "Number.of.Diseased.Coronary.Vessels"=NodCV, "Left.Main.Disease"=X_binvar[,c("Left.Main.Disease")], 
                         "Ejection.Fraction"=eject.frac, "Status"=df$Status, hyp.dys.liplo[,c("Dyslipidemia","Lipid.Lowering")],
                         prev.valve.prev.cabg, "Year.CABG"=year)
  sim_data$Status <- factor(sim_data$Status)
  levels(sim_data$Status) <- list(Elective="0", Emergent="1", Urgent="2")
  sim_data$Angina.Type <- factor(sim_data$Angina.Type)
  levels(sim_data$Angina.Type) <- list("0"="0",Stable="1", Unstable="2")
  sim_data$Number.of.Diseased.Coronary.Vessels <- factor(sim_data$Number.of.Diseased.Coronary.Vessels)
  levels(sim_data$Number.of.Diseased.Coronary.Vessels) <- list(One="0", Two="1", Three="2")
  sim_data$Diabetes.Control <- factor(sim_data$Diabetes.Control)
  levels(sim_data$Diabetes.Control) <- list(NoDiabetes="0", Oral="1", Insulin="2", Diet="3", None="4")
  
  
  ### Generate error terms u (Treatment) & v (Outcome)
  m.sigma <- matrix(c(1,0.5,0.5,1), 2, 2)
  error_terms <- data.frame(rmvnorm(N, mean = rep(0,2), sigma = m.sigma, method = "chol"))
  colnames(error_terms) <- c("u", "v")
  
  #### Simulate the treatment CT 

  generate_treatment <- function(df,avsu) {
    if(avsu == T){
      coeff <- c(c(0.0002, 1.8e-05, 5.7607462e-05, 0.0001552214, 7.9430717e-05)*2, 0.39362988, 0.011369225, 0.026528193, 0.015158967, 
               c(2.6477438, 0.18280118)*2, 0.098673493 * 2, c(0.16173817, 0.057619453, 0.015288401, 0.18238363) * 2, c(0.0028877989, 2.6629914e-07,
               7.9335682e-07, 1.7208013e-06, 6.6114529e-07) * 2, 0.58661935 * 2, 0.76895577 * 2, 0.12433038 * 2, 0.30601205 * 2, 1.2538594 * 2, 
               0.6998369 * 2, 0.25876773 * 2)
      intercept <- -793
    } else{
      coeff <- c(0.0002, 1.8e-05, 5.7607462e-05, 0.0001552214, 7.9430717e-05, 0.39362988, 0.011369225, 0.026528193, 0.015158967, 
                 2.6477438, 0.18280118, 0.098673493, 0.16173817, 0.057619453, 0.015288401, 0.18238363, 0.0028877989, 2.6629914e-07,
                 7.9335682e-07, 1.7208013e-06, 6.6114529e-07, 0.58661935, 0.76895577, 0.12433038, 0.30601205, 1.2538594, 0.6998369,
                 0.25876773)
      intercept <- -793      
    }

    formula_treat <- function(coeff, intercept) {
      tmp <- intercept + coeff[1]*df$Age+coeff[2]*pmax(df$Age-49,0)^3+coeff[3]*pmax(df$Age-63.635892,0)^3-coeff[4]*pmax(df$Age-71.318767,0)^3+
        coeff[5]*pmax(df$Age-82,0)^3+coeff[6]*df$Year.CABG+coeff[7]*pmax(df$Year.CABG-2009,0)^3-coeff[8]*pmax(df$Year.CABG-2013,0)^3+
        coeff[9]*pmax(df$Year.CABG-2016,0)^3-coeff[10]*(df$Status=="Emergent")-coeff[11]*(df$Status=="Urgent")+coeff[12]*(df$Smoker=="1")-
        coeff[13]*(df$Diabetes.Control=="Insulin")-coeff[14]*(df$Diabetes.Control=="NoDiabetes")-coeff[15]*(df$Diabetes.Control=="None")-
        coeff[16]*(df$Diabetes.Control=="Oral")+coeff[17]*df$CreaCl-coeff[18]*pmax(df$CreaCl-36.2025,0)^3-
        coeff[19]*pmax(df$CreaCl-67.709583,0)^3+coeff[20]*pmax(df$CreaCl-90.909091,0)^3-coeff[21]*pmax(df$CreaCl-140.78286,0)^3+
        coeff[22]*(df$Dialysis=="1")+coeff[23]*(df$Peripheral.Vascular.Disease=="1")+coeff[24]*(df$Cerebrovascular.Disease=="1")+
        coeff[25]*(df$Cerebrovascular.Accident=="1")+coeff[26]*(df$Previous.Coronary.Artery.Bypass=="1")+coeff[27]*(df$Previous.Valve=="1")-
        coeff[28]*df$Hypertension
      tmp
      return(tmp)
      }
    tmp <- formula_treat(coeff, intercept)
    prob = exp(tmp)/(1+exp(tmp))
    prob
    return(prob)}
  
  prob_CT <- generate_treatment(sim_data, avsu = T)
  CT_tmp = rbinom(N,1,prob_CT)
  
  sum(CT_tmp==1) #1132
  sum(CT_tmp==0) #1368
  
  CT <- ifelse(CT_tmp==1, CT <- 0.5, CT <- -0.5)
  sim_data <- cbind(CT,sim_data)
  
  
  ### Simulate the outcome stroke 
  generate_outcome <- function(df, betaTr.zero, avsu){
    if(avsu == TRUE) {
      if(betaTr.zero==TRUE) {
        coeff <- c(-3.131/2, "beta.Treat" = 0, 0.004, 0.601, 0.153, 0.438, 0.455, 0.137, 0.024,
                 0.672, -0.031, -1.130, -0.015,  0.682,  0.215, -0.050, 1.024, 1.754, 1.754, -0.442)*2
       }else{coeff <- c(-3.131/2, "beta.Treat" = -0.1, 0.004, 0.601, 0.153, 0.438, 0.455, 0.137, 0.024,
                 0.672, -0.031, -1.130, -0.015,  0.682,  0.215, -0.050, 1.024, 1.754, 1.754, -0.442) * 2}
    }else{ #avsu=FALSE
      if(betaTr.zero==TRUE) {
        coeff <- c(-3.131, "beta.Treat" = 0, 0.004, 0.601, 0.153, 0.438, 0.455, 0.137, 0.024,
                   0.672, -0.031, -1.130, -0.015,  0.682,  0.215, -0.050, 1.024, 1.754, 1.754, -0.442)
      }else{coeff <- c(-3.131, "beta.Treat" = -0.890, 0.004, 0.601, 0.153, 0.438, 0.455, 0.137, 0.024,
                      0.672, -0.031, -1.130, -0.015,  0.682,  0.215, -0.050, 1.024, 1.754, 1.754, -0.442)}
    }
    X <- model.matrix(~ CT + Age + Status + Gender + Smoker + Diabetes.Control + Afib.flutter +
                          CreaCl + Dialysis + Peripheral.Vascular.Disease + Cerebrovascular.Disease + 
                          Cerebrovascular.Accident + Previous.Coronary.Artery.Bypass + Previous.Valve + Hypertension, data=df)
    logprob <- apply(coeff * t(X), 2, sum)
    prob <- exp(logprob)/(1+exp(logprob))
    return(prob)}
  
  prob_stroke <- generate_outcome(sim_data, betaTr.zero, avsu)
  stroke = rbinom(N,1,prob_stroke)
  #stroke <- factor(stroke)

  sim_data <- cbind("Postoperative.stroke"=stroke, sim_data)
  sim_data$CT <- CT_tmp
  return(sim_data)
}


# m.cor <- cor(model.matrix(~Age+Gender+Smoker+Diabetes.Control+CreaCl+Dialysis+Hypertension+
#                             Peripheral.Vascular.Disease+Cerebrovascular.Disease+Cerebrovascular.Accident+
#                             Myocardial.Infarction+Congestive.Heart.Failure+Angina.Type+Afib.flutter+
#                             Number.of.Diseased.Coronary.Vessels+Left.Main.Disease+Ejection.Fraction+
#                             Status+Dyslipidemia+Lipid.Lowering+Previous.Coronary.Artery.Bypass+Previous.Valve+
#                             Year.CABG, data=sim_data))
# m.cor.true <- cor(model.matrix(~Age+Gender+Smoker+Diabetes.Control+CreaCl+Dialysis+Hypertension+
#                                  Peripheral.Vascular.Disease+Cerebrovascular.Disease+Cerebrovascular.Accident+
#                                  Myocardial.Infarction+Congestive.Heart.Failure+Angina.Type+Afib.flutter+
#                                  Number.of.Diseased.Coronary.Vessels+Left.Main.Disease+Ejection.Fraction+
#                                  Status+Dyslipidemia+Lipid.Lowering+Previous.Coronary.Artery.Bypass+Previous.Valve+
#                                  Year.CABG, data=dataCABG))
