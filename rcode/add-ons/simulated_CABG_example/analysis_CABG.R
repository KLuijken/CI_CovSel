#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Perform CABG example analysis
#------------------------------------------------------------------------------#

# This script describes the analysis on simulated data of the example CABG study
# presented in manuscript section 3. The following analyses are presented:
## conditional odds ratio of a full model, estimated using FLIC
## conditional odds ratio of a full model, estimated using FLIC and backward 
#     elimination
## marginal odds ratio of a full model, estimated using predicted potential 
#     outcomes of a FLIC model
## marginal odds ratio of a full model, estimated using predicted potential 
#     outcomes of a FLIC model and backward elimination


# Load librairies + source code ----
#------------------------------------------------------------------------------#
library(logistf)

source(file = "./rcode/add-ons/simulated_CABG_example/simulate_data_Gregorich.R")

# Generate simulated data ----
#------------------------------------------------------------------------------#

set.seed(20200618)
CABG_data<- generate_data(N = 2266,
              betaTr.zero = FALSE,
              avsu = FALSE)

# Conditional OR, FLIC, full model ----
#------------------------------------------------------------------------------#

# Estimate model
Firth_full   <- logistf(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                         Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                         Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                         Cerebrovascular.Disease + Myocardial.Infarction + 
                         Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                         Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                         Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                         Previous.Valve + Previous.Coronary.Artery.Bypass +
                         Year.CABG,
                       data = CABG_data,
                       firth = TRUE, # set firth = FALSE for maximum likelihood estimation
                       pl = TRUE)    # compute profile likelihood confidence intervals

# Obtain coefficient + CI
cOR_full_Firth_Est <- round( exp( summary( Firth_full)$coefficients["CT"]), digits=2)
cOR_full_Firth_Low <- round( exp( summary( Firth_full)$ci.lower["CT"]), digits=2)
cOR_full_Firth_Up  <- round( exp( summary( Firth_full)$ci.upper["CT"]), digits=2)
paste0(cOR_full_Firth_Est,
       "(95% CI, ", cOR_full_Firth_Low,
       "; ", cOR_full_Firth_Up,")")

# Conditional OR, FLIC, backward elimination ----
#------------------------------------------------------------------------------#

# Estimate model
Firth_selected      <- backward(Firth_full,
                                scope = c("Age", "Gender", "Smoker", 
                                          "Diabetes.Control", "CreaCl", "Dialysis", 
                                          "Hypertension", "Peripheral.Vascular.Disease",
                                          "Cerebrovascular.Accident", 
                                          "Cerebrovascular.Disease", "Myocardial.Infarction", 
                                          "Congestive.Heart.Failure", "Angina.Type", 
                                          "Afib.flutter", "Number.of.Diseased.Coronary.Vessels",
                                          "Left.Main.Disease", "Ejection.Fraction",
                                          "Status", "Dyslipidemia", "Lipid.Lowering",
                                          "Previous.Valve", "Previous.Coronary.Artery.Bypass",
                                          "Year.CABG"),
                                slstay = 0.157,
                                trace = FALSE)

# Obtain coefficient + CI
cOR_selected_Firth_Est <- round( exp( summary( Firth_selected)$coefficients["CT"]), digits=2)
cOR_selected_Firth_Low <- round( exp( summary( Firth_selected)$ci.lower["CT"]), digits=2)
cOR_selected_Firth_Up  <- round( exp( summary( Firth_selected)$ci.upper["CT"]), digits=2)
paste0(cOR_selected_Firth_Est,
       "(95% CI, ", cOR_selected_Firth_Low,
       "; ", cOR_selected_Firth_Up,")")


# Marginal OR, FLIC, full model ----
#------------------------------------------------------------------------------#

# Use Firth_full model previous step.
# Re-estimate the intercept
pred_intercept       <- Firth_full$linear.predictors
intercept_Firth_full <- glm(CABG_data$Postoperative.stroke ~ offset(pred_intercept), 
                            family = binomial)$coefficients[1]

# Obtain predicted potential outcomes
mod_mat <- model.matrix(as.formula(
  Firth_full$call[["formula"]]),
  data = CABG_data)

mod_mat[,"CT"] <- 0
PredA0 <- plogis(as.vector(mod_mat %*% Firth_full$coefficients + intercept_Firth_full))
mod_mat[,"CT"] <- 1
PredA1 <- plogis(as.vector(mod_mat %*% Firth_full$coefficients + intercept_Firth_full))

# Estimate Marginal OR
MOR_full_Firth_Est <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))

# Bootstrap confidence interval:
b_rep <- 100
seeds <- sample(1:100000, size = b_rep, replace = FALSE)
CI_marginal_OR <- matrix(NA, nrow = b_rep, ncol = 1)

for(i in 1:b_rep){
  # Sample bootstrap data using prespecified seeds
  set.seed(seeds[i])
  sampled_rows <- sample(1:nrow(CABG_data), size = nrow(CABG_data), replace = TRUE)
  bs_sample <- CABG_data[sampled_rows,]
  
  # Estimate Firth model in bootstrap sample
  Firth_full_bs   <- logistf(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                            Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                            Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                            Cerebrovascular.Disease + Myocardial.Infarction + 
                            Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                            Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                            Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                            Previous.Valve + Previous.Coronary.Artery.Bypass +
                            Year.CABG,
                          data = bs_sample,
                          firth = TRUE, # set firth = FALSE for maximum likelihood estimation
                          pl = FALSE)   # do not compute profile likelihood confidence intervals
  
  
  # Re-estimate intercept in bootrap sample
  pred_intercept_bs       <- Firth_full_bs$linear.predictors
  intercept_Firth_full_bs <- glm(bs_sample$Postoperative.stroke ~ offset(pred_intercept_bs), 
                              family = binomial)$coefficients[1]
  
  mod_mat <- model.matrix(as.formula(
    Firth_full_bs$call[["formula"]]),
    data = bs_sample)
  
  # Estimate MOR in bootrap sample
  mod_mat[,"CT"] <- 0
  PredA0 <- plogis(as.vector(mod_mat %*% Firth_full_bs$coefficients + intercept_Firth_full_bs))
  mod_mat[,"CT"] <- 1
  PredA1 <- plogis(as.vector(mod_mat %*% Firth_full_bs$coefficients + intercept_Firth_full_bs))
  
  CI_marginal_OR[i,] <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))
  
}

# Obtain MOR + CI
MOR_full_Firth_Low <- round(quantile(CI_marginal_OR, 0.025), digits = 2)
MOR_full_Firth_Up  <- round(quantile(CI_marginal_OR, 0.975), digits = 2)
paste0(round(MOR_full_Firth_Est,digits=2),
       "(95% CI, ", MOR_full_Firth_Low,
       "; ", MOR_full_Firth_Up,")")


# Marginal OR, FLIC, backward elimination ----
#------------------------------------------------------------------------------#

# Use Firth_selected model previous step.
# Re-estimate the intercept
pred_intercept       <- Firth_selected$linear.predictors
intercept_Firth_selected <- glm(CABG_data$Postoperative.stroke ~ offset(pred_intercept), 
                            family = binomial)$coefficients[1]

# Obtain predicted potential outcomes
mod_mat <- model.matrix(as.formula(
  Firth_selected$call[["formula"]]),
  data = CABG_data)

mod_mat[,"CT"] <- 0
PredA0 <- plogis(as.vector(mod_mat %*% Firth_selected$coefficients + intercept_Firth_selected))
mod_mat[,"CT"] <- 1
PredA1 <- plogis(as.vector(mod_mat %*% Firth_selected$coefficients + intercept_Firth_selected))

# Estimate Marginal OR
MOR_selected_Firth_Est <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))

# Bootstrap confidence interval:
b_rep <- 100
seeds <- sample(1:100000, size = b_rep, replace = FALSE)
CI_marginal_OR <- matrix(NA, nrow = b_rep, ncol = 1)

for(i in 1:b_rep){
  # Sample bootstrap data using prespecified seeds
  set.seed(seeds[i])
  sampled_rows <- sample(1:nrow(CABG_data), size = nrow(CABG_data), replace = TRUE)
  bs_sample <- CABG_data[sampled_rows,]
  
  # Estimate Firth model in bootstrap sample
  Firth_full_bs   <- logistf(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                            Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                            Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                            Cerebrovascular.Disease + Myocardial.Infarction + 
                            Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                            Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                            Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                            Previous.Valve + Previous.Coronary.Artery.Bypass +
                            Year.CABG,
                          data = bs_sample,
                          firth = TRUE, # set firth = FALSE for maximum likelihood estimation
                          pl = FALSE)   # do not compute profile likelihood confidence intervals
  
  # Perform backward elimination in bootstrap sample
  Firth_selected_bs      <- backwardf(Firth_full_bs,
                                   scope = c("Age", "Gender", "Smoker", 
                                             "Diabetes.Control", "CreaCl", "Dialysis", 
                                             "Hypertension", "Peripheral.Vascular.Disease",
                                             "Cerebrovascular.Accident", 
                                             "Cerebrovascular.Disease", "Myocardial.Infarction", 
                                             "Congestive.Heart.Failure", "Angina.Type", 
                                             "Afib.flutter", "Number.of.Diseased.Coronary.Vessels",
                                             "Left.Main.Disease", "Ejection.Fraction",
                                             "Status", "Dyslipidemia", "Lipid.Lowering",
                                             "Previous.Valve", "Previous.Coronary.Artery.Bypass",
                                             "Year.CABG"),
                                   slstay = 0.157,
                                   trace = FALSE,
                                   data = bs_sample)
  
  # Re-estimate intercept in bootrap sample
  pred_intercept_bs       <- Firth_selected_bs$linear.predictors
  intercept_Firth_selected_bs <- glm(bs_sample$Postoperative.stroke ~ offset(pred_intercept_bs), 
                                 family = binomial)$coefficients[1]
  
  mod_mat <- model.matrix(as.formula(
    Firth_selected_bs$call[["formula"]]),
    data = bs_sample)
  
  # Estimate MOR in bootrap sample
  mod_mat[,"CT"] <- 0
  PredA0 <- plogis(as.vector(mod_mat %*% Firth_selected_bs$coefficients + intercept_Firth_selected_bs))
  mod_mat[,"CT"] <- 1
  PredA1 <- plogis(as.vector(mod_mat %*% Firth_selected_bs$coefficients + intercept_Firth_selected_bs))
  
  CI_marginal_OR[i,] <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))
  
}

# Obtain MOR + CI
MOR_selected_Firth_Low <- round(quantile(CI_marginal_OR, 0.025), digits = 2)
MOR_selected_Firth_Up  <- round(quantile(CI_marginal_OR, 0.975), digits = 2)
paste0(round(MOR_selected_Firth_Est,digits=2),
       "(95% CI, ", MOR_selected_Firth_Low,
       "; ", MOR_selected_Firth_Up,")")


# Investigate individual patient risks ----
#------------------------------------------------------------------------------#

# Generate artificial high and low risk patients
patient_high_exp <- data.frame(CT = 1, Age = 70, Gender = 1, Smoker = 1,
                               Diabetes.Control = "Insulin", CreaCl = 45,
                               Dialysis = 1, Hypertension = 1,
                               Peripheral.Vascular.Disease = 1,
                               Cerebrovascular.Accident = 1,
                               Cerebrovascular.Disease = 1,
                               Myocardial.Infarction =1,
                               Congestive.Heart.Failure = 1,
                               Angina.Type = "Unstable",
                               Afib.flutter = 1,
                               Number.of.Disease.Cornary.Vessels = "Three",
                               Left.Main.Disease = 1,
                               Ejection.Fraction = 21,
                               Status = "Urgent",
                               Dyslipidemia = 1,
                               Lipid.Lowering = 1,
                               Previous.Valve = 1,
                               Previous.Coronary.Artery.Bypass = 1,
                               Year.CABG = 2012)
patient_high_unexp <- data.frame(CT = 0, Age = 70, Gender = 1, Smoker = 1,
                                 Diabetes.Control = "Insulin", CreaCl = 45,
                                 Dialysis = 1, Hypertension = 1,
                                 Peripheral.Vascular.Disease = 1,
                                 Cerebrovascular.Accident = 1,
                                 Cerebrovascular.Disease = 1,
                                 Myocardial.Infarction =1,
                                 Congestive.Heart.Failure = 1,
                                 Angina.Type = "Unstable",
                                 Afib.flutter = 1,
                                 Number.of.Disease.Cornary.Vessels = "Three",
                                 Left.Main.Disease = 1,
                                 Ejection.Fraction = 21,
                                 Status = "Urgent",
                                 Dyslipidemia = 1,
                                 Lipid.Lowering = 1,
                                 Previous.Valve = 1,
                                 Previous.Coronary.Artery.Bypass = 1,
                                 Year.CABG = 2012)

patient_low_exp <- data.frame(CT = 1, Age = 45, Gender = 0, Smoker = 0,
                               Diabetes.Control = "NoDiabetes", CreaCl = 95,
                               Dialysis = 0, Hypertension = 0,
                               Peripheral.Vascular.Disease = 0,
                               Cerebrovascular.Accident = 0,
                               Cerebrovascular.Disease = 0,
                               Myocardial.Infarction = 0,
                               Congestive.Heart.Failure = 0,
                               Angina.Type = "0",
                               Afib.flutter = 0,
                               Number.of.Disease.Cornary.Vessels = "One",
                               Left.Main.Disease = 0,
                               Ejection.Fraction = 65,
                               Status = "Elective",
                               Dyslipidemia = 0,
                               Lipid.Lowering = 0,
                               Previous.Valve = 0,
                               Previous.Coronary.Artery.Bypass = 0,
                               Year.CABG = 2015)

patient_low_unexp <- data.frame(CT = 0, Age = 45, Gender = 0, Smoker = 0,
                              Diabetes.Control = "NoDiabetes", CreaCl = 95,
                              Dialysis = 0, Hypertension = 0,
                              Peripheral.Vascular.Disease = 0,
                              Cerebrovascular.Accident = 0,
                              Cerebrovascular.Disease = 0,
                              Myocardial.Infarction = 0,
                              Congestive.Heart.Failure = 0,
                              Angina.Type = "0",
                              Afib.flutter = 0,
                              Number.of.Disease.Cornary.Vessels = "One",
                              Left.Main.Disease = 0,
                              Ejection.Fraction = 65,
                              Status = "Elective",
                              Dyslipidemia = 0,
                              Lipid.Lowering = 0,
                              Previous.Valve = 0,
                              Previous.Coronary.Artery.Bypass = 0,
                              Year.CABG = 2015)
