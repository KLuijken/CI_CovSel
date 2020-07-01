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
cOR_full_Firth_Est <- round( exp( coef( Firth_full)["CT"]), digits = 2)
cOR_full_Firth_Low <- round( exp( confint(Firth_full)[2,1]), digits=2)
cOR_full_Firth_Up  <- round( exp( confint(Firth_full)[2,2]), digits=2)
paste0(cOR_full_Firth_Est,
       "(95% CI, ", cOR_full_Firth_Low,
       "; ", cOR_full_Firth_Up,")")

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

### NB: confidence intervals are bootstrapped below, such that they can be combined
### with bootstrap estimation of confidence intervals of the selected model

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
cOR_selected_Firth_Est <- round( exp( coef( Firth_selected)["CT"]), digits = 2)
cOR_selected_Firth_Low <- round( exp( confint(Firth_selected)[2,1]), digits=2)
cOR_selected_Firth_Up  <- round( exp( confint(Firth_selected)[2,2]), digits=2)
paste0(cOR_selected_Firth_Est,
       "(95% CI, ", cOR_selected_Firth_Low,
       "; ", cOR_selected_Firth_Up,")")


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
b_rep <- 1
seeds <- 1:b_rep
coefficients_full <-
  coefficients_selected <- matrix(0,
                                  ncol = length(coef(Firth_full)),
                                  nrow = b_rep,
                               dimnames = list(NULL, names(coef(Firth_full))))

CI_marginal_OR_full     <- 
  CI_marginal_OR_selected <- matrix(NA, nrow = b_rep, ncol = 1)

# Warning message:
#   In if (working$df == 1 | working$df == mat[mat[, 3] == max(mat[,  :
#     the condition has length > 1 and only the first element will be used

for(i in 1:b_rep){
  # Sample bootstrap data using prespecified seeds
  set.seed(seeds[i])
  # Stratified sampling, events
  events            <- CABG_data[CABG_data$Postoperative.stroke ==1,]
  sampled_events    <- sample(1:nrow(events),
                              size = nrow(events), replace = TRUE)
  bs_events         <- events[sampled_events,]
  # Sample non-events
  nonevents         <- CABG_data[CABG_data$Postoperative.stroke == 0,]
  sampled_nonevents <- sample(1:nrow(nonevents),
                              size = nrow(nonevents), replace = TRUE)
  bs_nonevents      <- nonevents[sampled_nonevents,]
  
  # Complete bootstrap sample
  bs_sample         <- rbind(bs_events,bs_nonevents)
  
  # Remove variables with no variance from the full model
  drop <- colnames(bs_sample)[sapply(bs_sample, function(x) ifelse(class(x) == "factor",
                                                                   all(duplicated(x)[-1L]),
                                                                   var(x)==0))]
  bs_sample <- bs_sample[,!(colnames(bs_sample) %in% drop)]
  
  #### Full model: estimate Firth model in bootstrap sample
  Firth_full_bs   <- logistf(as.formula(paste0("bs_sample$Postoperative.stroke ~ ",
                                               paste(colnames(bs_sample[,2:ncol(bs_sample)]),
                                                     collapse = "+"))),
                             data = bs_sample,
                             firth = TRUE, # set firth = FALSE for maximum likelihood estimation
                             pl = FALSE)   # do not compute profile likelihood confidence intervals
  
  
  # Re-estimate intercept in bootrap sample
  pred_intercept_bs       <- Firth_full_bs$linear.predictors
  intercept_Firth_full_bs <- glm(bs_sample$Postoperative.stroke ~ offset(pred_intercept_bs), 
                                 family = binomial)$coefficients[1]
  
  # Store coefficients
  coefficients_full[,"(Intercept)"]  <- intercept_Firth_full_bs
  coefficients_full[i, names(
    coef(Firth_full_bs)[-1])]        <- coef(Firth_full_bs)[-1]
  
  # Estimate MOR in bootrap sample
  mod_mat <- model.matrix(as.formula(
    Firth_full_bs$call[["formula"]]),
    data = bs_sample)
  
  mod_mat[,"CT"] <- 0
  PredA0 <- plogis(as.vector(mod_mat %*% Firth_full_bs$coefficients + intercept_Firth_full_bs))
  mod_mat[,"CT"] <- 1
  PredA1 <- plogis(as.vector(mod_mat %*% Firth_full_bs$coefficients + intercept_Firth_full_bs))
  
  CI_marginal_OR_full[i,] <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))
  
  
  # #### Selected model: perform backward elimination in the bootstrap sample
  # Firth_selected_bs      <- backward(Firth_full_bs,
  #                                     scope = colnames(bs_sample)[
  #                                       !(colnames(bs_sample) %in% c("Postoperative.stroke","CT"))],
  #                                     slstay = 0.157,
  #                                     trace = FALSE,
  #                                     data = bs_sample,
  #                                     pl = F,
  #                                     analysis_scenario = NULL)
  # 
  # # Re-estimate intercept in bootrap sample
  # pred_intercept_bs           <- Firth_selected_bs$linear.predictors
  # intercept_Firth_selected_bs <- glm(bs_sample$Postoperative.stroke ~ offset(pred_intercept_bs),
  #                                    family = binomial)$coefficients[1]
  # 
  # # Store coefficients
  # coefficients_selected[,"(Intercept)"]  <- intercept_Firth_selected_bs
  # coefficients_selected[i, names(
  #   coef(Firth_selected_bs)[-1])]        <- coef(Firth_selected_bs)[-1]
  # 
  # # Estimate MOR in bootrap sample
  # mod_mat <- model.matrix(as.formula(
  #   Firth_selected_bs$call[["formula"]]),
  #   data = bs_sample)
  # 
  # mod_mat[,"CT"] <- 0
  # PredA0 <- plogis(as.vector(mod_mat %*% Firth_selected_bs$coefficients + intercept_Firth_selected_bs))
  # mod_mat[,"CT"] <- 1
  # PredA1 <- plogis(as.vector(mod_mat %*% Firth_selected_bs$coefficients + intercept_Firth_selected_bs))
  # 
  # CI_marginal_OR_selected[i,] <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))


}

# Obtain MOR + CI
MOR_full_Firth_Low <- round(quantile(CI_marginal_OR_full, 0.025, na.rm = T), digits = 2)
MOR_full_Firth_Up  <- round(quantile(CI_marginal_OR_full, 0.975, na.rm = T), digits = 2)
paste0(round(MOR_full_Firth_Est,digits=2),
       "(95% CI, ", MOR_full_Firth_Low,
       "; ", MOR_full_Firth_Up,")")

# Obtain number of variables
boot_01_full <- (coefficients_full != 0) * 1
boot_inclusion_full <- apply(boot_01_full, 2, function(x) sum(x) / length(x) * 100)

# Obtain MOR + CI model with backward elimination
MOR_selected_Firth_Low <- round(quantile(CI_marginal_OR_selected, 0.025), digits = 2)
MOR_selected_Firth_Up  <- round(quantile(CI_marginal_OR_selected, 0.975), digits = 2)
paste0(round(MOR_selected_Firth_Est,digits=2),
       "(95% CI, ", MOR_selected_Firth_Low,
       "; ", MOR_selected_Firth_Up,")")

# Obtain number of variables
boot_01_selected <- (coefficients_selected != 0) * 1
boot_inclusion_selected <- apply(boot_01_selected, 2, function(x) sum(x) / length(x) * 100)


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
