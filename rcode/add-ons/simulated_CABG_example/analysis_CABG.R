#------------------------------------------------------------------------------#
# Manuscript: Causal inference when selection of confounders is partly based on 
#   backward elimination: likely biased, not often more efficient
# Authors: K Luijken, R H H Groenwold, M van Smeden, S Strohmaier, G Heinze
# Author script: K Luijken
#
# Perform CABG example analysis
#------------------------------------------------------------------------------#

# This script describes the analysis on simulated data of the example CABG study
# presented in manuscript section 2. The following analyses are presented:
## conditional odds ratio of a full model, estimated using FLIC
## marginal odds ratio of a full model, estimated using predicted potential 
#     outcomes of a FLIC model
## conditional odds ratio of a full model, estimated using FLIC and backward 
#     elimination
## marginal odds ratio of a full model, estimated using predicted potential 
#     outcomes of a FLIC model and backward elimination
## confidence intervals for marginal odds ratios are obtained using bootstrap

# All analyses accompanying the manuscript can be found on https://github.com/KLuijken/CI_CovSel

# Load librairies ----
#------------------------------------------------------------------------------#
devtools::install_github("georgheinze/logistf", ref="develop")
library(logistf)


# Generate simulated data ----
#------------------------------------------------------------------------------#
source(file = "./rcode/add-ons/simulated_CABG_example/simulate_data_Gregorich.R")
set.seed(20200618)
CABG_data<- generate_data(N = 2266,
              betaTr.zero = FALSE,
              avsu = FALSE)
CABG_data$Age <- scale(CABG_data$Age) * 0.5
CABG_data$CreaCl <- scale(CABG_data$CreaCl) * 0.5
CABG_data$Ejection.Fraction <- scale(CABG_data$Ejection.Fraction) * 0.5
CABG_data$Year.CABG <- scale(CABG_data$Year.CABG) * 0.5

# Conditional OR, FLIC, full model ----
#------------------------------------------------------------------------------#

# Estimate model
Firth_full   <- logistf(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                         Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                         Peripheral.Vascular.Disease + Cerebrovascular.Accident+ 
                         Cerebrovascular.Disease + Myocardial.Infarction + 
                         Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                         Number.of.Diseased.Coronary.Vessels + 
                         Left.Main.Disease + Ejection.Fraction + Status + 
                         Dyslipidemia + Lipid.Lowering + Previous.Valve + 
                         Previous.Coronary.Artery.Bypass + Year.CABG,
                       data = CABG_data,
                       control = logistf.control(maxit = 200, maxstep = 5),
                       firth = TRUE, # set firth = FALSE for max.likelihood est.
                       pl = TRUE,    # compute profile likelihood CIs
                       flic = TRUE)  # Perform intercept correction

# Obtain coefficients + CI
cOR_full_Firth_Est <- round( exp( coef( Firth_full)["CT"]), digits = 2)
cOR_full_Firth_Low <- round( exp( confint(Firth_full)[2,1]), digits=2)
cOR_full_Firth_Up  <- round( exp( confint(Firth_full)[2,2]), digits=2)
paste0(cOR_full_Firth_Est,
       "(95% CI, ", cOR_full_Firth_Low,
       "; ", cOR_full_Firth_Up,")")

# Marginal OR, FLIC, full model ----
#------------------------------------------------------------------------------#

# Create potential outcome datasets
All_unexposed <- CABG_data
All_unexposed$CT <- 0
All_exposed <- CABG_data
All_exposed$CT <- 1

# Obtain predicted potential outcomes
PredA0 <- as.vector(predict(Firth_full, 
                            newdata = All_unexposed, type = "response"))
PredA1 <- as.vector(predict(Firth_full, 
                            newdata = All_exposed, type = "response"))

# Estimate Marginal OR
MOR_full_Firth_Est <- 
  (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))

### NB: CIs are bootstrapped below, such that they can be combined with 
### bootstrap estimation of CIs of the selected model

# Conditional OR, FLIC, backward elimination ----
#------------------------------------------------------------------------------#

# Estimate model
Firth_selected  <- backward(Firth_full,
                            scope = c("Age", "Gender", "Smoker", 
                                      "Diabetes.Control", "CreaCl", "Dialysis", 
                                      "Hypertension", 
                                      "Peripheral.Vascular.Disease",
                                      "Cerebrovascular.Accident", 
                                      "Cerebrovascular.Disease", 
                                      "Myocardial.Infarction", 
                                      "Congestive.Heart.Failure", "Angina.Type", 
                                      "Afib.flutter", 
                                      "Number.of.Diseased.Coronary.Vessels",
                                      "Left.Main.Disease", "Ejection.Fraction",
                                      "Status", "Dyslipidemia","Lipid.Lowering",
                                      "Previous.Valve", 
                                      "Previous.Coronary.Artery.Bypass",
                                      "Year.CABG"),
                            slstay = 0.157,
                            trace = TRUE,
                            control = logistf.control(maxit = 200, maxstep = 5),
                            pl = TRUE)
Firth_selected      <- flic(Firth_selected)

# Obtain coefficient + CI
cOR_selected_Firth_Est <- round( exp( coef( Firth_selected)["CT"]), digits = 2)
cOR_selected_Firth_Low <- round( exp( confint(Firth_selected)[2,1]), digits=2)
cOR_selected_Firth_Up  <- round( exp( confint(Firth_selected)[2,2]), digits=2)
paste0(cOR_selected_Firth_Est,
       "(95% CI, ", cOR_selected_Firth_Low,
       "; ", cOR_selected_Firth_Up,")")


# Marginal OR, FLIC, backward elimination ----
#------------------------------------------------------------------------------#

# Obtain predicted potential outcomes
PredA0 <- as.vector(predict(Firth_selected, 
                            newdata = All_unexposed, type = "response"))
PredA1 <- as.vector(predict(Firth_selected, 
                            newdata = All_exposed, type = "response"))

# Estimate Marginal OR
MOR_selected_Firth_Est <- 
  (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))


# Bootstrap confidence interval ----
#------------------------------------------------------------------------------#

b_rep <- 10
seeds <- 1:b_rep
coefficients_full <-
  coefficients_selected <- matrix(NA,
                                  ncol = length(coef(Firth_full)),
                                  nrow = b_rep,
                               dimnames = list(NULL, names(coef(Firth_full))))

CI_marginal_OR_full     <- 
  CI_marginal_OR_selected <- matrix(NA, nrow = b_rep, ncol = 1)

# Run for 500 bootstrap resamplings
for(i in 1:b_rep){
  
  # Sample bootstrap data using prespecified seeds
  set.seed(seeds[i])
  
  # Bootstrap sampling ----
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
  drop <- colnames(bs_sample)[sapply(bs_sample, function(x) 
    ifelse(class(x) == "factor",
           all(duplicated(x)[-1L]),
           var(x)==0))]
  bs_sample <- bs_sample[,!(colnames(bs_sample) %in% drop)]
  
  # Create potential outcome datasets
  All_unexposed_bs <- bs_sample
  All_unexposed_bs$CT <- 0
  All_exposed_bs <- bs_sample
  All_exposed_bs$CT <- 1
  
  # Full model ----
  # Estimate Firth model in bootstrap sample
  Firth_full_bs <- logistf(as.formula(paste0("Postoperative.stroke ~ ",
                                  paste(colnames(bs_sample[,2:ncol(bs_sample)]),
                                  collapse = "+"))),
                             data = bs_sample,
                             control = logistf.control(maxit = 200,maxstep = 5),
                             firth = TRUE, # firth = FALSE for ML estimation
                             pl = FALSE,   # no profile likelihood CIs
                             flic = TRUE)

  # Store coefficients
  coefficients_full[i, names(coef(Firth_full_bs))] <- coef(Firth_full_bs)
  
  # Obtain predicted potential outcomes
  PredA0 <- as.vector(predict(Firth_full_bs, 
                              newdata = All_unexposed_bs, type = "response"))
  PredA1 <- as.vector(predict(Firth_full_bs, 
                              newdata = All_exposed_bs, type = "response"))
  
  # Estimate MOR in bootstrap sample
  CI_marginal_OR_full[i,] <- 
    (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))
  
  
  # Selected model ----
  # Perform backward elimination in the bootstrap sample
  Firth_selected_bs <- backward(Firth_full_bs,
                                scope = colnames(bs_sample)[
                                  !(colnames(bs_sample) %in% 
                                      c("Postoperative.stroke","CT",drop))],
                                slstay = 0.157,
                                trace = FALSE,
                                pl = FALSE,    # no profile likelihood CIs
                                data = bs_sample,
                                control = logistf.control(maxit = 200,
                                                          maxstep = 5))
  
  # Re-estimate intercept in bootstrap sample
  Firth_selected_bs      <- flic(Firth_selected_bs)

  # Store coefficients
  coefficients_selected[i, names(coef(Firth_selected_bs))] <- 
    coef(Firth_selected_bs)

  # Obtain predicted potential outcomes
  PredA0 <- as.vector(predict(Firth_selected_bs, 
                              newdata = All_unexposed_bs, type = "response"))
  PredA1 <- as.vector(predict(Firth_selected_bs, 
                              newdata = All_exposed_bs, type = "response"))
  
  # Estimate MOR in bootrap sample
  CI_marginal_OR_selected[i,] <- 
    (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))


}



# Summarize results
# Obtain MOR + CI
MOR_full_Firth_Low <- round(quantile(CI_marginal_OR_full, 0.025, na.rm = T),
                            digits = 2)
MOR_full_Firth_Up  <- round(quantile(CI_marginal_OR_full, 0.975, na.rm = T),
                            digits = 2)
paste0(round(MOR_full_Firth_Est,digits=2),
       "(95% CI, ", MOR_full_Firth_Low,
       "; ", MOR_full_Firth_Up,")")

# Obtain number of variables: stability measures
boot_01_full <- (coefficients_full != 0) * 1
boot_inclusion_full <- apply(boot_01_full, 2, function(x) sum(x, na.rm=T)/length(x)*100)

# Obtain MOR + CI model with backward elimination
MOR_selected_Firth_Low <- round(quantile(CI_marginal_OR_selected, 0.025),
                                digits = 2)
MOR_selected_Firth_Up  <- round(quantile(CI_marginal_OR_selected, 0.975), 
                                digits = 2)
paste0(round(MOR_selected_Firth_Est,digits=2),
       "(95% CI, ", MOR_selected_Firth_Low,
       "; ", MOR_selected_Firth_Up,")")

# Obtain number of variables: stability measures
boot_01_selected <- (coefficients_selected != 0) * 1
boot_inclusion_selected <- apply(boot_01_selected, 2, 
                                 function(x) sum(x, na.rm=T)/length(x)*100)

