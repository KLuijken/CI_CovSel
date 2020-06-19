#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Perform CABG example analysis
#------------------------------------------------------------------------------#


# Load librairies + source code ----
#------------------------------------------------------------------------------#
library(logistf)

source(file = "./rcode/add-ons/simulated_CABG_example/simulate_data_Gregorich.R")

# Helper function ----
#------------------------------------------------------------------------------#

backwardf <- function (object, scope, steps = 1000, slstay, trace = TRUE, 
                       printwork = FALSE, ...){
  data <-  object$data
  
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
  while (istep < steps & working$df >= 2){
    istep <- istep + 1
    mat <- drop1(working)
    inscope <- match(scope, rownames(mat))  
    inscope <- inscope[!is.na(inscope)]
    if (all(mat[inscope, 3] < slstay)){      
      break}
    removal <- rownames(mat)[mat[, 3] == max(mat[inscope, 
                                                 3])]
    newform = as.formula(paste("~.-", removal))
    if (working$df == 2 | working$df == mat[mat[, 3] == max(mat[, 3]), 2]){
      working <- update(working, formula = newform, pl = FALSE)}else{ 
        working <- update(working, formula = newform)}
    if (trace) {
      cat("drop1:\n")       
      print(mat)           
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

# Generate simulated data ----
#------------------------------------------------------------------------------#

set.seed(20200618)
CABG_data<- generate_data(N = 2266,
              betaTr.zero = FALSE,
              avsu = FALSE)

# Conditional OR ----
#------------------------------------------------------------------------------#

# --- Full model --- #

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
paste0(cOR_full_Firth_Est, "(95% CI, ", cOR_full_Firth_Low,"; ", cOR_full_Firth_Up,")")

# --- Selected model --- #

# Estimate model
Firth_selected      <- backwardf(Firth_full,
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



# Marginal OR ----
#------------------------------------------------------------------------------#

# --- Full model --- #

# Use Firth_full model previous step.
# Re-estimate the intercept
pred_intercept       <- Firth_full$linear.predictors
intercept_Firth_full <- glm(CABG_data$Postoperative.stroke ~ offset(pred_intercept), 
                            family = binomial)$coefficients[1]

# Obtain predicted potential outcomes
mod_mat <- model.matrix(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                          Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                          Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                          Cerebrovascular.Disease + Myocardial.Infarction + 
                          Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                          Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                          Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                          Previous.Valve + Previous.Coronary.Artery.Bypass +
                          Year.CABG,
                        data = CABG_data)

mod_mat[,"CT"] <- 0
PredA0 <- plogis(as.vector(test %*% Firth_full$coefficients + intercept_Firth_full))
mod_mat[,"CT"] <- 1
PredA1 <- plogis(as.vector(test %*% Firth_full$coefficients + intercept_Firth_full))

# Estimate Marginal OR
marginal_OR <- (mean(PredA1) * (1- mean(PredA0)))/((1-mean(PredA1)) * mean(PredA0))

# Bootstrap confidence interval:


# Conditional RR ----
#------------------------------------------------------------------------------#

# --- Full model --- #

# Estimate model
Poisson_full   <- glm(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                          Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                          Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                          Cerebrovascular.Disease + Myocardial.Infarction + 
                          Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                          Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                          Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                          Previous.Valve + Previous.Coronary.Artery.Bypass +
                          Year.CABG,
                        data = CABG_data,
                        family = poisson)

# Obtain coefficient + CI
cRR_full_Poisson_Est <- round( exp( summary( Poisson_full)$coefficients["CT", "Estimate"]), digits=2)
cRR_full_Poisson_Low <- round( exp( summary( Poisson_full)$coefficients["CT", "Estimate"] -
                                      1.96*summary( Poisson_full)$coefficients["CT", "Std. Error"]), digits=2)
cRR_full_Poisson_Up  <- round( exp( summary( Poisson_full)$coefficients["CT", "Estimate"] +
                                      1.96*summary( Poisson_full)$coefficients["CT", "Std. Error"]), digits=2)
paste0(cRR_full_Poisson_Est, "(95% CI, ", cRR_full_Poisson_Low,"; ", cRR_full_Poisson_Up,")")
