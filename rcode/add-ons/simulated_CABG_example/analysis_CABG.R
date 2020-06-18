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

# Perform analysis ----
#------------------------------------------------------------------------------#

# Generate data
CABG_data<- generate_data(N = 2266,
              betaTr.zero = FALSE,
              avsu = FALSE)

# Models
ML_full     <- glm(Postoperative.stroke ~ CT + Age + Gender + Smoker + 
                     Diabetes.Control + CreaCl+ Dialysis + Hypertension + 
                     Peripheral.Vascular.Disease + Cerebrovascular.Accident + 
                     Cerebrovascular.Disease + Myocardial.Infarction + 
                     Congestive.Heart.Failure + Angina.Type + Afib.flutter + 
                     Number.of.Diseased.Coronary.Vessels + Left.Main.Disease +
                     Ejection.Fraction + Status + Dyslipidemia + Lipid.Lowering +
                     Previous.Valve + Previous.Coronary.Artery.Bypass +
                     Year.CABG,
                 data = CABG_data,
                 family = binomial)

ML_selected <- step(ML_full, 
                    scope = list(lower = Postoperative.stroke ~ CT),
                    direction = "backward",
                    k = 2,
                    trace = FALSE)

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
                       firth = TRUE,
                       pl = TRUE)

# Re-estimate the intercept
pred_intercept       <- t(Firth_full$coefficients[-1]) %*% 
                        data.matrix(CABG_data[,!(colnames(CABG_data) %in% "Postoperative.stroke")])
intercept_Firth_full <- glm(CABG_data$Postoperative.stroke ~ offset(pred_intercept), 
                            family = binomial) 

FLIC_allcoefficients <- c(Intercept_Firth_full$coefficients[1],
                          Firth_full$coefficients[-1])

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

# Conditional odds ratio, full model
cOR_full_ML_Est <- summary(ML_full)$coefficients["CT","Estimate"]
cOR_full_ML_Low <- summary(ML_full)$coefficients["CT","Estimate"] - 1.96 *
  summary(ML_full)$coefficients["CT","Std. Error"]
cOR_full_ML_Up  <- summary(ML_full)$coefficients["CT","Estimate"] + 1.96 *
  summary(ML_full)$coefficients["CT","Std. Error"]

cOR_full_Firth_Est <- summary(Firth_full, print = F)$coefficients["CT"]
cOR_full_Firth_Low <- summary(Firth_full, print = F)$ci.lower["CT"]
cOR_full_Firth_Up  <- summary(Firth_full, print = F)$ci.upper["CT"]
