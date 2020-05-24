#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper function to perform backward elimination after FLIC estimation
#------------------------------------------------------------------------------#


backwardf <- function (object, scope, steps = 1000, slstay, trace = TRUE, 
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
  while (istep < steps & working$df >= 2){
    istep <- istep + 1
    mat <- drop1(working)
    inscope <- match(scope, rownames(mat))  
    inscope <- inscope[!is.na(inscope)]
    if (all(mat[inscope, 3] < slstay)){      ### Bug-fix from original backward() function in logistf package
                                            ### added inscope to subset the p-values
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
