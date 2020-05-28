#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Analyse data: estimate four models and store results
#------------------------------------------------------------------------------#


# Load libraries  ----
#------------------------------------------------------------------------------#
source("./rcode/analyses/helpers/be_firth.R")
source("./rcode/analyses/helpers/MRR.R")
source("./rcode/analyses/helpers/process_models.R")
source("./rcode/analyses/helpers/unadjusted.R")


# Analyse data  ----
#------------------------------------------------------------------------------#

# Output datatable with results of full and selected model for a single method
analyse_data <- function(analysis_scenario,
                         datagen_scenario,
                         data,
                         seed){
  scen_num <- datagen_scenario[['scen_num']]
  
  # Unadjusted analysis
  unadjusted <- analyse_unadjusted(data=data)
  
  # Store results of unadjusted
  results_unadj      <- data.table(scen_num, seed,"Unadjusted",
                                   analysis_scenario[['method']],
                                  t(unadjusted), t(rep(NA, times = (2 * ncol(data)))),
                                  "NULL", "NULL")
  
  # Estimate full model using maximum likelihood or FLIC, based on analysis scenario
  full <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L",(1:datagen_scenario[['nL']]))),collapse = "+")),
                               data = data,
                               firth = isTRUE(analysis_scenario[['method']] == "FLIC"),
                               pl = F))
  
  # Re-estimate intercept and store try_catch values and warnings in objects
  full <- pre_model(inputmodel = full,
                    datagen_scenario = datagen_scenario,
                    data = data)

  # Obtain marginal risk ratio and marginal odds ratio
  marginals_full    <- estimate_marginals(data =data,
                                          int = full$M_int$coefficients[1],
                                          modelcoefs = full$M$coefficients[-1])
  
  # Obtain model coefficients and standard errors
  coefficients_full <- obtain_coefficients(model = full$M,
                                           intmodel = full$M_int,
                                           datagen_scenario = datagen_scenario)
  
  # Obtain warnings model estimation
  warnings_full     <- obtain_warnings(premodel = full$preM,
                                       preintmodel = full$preM_int)
  
  # Store results of full model
  results_full      <- data.table(scen_num, seed,"Full",
                                  analysis_scenario[['method']],
                                  t(marginals_full),
                                  t(coefficients_full),
                                  t(warnings_full))
                                  
  # Use backward elimination on full model (either ML or FLIC)
  selected <- tryCatch.W.E(backwardf(object = full$M,
                                       slstay = analysis_scenario[['pcutoff']],
                                       trace = FALSE,
                                       scope = c(paste0("L",(1:datagen_scenario[['nL']]))),
                                      analysis_scenario = analysis_scenario)
                            )
   

   # Re-estimate intercept and store try_catch values and warnings in objects
   selected <- pre_model(inputmodel = selected,
                     datagen_scenario = datagen_scenario,
                     data = data)
   
   # Obtain marginal risk ratio and marginal odds ratio
   marginals_sel    <- estimate_marginals(data =data,
                                          int = selected$M_int$coefficients[1], 
                                          modelcoefs = selected$M$coefficients[-1])
   
   # Obtain model coefficients and standard errors
   coefficients_sel <- obtain_coefficients(model = selected$M,
                                           intmodel = selected$M_int, 
                                           datagen_scenario = datagen_scenario)
   
   # Obtain warnings model estimation
   warnings_sel     <- obtain_warnings(premodel = selected$preM,
                                       preintmodel = selected$preM_int)
   
   # Store results of selected model
   results_sel      <- data.table(scen_num, seed,
                                  paste0("Selected_",analysis_scenario[['pcutoff']]),
                                  analysis_scenario[['method']],
                                  t(marginals_sel),
                                  t(coefficients_sel),
                                  t(warnings_sel))

   # Set colnames equal
   colnames(results_unadj)  <- 
     colnames(results_full) <- 
     colnames(results_sel)  <- c("Scennum","Seed","Model","Method",
                              "MRR",
                              "MOR",
                              "(Intercept)",
                              names(data)[-1],
                              "se(Intercept)",
                              paste0("se(",names(data),")")[-1],
                              "Model warning",
                              "Intercept model warning")
  
  # Combine results in output matrix
  results <- rbind(results_unadj,results_full, results_sel)
  dir <- paste(analysis_scenario[['method']],
               analysis_scenario[['pcutoff']],
               sep="_")
  file <- paste0("S", scen_num)
  results <- cbind(results, filepath = paste0("./data/raw/",dir,"/", file, ".rds"))
  
  return(results)
}
