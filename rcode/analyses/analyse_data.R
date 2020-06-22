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
  # Unadjusted analysis
  unadjusted <- analyse_unadjusted(data=data)
  
  # Store results of unadjusted
  results_unadj      <- data.table(scen_num = datagen_scenario[['scen_num']], 
                                   seed = seed, 
                                   model = "unadjusted",
                                   method = analysis_scenario[['method']],
                                   t(unadjusted),
                                   mod_warning = NA, intmod_warning = NA)
  
  # Estimate full model using maximum likelihood or FLIC, based on analysis scenario
  full <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L",(1:datagen_scenario[['nL']]))),collapse = "+")),
                               data = data,
                               firth = isTRUE(analysis_scenario[['method']] == "FLIC"),
                               pl = F))
  
  # Re-estimate intercept and store try_catch values and warnings in objects
  full <- pre_model(inputmodel = full,
                    datagen_scenario = datagen_scenario,
                    data = data)

  # Obtain warnings model estimation
  warnings_full     <- obtain_warnings(premodel = full$preM,
                                       preintmodel = full$preM_int)
  
  # Obtain marginal risk ratio and marginal odds ratio
  marginals_full    <- estimate_marginals(warning = warnings_full,
                                          data =data,
                                          int = full$M_int$coefficients[1],
                                          modelcoefs = full$M$coefficients[-1])
  
  # Obtain model coefficients and standard errors
  coefficients_full <- obtain_coefficients(warning = warnings_full,
                                           model = full$M,
                                           intmodel = full$M_int,
                                           datagen_scenario = datagen_scenario)
  
  # Store results of full model
  results_full      <- data.table(datagen_scenario[['scen_num']],
                                  seed,
                                  "full",
                                  analysis_scenario[['method']],
                                  t(marginals_full),
                                  t(coefficients_full),
                                  t(warnings_full))
                                  
  # Use backward elimination on full model (either ML or FLIC)
  selected <- tryCatch.W.E(backwardf(object = full$M,
                                       slstay = analysis_scenario[['pcutoff']],
                                       trace = FALSE,
                                       scope = c(paste0("L",(1:datagen_scenario[['nL']]))),
                                       analysis_scenario = analysis_scenario,
                                       pl = F)
                            )

   # Re-estimate intercept and store try_catch values and warnings in objects
   selected <- pre_model(inputmodel = selected,
                     datagen_scenario = datagen_scenario,
                     data = data)
   
   # Obtain warnings model estimation
   warnings_sel     <- obtain_warnings(premodel = selected$preM,
                                       preintmodel = selected$preM_int)
   
   # Obtain marginal risk ratio and marginal odds ratio
   marginals_sel    <- estimate_marginals(warning = warnings_sel,
                                          data =data,
                                          int = selected$M_int$coefficients[1], 
                                          modelcoefs = selected$M$coefficients[-1])
   
   # Obtain model coefficients and standard errors
   coefficients_sel <- obtain_coefficients(warning = warnings_sel,
                                           model = selected$M,
                                           intmodel = selected$M_int, 
                                           datagen_scenario = datagen_scenario)
   
   # Store results of selected model
   results_sel      <- data.table(datagen_scenario[['scen_num']], 
                                  seed,
                                  paste0("selected_",analysis_scenario[['pcutoff']]),
                                  analysis_scenario[['method']],
                                  t(marginals_sel),
                                  t(coefficients_sel),
                                  t(warnings_sel))

   # Set colnames equal
   colnames(results_full) <- 
     colnames(results_sel)  <- c("scen_num","seed","model","method",
                              "MRR",
                              "MOR",
                              "(Intercept)",
                              names(data)[-1],
                              "se(Intercept)",
                              paste0("se(",names(data),")")[-1],
                              "mod_warning",
                              "intmod_warning")
  
  # Combine results in output matrix
  results <- rbind(results_unadj,results_full, results_sel, fill=T)
  dir <- paste(analysis_scenario[['method']],
               analysis_scenario[['pcutoff']],
               sep="_")
  file <- paste0("S", datagen_scenario[['scen_num']])
  results <- cbind(results, filepath = paste0("./data/raw/",dir,"/", file, ".rds"))
  
  return(results)
}
