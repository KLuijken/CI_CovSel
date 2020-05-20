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
                                  NA,
                                  t(unadjusted), t(rep(NA, times = (2 * ncol(data)))),
                                  freq_Y = sum(data$Y)/nrow(data),
                                  freq_A = sum(data$A)/nrow(data))
  
  # Estimate full model using maximum likelihood or FLIC, based on analysis scenario
  full <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L",(1:datagen_scenario[['nL']]))),collapse = "+")),
                               data = data,
                               firth = isTRUE(analysis_scenario[['method']] == "FLIC"),
                               pl = F))
  
  # Obtain warnings about model estimation
  full <- pre_model(inputmodel = full,
                    datagen_scenario = datagen_scenario,
                    data = data)

  # Obtain marginal risk ratio and marginal odds ratio
  marginals_full    <- unlist(estimate_marginals(data =data,
                                          int = full$M_int$coefficients[1],
                                          modelcoefs = full$M$coefficients[-1]))
  
  # Obtain model coefficients and standard errors
  coefficients_full <- unlist(obtain_coefficients(model = full$M,
                                           intmodel = full$M_int,
                                           datagen_scenario = datagen_scenario))
  
  # Store results of full model
  results_full      <- data.table(scen_num, seed,"Full",
                                  analysis_scenario[['method']],
                                  t(marginals_full), t(coefficients_full), NA, NA)
  
  # # Use backward elimination on full model (either ML or FLIC)
  # selected <- tryCatch.W.E(backwardf(object = full$M,
  #                                     sls = analysis_scenario[['pcutoff']],
  #                                     trace = FALSE,
  #                                     scope = c(paste0("L",(1:datagen_scenario[['nL']])))))
  # 
  # # Obtain warnings about model estimation
  # selected <- pre_model(inputmodel = selected,
  #                   datagen_scenario = datagen_scenario,
  #                   data = data)
  # 
  # # Obtain marginal risk ratio and marginal odds ratio
  # marginals_sel    <- unlist(estimate_marginals(data =data,
  #                                        int = selected$M_int$coefficients[1], 
  #                                        modelcoefs = selected$M$coefficients[-1]))
  # 
  # # Obtain model coefficients and standard errors
  # coefficients_sel <- unlist(obtain_coefficients(model = selected$M,
  #                                         intmodel = selected$M_int, 
  #                                         datagen_scenario = datagen_scenario))
  # 
  # # Store results of selected model
  # results_sel      <- data.table(scen_num, seed,
  #                                paste0("Selected_",analysis_scenario[['pcutoff']]),
  #                                analysis_scenario[['method']],
  #                                t(marginals_sel),
  #                                t(coefficients_sel), NA, NA)

  # Set colnames equal
  # colnames(results_unadj)  <- 
  #   colnames(results_full) <- 
  #   colnames(results_sel)  <- c("Scennum","Seed","Model","Method",
  #                            "MRR",
  #                            "MOR",
  #                            "(Intercept)",
  #                            names(data)[-1],
  #                            "se(Intercept)",
  #                            paste0("se(",names(data),")")[-1],
  #                            "freq_Y",
  #                            "freq_A")

  colnames(results_unadj)  <-
    colnames(results_full) <- c("scennum","seed","model","method",
                             "MRR",
                             "MOR",
                             "(Intercept)",
                             names(data)[-1],
                             "se(Intercept)",
                             paste0("se(",names(data),")")[-1],
                             "freq_Y",
                             "freq_A")
  
  # Combine results in output matrix
  #results <- rbind(results_unadj,results_full, results_sel)
  results <- rbind(results_unadj,results_full)
  dir <- paste(analysis_scenario[['method']],
               analysis_scenario[['pcutoff']],
               sep="_")
  file <- paste0("S",unique(results[["Scennum"]]))
  results <- cbind(results, filepath = paste0("./data/raw/",dir,"/",file,".rds"))
  
  return(results)
}