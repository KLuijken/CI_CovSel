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


# Analyse data  ----
#------------------------------------------------------------------------------#

# Output datatable with results of full and selected model for a single method
analyse_data <- function(analysis_scenario,
                               datagen_scenario,
                               df,
                               seed){
  # Unadjusted model ----
  # Estimate unadjusted model using maximum likelihood or FLIC, based on analysis scenario
  unadjusted        <- tryCatch.W.E(logistf(Y ~ A,
                               data = df,
                               pl = FALSE,
                               control = logistf.control(maxit = 50, maxstep = 2),
                               firth = isTRUE(analysis_scenario[['method']] == "FLIC"),
                               flic = TRUE))
  
  # If model fitting succeeds, continue, else, in case simple error occurs in 
  # model fitting, store error and NA for all values
  if(class(unadjusted$value)[1] == "logistf"){
    # Obtain warnings model estimation
    warnings_unadj    <- obtain_warnings(premodel = unadjusted)
    
    # Obtain marginal risk ratio and marginal odds ratio
    marginals_unadj   <- estimate_marginals(data = df,
                                    modelcoefs = unadjusted$value$coefficients)
    
    # Obtain model coefficients and standard errors
    coefs_unadj       <- obtain_coefficients(model = unadjusted$value,
                                    data = df,
                                    datagen_scenario = datagen_scenario)
    
    # Store results of unadjusted
    results_unadj     <- data.table(datagen_scenario[['scen_num']],
                                    seed,
                                    "unadjusted",
                                    analysis_scenario[['method']],
                                    t(marginals_unadj),
                                    t(coefs_unadj),
                                    warnings_unadj)}else{
    # Store results in case simple error occurs
    results_unadj     <- data.table(datagen_scenario[['scen_num']],
                                    seed,
                                    "unadjusted",
                                    analysis_scenario[['method']],
                                    t(rep(NA, times = ((datagen_scenario[['nL']] + 2)*2) + 3)), # empty values for coefficients (int, A, L, R_squared), se of coefficients (int, A, L, R_squared), and marginals
                                    mod_warning = paste(unadjusted$value)) # SimpleErrors are stored in $value
                                    }
  
  # Full model ----
  # Estimate full model using maximum likelihood or FLIC, based on analysis scenario
  full              <- tryCatch.W.E(logistf(as.formula(paste(c("Y~A" ,paste0("L",(1:datagen_scenario[['nL']]))),collapse = "+")),
                               data = df,
                               pl = FALSE,
                               control = logistf.control(maxit = 50, maxstep = 2),
                               firth = isTRUE(analysis_scenario[['method']] == "FLIC"),
                               flic = TRUE))
  
  # If model fitting succeeds, continue, else, in case simple error occurs in 
  # model fitting, store error and NA for all values
  if(class(full$value)[1] == "logistf"){
    # Obtain warnings model estimation
    warnings_full     <- obtain_warnings(premodel = full)
    
    # Obtain marginal risk ratio and marginal odds ratio
    marginals_full    <- estimate_marginals(data = df,
                                    modelcoefs = full$value$coefficients)
    
    # Obtain model coefficients and standard errors
    coefs_full        <- obtain_coefficients(model = full$value,
                                    data = df,
                                    datagen_scenario = datagen_scenario)
    
    # Store results of full model
    results_full      <- data.table(datagen_scenario[['scen_num']],
                                    seed,
                                    "full",
                                    analysis_scenario[['method']],
                                    t(marginals_full),
                                    t(coefs_full),
                                    warnings_full)}else{
    # Store results in case simple error occurs
    results_full      <- data.table(datagen_scenario[['scen_num']],
                                    seed,
                                    "full",
                                    analysis_scenario[['method']],
                                    t(rep(NA, times = ((datagen_scenario[['nL']] + 2)*2) + 3)),  # empty values for coefficients (int, A, L, R_squared), se of coefficients (int, A, L, R_squared), and marginals
                                    mod_warning = paste(full$value)) # SimpleErrors are stored in $value
                                    }

  
  # Selected model ----
  # Use backward elimination on full model (either ML or FLIC)
  selected          <- tryCatch.W.E(backwardf(object = full$value,
                               slstay = as.numeric(analysis_scenario[['pcutoff']]),
                               trace = FALSE,
                               scope = c(paste0("L",(1:datagen_scenario[['nL']]))),
                               analysis_scenario = analysis_scenario))
  
  # If model fitting succeeds, continue, else, in case simple error occurs in 
  # model fitting, store error and NA for all values
  if(class(selected$value)[1] == "logistf"){
    # Obtain warnings model estimation
    warnings_sel     <- obtain_warnings(premodel = selected)
    
    # Apply flic
    selected_flic    <- flic(selected$value)
    
    # Obtain marginal risk ratio and marginal odds ratio
    marginals_sel    <- estimate_marginals(data = df,
                                    modelcoefs = selected_flic$coefficients)
    
    # Obtain model coefficients and standard errors
    # NOTE: this is suboptimal, as this does not save the se(Intercept) after flic,
    # but the function "flic()" is not working properly yet since logistf package update
    # and intercept coefficient output is not the main focus.
    coefs_sel        <- obtain_coefficients(model = selected$value,
                                    data = df, 
                                    datagen_scenario = datagen_scenario)
    coefs_sel["(Intercept)"] <- selected_flic$coefficients[1]
    
    # Store results of selected model
    results_sel      <- data.table(datagen_scenario[['scen_num']], 
                                    seed,
                                    paste0("selected_",analysis_scenario[['pcutoff']]),
                                    analysis_scenario[['method']],
                                    t(marginals_sel),
                                    t(coefs_sel),
                                    warnings_sel)}else{
    # Store results in case simple error occurs
    results_sel      <- data.table(datagen_scenario[['scen_num']],
                                    seed,
                                    paste0("selected_",analysis_scenario[['pcutoff']]),
                                    analysis_scenario[['method']],
                                    t(rep(NA, times = ((datagen_scenario[['nL']] + 2)*2) + 3)),  # empty values for coefficients (int, A, L, R_squared), se of coefficients (int, A, L, R_squared), and marginals
                                    mod_warning = paste(selected$value)) # SimpleErrors are stored in $value 
                                    }
  
  
  
  # Set colnames equal
  colnames(results_unadj) <- 
    colnames(results_full) <- 
    colnames(results_sel)  <- c("scen_num","seed","model","method",
                                "MRR",
                                "MOR",
                                "(Intercept)",
                                names(df)[-1],
                                "se(Intercept)",
                                paste0("se(",names(df),")")[-1],
                                "R_squared",
                                "mod_warning")
  
  # Combine results in output matrix
  results <- rbind(results_unadj,results_full, results_sel, fill=T)
  dir <- paste(analysis_scenario[['method']],
               analysis_scenario[['pcutoff']],
               sep="_")
  file <- paste0("S", datagen_scenario[['scen_num']])
  results <- cbind(results, filepath = paste0("./data/raw/",dir,"/", file, ".rds"))
  
  return(results)
}
