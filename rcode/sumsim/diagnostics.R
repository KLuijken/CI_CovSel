#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Perfrom diagnostic checks of simulation results
#------------------------------------------------------------------------------#


# Simulated data ----
#------------------------------------------------------------------------------#
# Check whether exposure and outcome frequency are around expected value:

diagn_freqs <- function(use_datagen_scenarios){
  frequency_A <- data.frame(matrix(NA, nrow = rep, ncol = nrow(use_datagen_scenarios)))
  colnames(frequency_A) <- paste0("freqA_S",use_datagen_scenarios[['scen_num']])
  frequency_Y <- data.frame(matrix(NA, nrow = rep, ncol = nrow(use_datagen_scenarios)))
  colnames(frequency_Y) <- paste0("freqY_S",use_datagen_scenarios[['scen_num']],
                                  "_Yint = ",use_datagen_scenarios[['Yint']])
  
  for(i in 1:nrow(use_datagen_scenarios)){
    results <- readRDS(paste0("./data/raw/descriptives/S",use_datagen_scenarios[i,"scen_num"],".rds"))
    frequency_A[,i] <- results$freq_A
    frequency_Y[,i] <- results$freq_Y
    
  }
  
  # Create output directory
  dir.create(file.path(".","data","diagnostics"), recursive = TRUE)
  
  # Save diagnostics
  print(dfSummary(frequency_A,
                  round.digits=3, varnumbers = F),
        file=file.path("./data/diagnostics/FrequencyA.html"),
        silent = T)
  
  print(dfSummary(frequency_Y,
                  round.digits=3, varnumbers = F),
        file=file.path("./data/diagnostics/FrequencyY.html"),
        silent = T)
}


# Model warnings ----
#------------------------------------------------------------------------------#

diagn_warnings_onemethod <- function(use_datagen_scenarios, method, pcutoff){
  for(i in 1:nrow(use_datagen_scenarios)){
    results <- readRDS(paste0("./data/raw/",method,"_",pcutoff,"/S",i,".rds"))
    
    # Store model warnings
    if(!identical(results$mod_warning[!is.na(results$mod_warning)],
                             character(0))){
      write.table(paste0(results$mod_warning[!is.na(results$mod_warning)],
                         " in Scenario ", 
                         results$scen_num[!is.na(results$mod_warning)],
                         " Seed ", 
                         results$seed[!is.na(results$mod_warning)]),
                  file = "./data/diagnostics/warningsmodel.txt", 
                  append = T,
                  row.names = F,
                  col.names = F)} 
    
    # Store intercept model warnings
      if(!identical(results$intmod_warning[!is.na(results$intmod_warning)],
                             character(0))){
      write.table(paste0(results$intmod_warning[!is.na(results$intmod_warning)],
                         " in Scenario ", 
                         results$scen_num[!is.na(results$intmod_warning)],
                         " Seed ", 
                         results$seed[!is.na(results$intmod_warning)]),
                  file = "./data/diagnostics/warningsintmodel.txt", 
                  append = T, 
                  row.names = F,
                  col.names = F)} 
    
  }
}

diagn_warnings <- function(){
  apply(analysis_scenarios(),
        MARGIN = 1,
        FUN = function(x) diagn_warnings_onemethod(use_datagen_scenarios = use_datagen_scenarios,
                                                method = x[['method']],
                                                pcutoff = x[['pcutoff']]))
}

