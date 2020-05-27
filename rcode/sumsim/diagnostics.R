#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Perfrom diagnostic checks of simulation results
#------------------------------------------------------------------------------#


# Simulated data ----
#------------------------------------------------------------------------------#
# Check whether exposure and outcome frequency are around expected value:

diagn_freqs_onemethod <- function(use_datagen_scenarios, method, pcutoff){
  frequency_A <- data.frame(matrix(NA, nrow = rep, ncol = nrow(use_datagen_scenarios)))
  colnames(frequency_A) <- paste0("freqA_S",1:nrow(use_datagen_scenarios))
  frequency_Y <- data.frame(matrix(NA, nrow = rep, ncol = nrow(use_datagen_scenarios)))
  colnames(frequency_Y) <- paste0("freqY_S",1:nrow(use_datagen_scenarios),
                                  "_Yint = ",use_datagen_scenarios[['Yint']])
  
  for(i in 1:nrow(use_datagen_scenarios)){
    results <- readRDS(paste0("./data/raw/",method,"_",pcutoff,"/S",i,".rds"))
    frequency_A[,i] <- results$freq_A[results$Model == "Unadjusted"]
    frequency_Y[,i] <- results$freq_Y[results$Model == "Unadjusted"]
    
  }
  
  print(dfSummary(frequency_A,
                  round.digits=3, varnumbers = F),
        file=file.path("./data/diagnostics",paste0(method,"_",pcutoff,"_FrequencyA.html")),
        silent = T)
  
  print(dfSummary(frequency_Y,
                  round.digits=3, varnumbers = F),
        file=file.path("./data/diagnostics",paste0(method,"_",pcutoff,"_FrequencyY.html")),
        silent = T)
}


diagn_freqs <- function(){
  apply(analysis_scenarios(),
        MARGIN = 1,
        FUN = function(x) diagn_freqs_onemethod(use_datagen_scenarios = use_datagen_scenarios,
                                                method = x[['method']],
                                                pcutoff = x[['pcutoff']]))
}


# Model warnings ----
#------------------------------------------------------------------------------#
