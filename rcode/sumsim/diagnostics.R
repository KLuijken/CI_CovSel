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
    results <- readRDS(paste0("./data/raw/descriptives/S", i,".rds"))
    frequency_A[,i] <- results$freq_A
    frequency_Y[,i] <- results$freq_Y
    
  }
  
  ### Overall frequencies
  hist_A <- unlist(c(frequency_A))
  pdf("./data/diagnostics/FrequencyA_hist.pdf")
  hist(hist_A, xlab = "frequency A")
  dev.off()
  
  hist_Y_0.2 <- unlist(c(frequency_Y[which(use_datagen_scenarios$eventrate == 0.2)]))
  hist_Y_0.03 <- unlist(c(frequency_Y[which(use_datagen_scenarios$eventrate == 0.03)]))
  pdf("./data/diagnostics/FrequencyY_hist.pdf")
  hist(hist_Y_0.2, xlab = "frequency Y")
  hist(hist_Y_0.03, xlab = "frequency Y")
  dev.off()
  
  
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
  write.table("warnings",
              file = paste0("./data/diagnostics/warningsmodel_",method,".txt"), 
              append = F,
              row.names = F,
              col.names = F)
  
  for(i in 1:nrow(use_datagen_scenarios)){
    results <- readRDS(paste0("./data/raw/",method,"_",pcutoff,"/S",
                              use_datagen_scenarios[i,"scen_num"],".rds"))
    
    # Store model warnings
    if(!identical(results$mod_warning[!is.na(results$mod_warning)],
                             logical(0))){
      write.table(paste0(results$mod_warning[!is.na(results$mod_warning)],
                         " in Scenario ", 
                         results$scen_num[!is.na(results$mod_warning)],
                         " Seed ", 
                         results$seed[!is.na(results$mod_warning)]),
                  file = paste0("./data/diagnostics/warningsmodel_",method,".txt"), 
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

