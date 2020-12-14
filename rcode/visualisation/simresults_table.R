#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Create table simulation results
#------------------------------------------------------------------------------#


# Load libraries + source code ----
#------------------------------------------------------------------------------#


# Helpers to generate results table ----
#------------------------------------------------------------------------------#

generate_overall_table <- function(method, pcutoff, estimator, use_datagen_scenarios){
  # Load results
  raw_results <- readRDS(paste0("./data/summarised/",
                                method, "_",
                                pcutoff, "_",
                                estimator, ".rds"))
  
  results <- data.table(scen_num = unique(raw_results[,scen_num]),
                        mse_full = raw_results[raw_results[,model]=="full", get(paste0("MSE_",estimator))],
                        mse_selected = raw_results[raw_results[,model]=="selected_0.157", get(paste0("MSE_",estimator))],
                        bias_full = raw_results[raw_results[,model]=="full", get(paste0("bias_",estimator))],
                        bias_selected = raw_results[raw_results[,model]=="selected_0.157", get(paste0("bias_",estimator))])
  
  # add datagenerating parameters to outcome frame
  results <- cbind(results,
                use_datagen_scenarios)
  
  # Store scenarios where selection results in lower mse
  mse_benefit <- results[results[,mse_selected]<results[,mse_full],]
  
  # Of the mse_benefit scenarios, how many have instruments in dgm?
  IVs <- mse_benefit[mse_benefit[, (bYL1 == 0 & bAL1 != 0) | 
                                   (bYL2 == 0 & bAL2 != 0) | 
                                   (bYL3 == 0 & bAL3 != 0) | 
                                   (bYL4 == 0 & bAL4 != 0)],]
  
  # Of the mse_benefit scenarios, how many have noise variables in dgm?
  noise <- mse_benefit[mse_benefit[,(bYL1 == 0 & bAL1 == 0) | 
                                     (bYL2 == 0 & bAL2 == 0) | 
                                     (bYL3 == 0 & bAL3 == 0) | 
                                     (bYL4 == 0 & bAL4 == 0)],]
  
  # Residual mse_benefit scenarios
  other <- mse_benefit[!(mse_benefit[,scen_num] %in% IVs[,scen_num] | mse_benefit[,scen_num] %in% noise[,scen_num]),]
  
  table   <- data.table(bias_full = mean(results[,bias_full]),
                        bias_selected = mean(results[,bias_selected]),
                        `Relative efficiency MSE(log(MRR)) BE approach`= median(results[,mse_selected]/results[,mse_full]), 
                        `Min relative efficiency MSE(log(MRR)) BE approach`= min(results[,mse_selected]/results[,mse_full]), 
                        `Max relative efficiency MSE(log(MRR)) BE approach`= max(results[,mse_selected]/results[,mse_full]), 
                        `Number of scenarios MSE(log(MRR)) BE < full` = nrow(mse_benefit),
                        `At least 3 IVs in dgm` = nrow(IVs),
                        `At least 3 noise vars in dgm` = nrow(noise),
                        `No IVs or noise in dgm` = nrow(other)
  )
  
  print(xtable(table),
        file=paste0("./results/tables/",method,"_",pcutoff,"_",estimator,"_extensive_sim_overall.txt"),
        include.rownames = FALSE)
  
}


generate_stratified_table <- function(method, pcutoff, estimator, use_datagen_scenarios){
  # Load results
  raw_results <- readRDS(paste0("./data/summarised/",
                                method, "_",
                                pcutoff, "_",
                                estimator, ".rds"))
  
  results <- data.table(scen_num = unique(raw_results[,scen_num]),
                        mse_full = raw_results[raw_results[,model]=="full", get(paste0("MSE_",estimator))],
                        mse_selected = raw_results[raw_results[,model]=="selected_0.157", get(paste0("MSE_",estimator))],
                        bias_full = raw_results[raw_results[,model]=="full", get(paste0("bias_",estimator))],
                        bias_selected = raw_results[raw_results[,model]=="selected_0.157", get(paste0("bias_",estimator))])
  
  results <- cbind(results,
                   use_datagen_scenarios)
  
  # Stratify the table
  levels_bYA <- unique(use_datagen_scenarios[,"bYA"])
  levels_eventrate <- unique(use_datagen_scenarios[,"eventrate"])
  levels_nevents <- unique(use_datagen_scenarios[,"nevents"])
  
  stratify <- expand.grid(nevents = levels_nevents, eventrate = levels_eventrate, bYA = levels_bYA)
  
  stratify <- data.table(bYA = stratify[,"bYA"],
                         eventrate = stratify[,"eventrate"],
                         nevents = stratify[,"nevents"])
  
  table <- data.frame(matrix(NA, nrow = nrow(stratify), ncol = 9))
  colnames(table) <- c("bias_full",
                       "bias_selected",
                       "Relative efficiency MSE(log(MRR)) BE approach", 
                       "Min relative efficiency MSE(log(MRR)) BE approach", 
                       "Max relative efficiency MSE(log(MRR)) BE approach", 
                       "Number of scenarios MSE(log(MRR)) BE < full",
                       "At least 3 IVs in dgm",
                       "At least 3 noise vars in dgm",
                       "No IVs or noise in dgm")
  
  for(i in 1:nrow(stratify)){
    # Store scenarios where selection results in lower mse
    mse_benefit <- results[results[, bYA == stratify[i,bYA] &
                                     eventrate == stratify[i,eventrate] &
                                     nevents == stratify[i,nevents] &
                                     mse_selected < mse_full],]
    
    # Of the mse_benefit scenarios, how many have instruments in dgm?
    IVs <- mse_benefit[mse_benefit[, bYA == stratify[i,bYA] &
                                     eventrate == stratify[i,eventrate] &
                                     nevents == stratify[i,nevents] & ( 
                                       bYL1 == 0 & bAL1 != 0 | 
                                         bYL2 == 0 & bAL2 != 0 | 
                                         bYL3 == 0 & bAL3 != 0 | 
                                         bYL4 == 0 & bAL4 != 0)],]
    
    # Of the mse_benefit scenarios, how many have noise variables in dgm?
    noise <- mse_benefit[mse_benefit[, bYA == stratify[i,bYA] &
                                       eventrate == stratify[i,eventrate] &
                                       nevents == stratify[i,nevents] & (
                                         bYL1 == 0 & bAL1 == 0 | 
                                           bYL2 == 0 & bAL2 == 0 | 
                                           bYL3 == 0 & bAL3 == 0 | 
                                           bYL4 == 0 & bAL4 == 0)],]
    
    # Residual mse_benefit scenarios
    other <- mse_benefit[!(mse_benefit[,scen_num] %in% IVs[,scen_num] | mse_benefit[,scen_num] %in% noise[,scen_num]),]
    
    # Fill final output table
    table[i,"bias_full"] <- mean(results[results[, bYA == stratify[i,bYA] &
                                                   eventrate == stratify[i,eventrate] &
                                                   nevents == stratify[i,nevents]],bias_full])
    table[i,"bias_selected"]<- mean(results[results[, bYA == stratify[i,bYA] &
                                                      eventrate == stratify[i,eventrate] &
                                                      nevents == stratify[i,nevents]],bias_selected])
    table[i,"Relative efficiency MSE(log(MRR)) BE approach"] <- median(results[results[, bYA == stratify[i,bYA] &
                                                                                         eventrate == stratify[i,eventrate] &
                                                                                         nevents == stratify[i,nevents]],mse_selected]/
                                                                         results[results[, bYA == stratify[i,bYA] &
                                                                                           eventrate == stratify[i,eventrate] &
                                                                                           nevents == stratify[i,nevents]],mse_full]) 
    table[i,"Min relative efficiency MSE(log(MRR)) BE approach"] <- min(results[results[, bYA == stratify[i,bYA] &
                                                                                          eventrate == stratify[i,eventrate] &
                                                                                          nevents == stratify[i,nevents]],mse_selected]/
                                                                          results[results[, bYA == stratify[i,bYA] &
                                                                                            eventrate == stratify[i,eventrate] &
                                                                                            nevents == stratify[i,nevents]],mse_full]) 
    table[i,"Max relative efficiency MSE(log(MRR)) BE approach"] <- max(results[results[, bYA == stratify[i,bYA] &
                                                                                          eventrate == stratify[i,eventrate] &
                                                                                          nevents == stratify[i,nevents]],mse_selected]/
                                                                          results[results[, bYA == stratify[i,bYA] &
                                                                                            eventrate == stratify[i,eventrate] &
                                                                                            nevents == stratify[i,nevents]],mse_full])
    table[i,"Number of scenarios MSE(log(MRR)) BE < full"] <- nrow(mse_benefit)
    table[i,"At least 3 IVs in dgm"] <- nrow(IVs)
    table[i,"At least 3 noise vars in dgm"] <- nrow(noise)
    table[i,"No IVs or noise in dgm"] <- nrow(other)
  }
  
  output <- cbind(stratify,table)
  
  print(xtable(output),
        file=paste0("./results/tables/",method,"_",pcutoff,"_",estimator,"_extensive_sim_stratified.txt"),
        include.rownames = FALSE)
  
}
