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

generate_overall_table <- function(method, pcutoff){
  # Load results
  raw_results <- readRDS(paste0("./data/summarised/",method,"_",pcutoff,".rds"))
  
  results <- data.table(scen_num = unique(raw_results[,scen_num]),
                        mse_full = raw_results[raw_results[,model]=="full",MSE_MRR],
                        mse_selected = raw_results[raw_results[,model]=="selected_0.157",MSE_MRR])
  
  results <- cbind(results,
                   bYA = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bYA"],
                   bAL1 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bAL1"],
                   bAL2 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bAL2"],
                   bAL3 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bAL3"],
                   bAL4 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bAL4"],
                   bYL1 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bYL1"],
                   bYL2 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bYL2"],
                   bYL3 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bYL3"],
                   bYL4 = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "bYL4"],
                   rhoL = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "rhoL"],
                   eventrate = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "eventrate"],
                   sd_UY = use_datagen_scenarios[use_datagen_scenarios$scen_num == results[,scen_num], "sd_UY"])
  
  # Store scenarios where selection results in lower mse
  mse_benefit <- results[results[,mse_full]<results[,mse_selected],]
  
  # Of the mse_benefit scenarios, how many have instruments in dgm?
  IVs <- mse_benefit[mse_benefit[,bYL1 == 0 & bAL1 != 0 | 
                                   bYL2 == 0 & bAL2 != 0 | 
                                   bYL3 == 0 & bAL3 != 0 | 
                                   bYL4 == 0 & bAL4 != 0],]
  
  # Of the mse_benefit scenarios, how many have noise variables in dgm?
  noise <- mse_benefit[mse_benefit[,bYL1 == 0 & bAL1 == 0 | 
                                     bYL2 == 0 & bAL2 == 0 | 
                                     bYL3 == 0 & bAL3 == 0 | 
                                     bYL4 == 0 & bAL4 == 0],]
  
  # Residual mse_benefit scenarios
  other <- mse_benefit[!(mse_benefit[,scen_num] %in% c(IVs[,scen_num],noise[,scen_num])),]
  
  table   <- data.table(`MSE(log(MRR)) BE < full` = nrow(mse_benefit),
                        `At least 3 IVs in dgm` = nrow(IVs),
                        `Similarity MSE(log(MRR)) between full and BE approach for scenarios with IVs`=
                          ifelse(nrow(IVs)!=0, 
                                 round(IVs[,mse_selected] / IVs[,mse_full], digits = 2),
                                 NA), # Need to think about similarity measure
                        `At least 3 noise vars in dgm` = nrow(noise),
                        `Similarity MSE(log(MRR)) between full and BE approach for scenarios with noise`=
                          ifelse(nrow(noise)!=0, 
                                 round(noise[,mse_selected] / noise[,mse_full], digits = 2),
                                 NA), 
                        `No IVs or noise in dgm` = nrow(mse_benefit) -  nrow(IVs) - nrow(noise),
                        `Similarity MSE(log(MRR)) between full and BE approach for other scenarios`=
                          ifelse(nrow(other)!=0, 
                                 round(other[,mse_selected] / other[,mse_full], digits = 2),
                                 NA))
  
  print(xtable(table),
        file=paste0("./results/tables/",method,"_",pcutoff,"_extensive_sim_overall.txt"))
  
}

