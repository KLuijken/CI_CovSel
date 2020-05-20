#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Summarize simulation results
#------------------------------------------------------------------------------#
# Output: bias, empirical SE, empirical variance, MSE for MRR or MOR


# Libraries  ----
#------------------------------------------------------------------------------#
source("./rcode/analyses/helpers/select_scenarios.R")

# Helpers  ----
#------------------------------------------------------------------------------#
Bias <- function(results, 
                 estimator, # string
                 truth){
  bias <- c(by(results, results[['model']],
             function(x) mean(log(unlist(x[[estimator]])))-log(truth)))
  
  return(bias)
}

EmpSE <- function(results,
                  estimator){
  empSE <- c(by(results,results[['model']], function(x) 
    sqrt(sum((log(unlist(x[[estimator]])) - mean(log(unlist(x[[estimator]]))))^2)/(nrow(x)-1))))
  
  return(empSE)
}

EmpVar <- function(results,
                   estimator){
  empvar <- c(by(results,results[['model']], function(x) 
    sum((log(unlist(x[[estimator]])) - mean(log(unlist(x[[estimator]]))))^2)/(nrow(x)-1)))
  
  return(empvar)
}

MSE <- function(results, 
                estimator, # string
                truth){
  mse <- c(by(results, results[['model']],
              function(x) mean((log(unlist(x[[estimator]]))-log(truth))^2)))
  
  return(mse)
}



# Summarize simulations  ----
#------------------------------------------------------------------------------#

# Function to summarize a single scenario. 
sum_one_scenario <- function(scen_num, method, pcutoff,
                estimator,
                truth){
  results <- readRDS(paste0("./data/results_raw/",method,"_",pcutoff,"/S",scen_num,".rds"))

  bias <- Bias(results = results, estimator = estimator, truth = truth)
  empSE <- EmpSE(results = results, estimator = estimator)
  empVar <- EmpVar(results = results, estimator = estimator)
  mse <- MSE(results = results, estimator = estimator, truth = truth)
  
  out <- data.table(unique(results[['Scennum']]),names(bias),
                    paste0(method,"_",pcutoff), cbind(bias,empSE,empVar,mse))
  colnames(out) <- c("scen_num","model","method",
                     paste0(c("bias_","empSE_","empvar_","MSE_"),estimator))
  
  return(out)
}

# test, remove this when all works
# sum_one_scenario(scen_num = 481,
#                  method = "ML",
#                  pcutoff = 0.157,
#                  estimator = 'MRR',
#                  truth = 1)



# Workhorse to summarize results from multiple scenarios. Output stored as data.table
# in .rds file in ./data/summarised
sum_multiple_scenarios <- function(use_simulation_scenarios,
                                   method,
                                   pcutoff,
                                   estimator, truth){
  output <- do.call(rbind, lapply(use_simulation_scenarios,
                                  FUN = function(x) sum_one_scenario(
                                    scen_num = use_simulation_scenarios[x],
                                                    method = method,
                                                    pcutoff = pcutoff,
                                                    estimator = estimator,
                                                    truth = truth)))
  
  saveRDS(output,file = paste0("./data/summarised/",method,"_",pcutoff,".rds"))
}

# test, remove this when all works

# sum_multiple_scenarios(select_scenario_numbers()[c(1:10)],
#                        method = "ML",
#                        pcutoff = 0.157,
#                        estimator = 'MRR',
#                        truth = 1)


