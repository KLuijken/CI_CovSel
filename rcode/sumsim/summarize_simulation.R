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
                 true_value){
  bias       <- c(by(results, results[['model']],
             function(x) mean(log(unlist(x[[estimator]])), na.rm = TRUE)-log(true_value)))
  
  return(bias)
}

EmpSE <- function(results,
                  estimator){
  count <- sum(!is.na(results[[estimator]]))
  empSE <- c(by(results,results[['model']], function(x) 
    sqrt(sum((log(unlist(x[[estimator]])) - mean(log(unlist(x[[estimator]])), na.rm = TRUE))^2)/(count-1))))
  
  return(empSE)
}

EmpVar <- function(results,
                   estimator){
  count <- sum(!is.na(results[[estimator]]))
  empvar <- c(by(results,results[['model']], function(x) 
    sum((log(unlist(x[[estimator]])) - mean(log(unlist(x[[estimator]])), na.rm = TRUE))^2)/(count-1)))
  
  return(empvar)
}

MSE <- function(results, 
                estimator,
                true_value){
  mse        <- c(by(results, results[['model']],
              function(x) mean((log(unlist(x[[estimator]]))-log(true_value))^2, na.rm = TRUE)))
  
  return(mse)
}

na_count <- function(results, 
                     estimator){
  na_count <- c(by(results,results[['model']],
                   function(x) sum(is.na(x[[estimator]]))))
}



# Summarize simulations  ----
#------------------------------------------------------------------------------#

# Function to summarize a single scenario. 
sum_one_scenario <- function(datagen_scenario, method, pcutoff,
                estimator){
  scen_num   <- datagen_scenario[['scen_num']]
  
  # load results
  results    <- readRDS(paste0("./data/raw/",method,"_",pcutoff,"/S",scen_num,".rds"))
  
  # store true value
  ifelse(estimator == "COR",{
    true_value <- datagen_scenario[['bYA']]
  },{
    truth      <- readRDS(paste0("./data/raw/truth/S",scen_num,".rds"))
    true_value <- truth[[estimator]]
  })
  
  bias <- Bias(results = results, estimator = estimator, true_value = true_value)
  empSE <- EmpSE(results = results, estimator = estimator)
  empVar <- EmpVar(results = results, estimator = estimator)
  mse <- MSE(results = results, estimator = estimator, true_value = true_value)
  na_count <- na_count(results = results, estimator = estimator)

  out <- data.table(scen_num, names(bias),
                    paste0(method,"_",pcutoff), cbind(bias, empSE, empVar, mse, na_count))
  colnames(out) <- c("scen_num", "model", "method",
                     paste0(c("bias_","empSE_","empvar_","MSE_"),estimator),"removed_replications")

  return(out)
}


# Workhorse to summarize results from multiple scenarios. Output stored as data.table
# in .rds file in ./data/summarised
sum_multiple_scenarios <- function(method,
                                   pcutoff,
                                   use_datagen_scenarios,
                                   estimator){
  output <- do.call(rbind, apply(use_datagen_scenarios,
                                 MARGIN = 1,
                                 FUN = function(x) sum_one_scenario(
                                    scen_num = x[['scen_num']],
                                                    method = method,
                                                    pcutoff = pcutoff,
                                                    estimator = estimator)))
  
  dir.create(file.path(".","data","summarised"), recursive = TRUE)
  saveRDS(output,file = paste0("./data/summarised/",method,"_",pcutoff,".rds"))
}

