#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Summarize simulation results
#------------------------------------------------------------------------------#
# Output: bias, empirical SE, empirical variance, MSE for MRR or MOR


# Libraries  ----
#------------------------------------------------------------------------------#

# Helpers  ----
#------------------------------------------------------------------------------#
# In the helper functions, we select na.rm = TRUE to deal with the NAs that are 
# introduced because in 10 datasets the effect could not be estimated using 
# Maximum Likelihood.

Value <- function(results, 
                  estimator, # string
                  transformation){
  bias       <- c(by(results, results[['model']],
                     function(x) mean(transformation(unlist(x[[estimator]])), na.rm = TRUE)))
  
  return(bias)
}

Bias <- function(results, 
                 estimator, # string
                 true_value,
                 transformation){
  bias       <- c(by(results, results[['model']],
             function(x) mean(transformation(unlist(x[[estimator]])), na.rm = TRUE)-transformation(true_value)))
  
  return(bias)
}

EmpSE <- function(results,
                  estimator,
                  transformation,
                  rep){
  empSE <- c(by(results,results[['model']], function(x) 
    sqrt(sum((transformation(unlist(x[[estimator]])) - mean(transformation(unlist(x[[estimator]])), na.rm = TRUE))^2, na.rm = TRUE)/(rep-1))))
  
  return(empSE)
}

EmpVar <- function(results,
                   estimator,
                   transformation,
                   rep){
  empvar <- c(by(results,results[['model']], function(x) 
    sum((transformation(unlist(x[[estimator]])) - mean(transformation(unlist(x[[estimator]])), na.rm = TRUE))^2, na.rm = TRUE)/(rep-1)))
  
  return(empvar)
}

MSE <- function(results, 
                estimator,
                true_value,
                transformation){
  mse        <- c(by(results, results[['model']],
              function(x) mean((transformation(unlist(x[[estimator]]))-transformation(true_value))^2, na.rm = TRUE)))
  
  return(mse)
}

count_missings <- function(results, 
                     estimator){
  count <- c(by(results,results[['model']],
                   function(x) sum(is.na(x[[estimator]]))))
  return(count)
}

count_warnings <- function(results){
  count <- c(by(results,results[['model']],
                function(x) sum(!is.na(x[['mod_warning']]))))
  return(count)
}

show_warnings <- function(results){
  paste(by(results,results[['model']],
                function(x) unique(x[['mod_warning']])))
}


# Summarize simulations  ----
#------------------------------------------------------------------------------#

# Function to summarize a single scenario. 
sum_one_scenario <- function(datagen_scenario, method, pcutoff,
                estimator, rep){
  scen_num   <- datagen_scenario[['scen_num']]
  
  # load results
  results    <- readRDS(paste0("./data/raw/",method,"_",pcutoff,"/S",scen_num,".rds"))
  
  # store true value and corresponding transformation
  ifelse(estimator == "A",{
    true_value <- datagen_scenario[['bYA']]
    trans      <- identity
  },{
    truth      <- readRDS(paste0("./data/raw/truth/S",scen_num,".rds"))
    true_value <- truth[[estimator]]
    trans      <- log
  })
  
  value <- Value(results = results, estimator = estimator, transformation = trans)
  bias <- Bias(results = results, estimator = estimator, true_value = true_value, transformation = trans)
  empSE <- EmpSE(results = results, estimator = estimator, rep = rep, transformation = trans)
  empVar <- EmpVar(results = results, estimator = estimator, rep = rep, transformation = trans)
  mse <- MSE(results = results, estimator = estimator, true_value = true_value, transformation = trans)
  na_count <- count_missings(results = results, estimator = estimator)
  war_count <- count_warnings(results = results)
  war_message <- show_warnings(results = results)

  out <- data.table(scen_num, names(bias),
                    paste0(method,"_",pcutoff),
                    value,
                    bias,
                    empSE,
                    empVar,
                    mse,
                    na_count,
                    war_count,
                    war_message)
  
  colnames(out) <- c("scen_num", "model", "method",
                     paste0(c("value_","bias_","empSE_","empvar_","MSE_"),estimator),
                     "removed_replications",
                     "warnings_count",
                     "warnings_message")

  return(out)
}


# Workhorse to summarize results from multiple scenarios. Output stored as data.table
# in .rds file in ./data/summarised
sum_multiple_scenarios <- function(method,
                                   pcutoff,
                                   use_datagen_scenarios,
                                   estimator, 
                                   rep){
  output <- do.call(rbind, apply(use_datagen_scenarios,
                                 MARGIN = 1,
                                 sum_one_scenario,
                                  method = method,
                                  pcutoff = pcutoff,
                                  estimator = estimator,
                                  rep = rep))
  
  dir.create(file.path(".","data","summarised"), recursive = TRUE)
  saveRDS(output,file = paste0("./data/summarised/",method,"_",pcutoff,"_",estimator,".rds"))
}

