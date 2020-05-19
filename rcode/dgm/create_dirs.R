#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Create directory structure to store sim results, run in exe
#------------------------------------------------------------------------------#



# Helpers ----
#------------------------------------------------------------------------------#


# Function to create a separate directory per method

create_dirpaths <- function(use_analysis_scenarios){
  levels <- paste(unique(use_analysis_scenarios[['method']]),
                  unique(use_analysis_scenarios[['pcutoff']]),
                  sep="_")
  dirpaths <- paste0("./data/results_raw/",levels)
  return(dirpaths)
}




# Function to create empty rds files to speed up saving during simruns

create_filepaths <- function(dirpaths, use_datagen_scenarios){
  filename <- paste0("S",unique(use_datagen_scenarios[['scen_num']]),".rds")
  filepath <- paste0(dirpaths,"/",filename)
  
  return(filepath)
}



