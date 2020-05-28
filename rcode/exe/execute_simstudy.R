#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Execute script to run + process + summarise simulation study
#------------------------------------------------------------------------------#


# Load librairies + source code ----
#------------------------------------------------------------------------------#
source(file = "./rcode/packages/packages.R")
source(file = "./rcode/dgm/create_dirs.R")
source(file = "./rcode/sim/run_sim.R")
source(file = "./rcode/sumsim/summarize_simulation.R")
source(file = "./rcode/sumsim/diagnostics.R")

# Select datagen_scenarios and analysis_scenarios to be used
use_datagen_scenarios <- datagen_scenarios()[1:20,]
use_analysis_scenarios <- analysis_scenarios()
use_simulation_scenarios <- 1:nrow(use_datagen_scenarios)
rep <- 10

# Create filestructure  ----
#------------------------------------------------------------------------------#

dirpaths  <- lapply(create_dirpaths(use_analysis_scenarios = analysis_scenarios()),
       function(x) dir.create(x,recursive = TRUE))
filepaths <- lapply(create_dirpaths(use_analysis_scenarios = analysis_scenarios()),
                    FUN = function(x) create_filepaths(dirpaths = x,
                                                       use_datagen_scenarios = use_datagen_scenarios))

invisible(lapply(unlist(filepaths),
                 FUN= function(x) saveRDS(NULL,file=x)))


# Run simulation study  ----
#------------------------------------------------------------------------------#
# run_sim()
run_sim(rep = rep,
        use_datagen_scenarios = use_datagen_scenarios,
        use_analysis_scenarios = use_analysis_scenarios)



# Summarize processed simulation output  ----
#------------------------------------------------------------------------------#
# summarize_sim()
invisible(lapply(analysis_scenarios()[['method']],
                 FUN = function(x) sum_multiple_scenarios(use_simulation_scenarios = use_simulation_scenarios,
                                                          method = analysis_scenarios()[['method']][x],                         
                                                          pcutoff = 0.157,
                                                          estimator = 'MRR',
                                                          truth = 1)))

# Diagnostic checks  ----
#------------------------------------------------------------------------------#
# frequency exposure outcome
diagn_freqs(use_datagen_scenarios)

# warnings in model fitting
diagn_warnings()
