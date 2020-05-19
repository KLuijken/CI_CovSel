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

# Select datagen_scenarios and analysis_scenarios to be used
use_datagen_scenarios <- datagen_scenarios()[1:10,]


# Create filestructure  ----
#------------------------------------------------------------------------------#

dirpaths  <- lapply(create_dirpaths(use_analysis_scenarios = analysis_scenarios()),
       dir.create)
filepaths <- lapply(dirpaths,
                    FUN = function(x) create_filepaths(dirpaths = dirpaths[x],
                                                       use_datagen_scenarios = use_datagen_scenarios))
invisible(lapply(unlist(filepaths),
                 FUN= function(x) saveRDS(NULL,file=x)))


# Run simulation study  ----
#------------------------------------------------------------------------------#
# run_sim()
run_sim(rep = 5,
        use_datagen_scenarios = use_datagen_scenarios,
        use_analysis_scenarios = use_analysis_scenarios)



# Summarize processed simulation output  ----
#------------------------------------------------------------------------------#
# summarize_sim()
invisible(lapply(analysis_scenarios()[['method']],
                 FUN = function(x) sum_multiple_scenarios(use_simulation_scenarios = select_scenario_numbers(),
                                                          method = analysis_scenarios()[['method']][x],                         
                                                          pcutoff = 0.157,
                                                          estimator = 'MRR',
                                                          truth = 1)))