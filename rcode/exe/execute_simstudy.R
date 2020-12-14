#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Workhorse to execute (run + process + summarise) simulation study
#------------------------------------------------------------------------------#
# Generate 1000 datasets for each of the 3960 scenarios, perform analyses described
# in manuscript and specified in ./rcode/dgm/sim_scen.R and ./rcode/dgm/analyse_data.R
# Save output per scenario under ./data/raw, processed data in ./data/summarised
# and tables and figures in ./results.
# Details on the functionalities and output can be found in helper scripts.

# Load librairies + source code ----
#------------------------------------------------------------------------------#
source(file = "./rcode/packages/packages.R")
source(file = "./rcode/dgm/create_dirs.R")
source(file = "./rcode/sim/run_sim.R")
source(file = "./rcode/sumsim/summarize_simulation.R")
source(file = "./rcode/sumsim/diagnostics.R")
source(file = "./rcode/visualisation/simresults_table.R")

# Select datagen_scenarios and analysis_scenarios to be used
all_datagen_scenarios <- datagen_scenarios()  # Redundant step, used for testing
use_datagen_scenarios <- all_datagen_scenarios
use_analysis_scenarios <- analysis_scenarios()
rep <- 1000

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



# Diagnostic checks  ----
#------------------------------------------------------------------------------#
# frequency exposure outcome
diagn_freqs(use_datagen_scenarios)

# warnings in model fitting
diagn_warnings()



# Summarize processed simulation output  ----
#------------------------------------------------------------------------------#
# summarize_sim()
invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) sum_multiple_scenarios(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  use_datagen_scenarios = use_datagen_scenarios,
                  estimator = "MRR",
                  rep = rep)))

invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) sum_multiple_scenarios(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  use_datagen_scenarios = use_datagen_scenarios,
                  estimator = "A",
                  rep = rep)))



# Generate output tables  ----
#------------------------------------------------------------------------------#
invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) generate_overall_table(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  estimator = "MRR",
                  use_datagen_scenarios = use_datagen_scenarios)))

invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) generate_stratified_table(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  estimator = "MRR",
                  use_datagen_scenarios = use_datagen_scenarios)))


invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) generate_overall_table(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  estimator = "A",
                  use_datagen_scenarios = use_datagen_scenarios)))

invisible(apply(use_analysis_scenarios,
                MARGIN = 1,
                FUN = function(x) generate_stratified_table(
                  method = x[['method']],
                  pcutoff = x[['pcutoff']],
                  estimator = "A",
                  use_datagen_scenarios = use_datagen_scenarios)))
