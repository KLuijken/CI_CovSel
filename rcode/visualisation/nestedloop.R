#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Visualisation simulation output
#------------------------------------------------------------------------------#

# Load libraries and data  ----
#------------------------------------------------------------------------------#
library(looplot)
library(ggplot2)

# Read in preferred data (plot in manuscript based on ....rds)
Comparison <- readRDS("./data/summarised/ML_0.157.rds")

# Select scenario numbers of preference (scenarios in manuscript ....)
use_scen_nums <- select_scenario_numbers()

# Pre-process data  ----
#------------------------------------------------------------------------------#

# Helper function for filtering data
plot_filter <- function(nobs, rhoL, sd_UY){
  plotinput <- input[input$nobs == nobs & input$rhoL == rhoL & input$sd_UY == sd_UY,]
  return(plotinput)
}



# Re-arrange data 
MSE_full       <- Comparison[Comparison$model=='Full','MSE_MRR']
MSE_selected   <- Comparison[Comparison$model=='Selected','MSE_MRR']
MSE_unadjusted <- Comparison[Comparison$model=='Unadjusted','MSE_MRR']
input          <- cbind(use_scen_nums, MSE_full, MSE_unadjusted)

input <- input[,c('bAL1','bAL2','bAL3','bYL1','bYL2','bYL3','MSE_full',
                  'MSE_selected','MSE_unadjusted')]

# Plot  ----
#------------------------------------------------------------------------------#

nested_loop_plot(input,
                 x = 'bAL1',
                 grid_rows = 'bYL1',
                 grid_cols = 'bYL2',
                 steps = c('bAL2','bAL3','bYL3')) # Make sure bYL3 is the outterstep
