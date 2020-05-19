#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Helper for summarizing simulation results
#------------------------------------------------------------------------------#

# Function to select scenario numbers based on parameter values. If not specified,
# all values for that parameter are selected.

select_scenario_numbers <- function(nobs = NA,
                                    nL = NA,
                                    bAL1 = NA,
                                    bAL2 = NA,
                                    bAL3 = NA,
                                    bYL1 = NA,
                                    bYL2 = NA,
                                    bYL3 = NA,
                                    sd_UY = NA,
                                    rhoL = NA){
  selection <- datagen_scenarios()[ifelse(is.na(nobs),TRUE,datagen_scenarios()[["nobs"]]%in% nobs) &
                                     ifelse(is.na(nL),TRUE,datagen_scenarios()[["nL"]]%in% nL) &
                                     ifelse(is.na(bAL1),TRUE,datagen_scenarios()[["bAL1"]]%in% bAL1) &
                                     ifelse(is.na(bAL2),TRUE,datagen_scenarios()[["bAL2"]]%in% bAL2) &
                                     ifelse(is.na(bAL3),TRUE,datagen_scenarios()[["bAL3"]]%in% bAL3) &
                                     ifelse(is.na(bYL1),TRUE,datagen_scenarios()[["bYL1"]]%in% bYL1) &
                                     ifelse(is.na(bYL2),TRUE,datagen_scenarios()[["bYL2"]]%in% bYL2) &
                                     ifelse(is.na(bYL3),TRUE,datagen_scenarios()[["bYL3"]]%in% bYL3) &
                                     ifelse(is.na(sd_UY),TRUE,datagen_scenarios()[["sd_UY"]]%in% sd_UY) &
                                     ifelse(is.na(rhoL),TRUE,datagen_scenarios()[["rhoL"]]%in% rhoL),]
  selected_scennums <- unique(selection[['scen_num']])
  
  return(selected_scennums)
}