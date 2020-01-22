####################################################
# Covariate selection project
#
# December 2019 
# Author: K Luijken
# File content: simulation output evaluation and summary
####################################################

allfilesS1 <- list.files(paste0(filepath,"/Simoutput/S1"), pattern = ".rds")
setwd(paste0(filepath,"/Simoutput/S1"))
allS1 <- do.call(list, lapply(1:length(allfilesS1),FUN= function(x) readRDS(allfilesS1[x])))

allfilesS2 <- list.files(paste0(filepath,"/Simoutput/S2"), pattern = ".rds")
setwd(paste0(filepath,"/Simoutput/S2"))
allS2 <- do.call(list, lapply(1:length(allfilesS2),FUN= function(x) readRDS(allfilesS2[x])))

allfilesS3 <- list.files(paste0(filepath,"/Simoutput/S3"), pattern = ".rds")
setwd(paste0(filepath,"/Simoutput/S3"))
allS3 <- do.call(list, lapply(1:length(allfilesS3),FUN= function(x) readRDS(allfilesS3[x])))


#-----------------------------------------------------------------------#
# Warnings
#-----------------------------------------------------------------------#


# Check warning functions
#------------------------------------------------------------------------

check.warning <- function(allscen, war.type, model){
  # allscen  = list of all scenarios loaded above
  # war.type = type of warning
  # model    = model of interst
  allwarnings  <- lapply(1:length(allscen), FUN = function(x) unlist(allscen[[x]][[war.type]][[model]], recursive = F))
  warningtable <- table(unlist(lapply(allwarnings, function(x) x == "NULL")))
  warn.text    <- unlist(allwarnings)[unlist(lapply(allwarnings, function(x) x != "NULL"))]
  
  
  return(list(warningtable=warningtable,warn.text=warn.text))
}

# S1, Warnings
#------------------------------------------------------------------------

# Firth-corrected logistic regression models in full data
S1warningfull    <- check.warning(allscen = allS1, war.type = "Warning", model = "Full")
S1warningsel0157 <- check.warning(allscen = allS1, war.type = "Warning", model = "Selectedp0.157")
S1warningsel0320 <- check.warning(allscen = allS1, war.type = "Warning", model = "Selectedp0.320")

# Intercept re-estimation models in full data
S1warningfull.int    <- check.warning(allscen = allS1, war.type = "Warning.int", model = "Full")
S1warningsel0157.int <- check.warning(allscen = allS1, war.type = "Warning.int", model = "Selectedp0.157")
S1warningsel0320.int <- check.warning(allscen = allS1, war.type = "Warning.int", model = "Selectedp0.320")

# Firth-corrected logistic regression models in bootstrap samples
S1bswarningfull    <- check.warning(allscen = allS1, war.type = "BSWarning", model = "Full")
S1bswarningsel0157 <- check.warning(allscen = allS1, war.type = "BSWarning", model = "Selectedp0.157")
S1bswarningsel0320 <- check.warning(allscen = allS1, war.type = "BSWarning", model = "Selectedp0.320")


# Intercept re-estimation models in bootstrap samples
S1bswarningfull.int    <- check.warning(allscen = allS1, war.type = "BSWarning.int", model = "Full")
S1bswarningsel0157.int <- check.warning(allscen = allS1, war.type = "BSWarning.int", model = "Selectedp0.157")
S1bswarningsel0320.int <- check.warning(allscen = allS1, war.type = "BSWarning.int", model = "Selectedp0.320")


# S2, Warnings
#------------------------------------------------------------------------

# Firth-corrected logistic regression models in full data
S2warningfull    <- check.warning(allscen = allS2, war.type = "Warning", model = "Full")
S2warningsel0157 <- check.warning(allscen = allS2, war.type = "Warning", model = "Selectedp0.157")
S2warningsel0320 <- check.warning(allscen = allS2, war.type = "Warning", model = "Selectedp0.320")

# Intercept re-estimation models in full data
S2warningfull.int    <- check.warning(allscen = allS2, war.type = "Warning.int", model = "Full")
S2warningsel0157.int <- check.warning(allscen = allS2, war.type = "Warning.int", model = "Selectedp0.157")
S2warningsel0320.int <- check.warning(allscen = allS2, war.type = "Warning.int", model = "Selectedp0.320")

# Firth-corrected logistic regression models in bootstrap samples
S2bswarningfull    <- check.warning(allscen = allS2, war.type = "BSWarning", model = "Full")
S2bswarningsel0157 <- check.warning(allscen = allS2, war.type = "BSWarning", model = "Selectedp0.157")
S2bswarningsel0320 <- check.warning(allscen = allS2, war.type = "BSWarning", model = "Selectedp0.320")

# Intercept re-estimation models in bootstrap samples
S2bswarningfull.int    <- check.warning(allscen = allS2, war.type = "BSWarning.int", model = "Full")
S2bswarningsel0157.int <- check.warning(allscen = allS2, war.type = "BSWarning.int", model = "Selectedp0.157")
S2bswarningsel0320.int <- check.warning(allscen = allS2, war.type = "BSWarning.int", model = "Selectedp0.320")

# S3, Warnings
#------------------------------------------------------------------------

# Firth-corrected logistic regression models in full data
S3warningfull    <- check.warning(allscen = allS3, war.type = "Warning", model = "Full")
S3warningsel0157 <- check.warning(allscen = allS3, war.type = "Warning", model = "Selectedp0.157")
S3warningsel0320 <- check.warning(allscen = allS3, war.type = "Warning", model = "Selectedp0.320")

# Intercept re-estimation models in full data
S3warningfull.int    <- check.warning(allscen = allS3, war.type = "Warning.int", model = "Full")
S3warningsel0157.int <- check.warning(allscen = allS3, war.type = "Warning.int", model = "Selectedp0.157")
S3warningsel0320.int <- check.warning(allscen = allS3, war.type = "Warning.int", model = "Selectedp0.320")

# Firth-corrected logistic regression models in bootstrap samples
S3bswarningfull    <- check.warning(allscen = allS3, war.type = "BSWarning", model = "Full")
S3bswarningsel0157 <- check.warning(allscen = allS3, war.type = "BSWarning", model = "Selectedp0.157")
S3bswarningsel0320 <- check.warning(allscen = allS3, war.type = "BSWarning", model = "Selectedp0.320")

# Intercept re-estimation models in bootstrap samples
S3bswarningfull.int    <- check.warning(allscen = allS3, war.type = "BSWarning.int", model = "Full")
S3bswarningsel0157.int <- check.warning(allscen = allS3, war.type = "BSWarning.int", model = "Selectedp0.157")
S3bswarningsel0320.int <- check.warning(allscen = allS3, war.type = "BSWarning.int", model = "Selectedp0.320")



