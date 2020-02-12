####################################################
# Covariate selection project
#
# December 2019 
# Author: K Luijken
# File content: simulation plots
####################################################

require(loopR)

#-----------------------------------------------------------------------#
# Scenario 1
#-----------------------------------------------------------------------#

# Load input
setup <- readRDS("Setup.rds")
for(i in 1:length(setup)){assign(names(setup)[i], unlist(setup[[i]]))}
setwd(paste0(filepath,"/Simoutput/S1"))
allS1 <- do.call(list, lapply(1:nrow(S1), FUN = function(x) readRDS(paste0(x,".rds"))))


#  S1, bias-------------------------------------------------------------#
bias <- matrix(NA, nrow=nrow(S1), ncol=length(Allmarginals), dimnames= list(NULL, Allmarginals))

# Extract information about bias
for(i in 1:nrow(S1)){
  bias[i,"Truth"]          <- mean(allS1[[i]]$Marginal["log(error)",,"Truth"])
  bias[i,"Unadjusted"]     <- mean(allS1[[i]]$Marginal["log(error)",,"Unadjusted"])
  bias[i,"Full"]           <- mean(allS1[[i]]$Marginal["log(error)",,"Full"])
  bias[i,"Selectedp0.157"] <- mean(allS1[[i]]$Marginal["log(error)",,"Selectedp0.157"])
  bias[i,"Selectedp0.320"] <- mean(allS1[[i]]$Marginal["log(error)",,"Selectedp0.320"])
}

S1bias <- data.frame(cbind(S1,bias))
S1bias$bL2AnL <- S1bias$bL2A*sqrt(S1bias$nL/2)

# Create nested loop plot
S1plotbias <- S1bias[,c("nobs", "nL","rhoL","bL2AnL","bAY",Allmodels)]
S1plotbias <- as.data.frame(S1plotbias)
PS1bias    <- nested_loop_plot(resdf = S1plotbias,
                            grid_cols = "nL",
                            grid_rows="nobs",
                            steps = c("rhoL","bAY"),
                            x = "bL2AnL",spu_x_shift = 0.3,
                            legend_name = "Model",
                            legend_labels = Allmodels,
                            ylim=c(-0.2,0.8),
                            steps_y_height = 0.1,
                            steps_y_base = 0.2,
                            y_name = "Bias log(MRR)",
                            x_name = "RhoL",
                            steps_names_annotate = T,
                            hline_intercept = 0,
                            colors = c("Red","Blue", "Green"))

#  S1, relative bias----------------------------------------------------#

# Compute relative bias by dividing by unadjusted bias
relbias <- matrix(NA, nrow=nrow(S1), ncol=3, dimnames = list(NULL,c("Full", "Selectedp0.157","Selectedp0.320")))
relbias[,"Full"] <- bias[,"Full"]/bias[,"Unadjusted"]
relbias[,"Selectedp0.157"] <- bias[,"Selectedp0.157"]/bias[,"Unadjusted"]
relbias[,"Selectedp0.320"] <- bias[,"Selectedp0.320"]/bias[,"Unadjusted"]

S1relbias <- data.frame(cbind(S1,relbias))
S1relbias$bL2AnL <- S1relbias$bL2A*sqrt(S1relbias$nL/2)

# Create nested loop plot
S1plotrelbias <- S1relbias[,c("nobs", "nL","rhoL","bL2AnL","bAY",Allmodels)]
S1plotrelbias <- as.data.frame(S1plotrelbias)
PS1relbias <- nested_loop_plot(resdf = S1plotrelbias,
                            grid_cols = "nL",
                            grid_rows="nobs",
                            steps = c("rhoL","bAY"),
                            x = "bL2AnL",spu_x_shift = 0.3,
                            legend_name = "Model",
                            legend_labels = c("Full","Selectedp0.157","Selectedp0.320"),
                            ylim=c(-0.2,1.8),
                            steps_y_height = 0.1,
                            steps_y_base = 0.2,
                            y_name = "relbias log(MRR)",
                            x_name = "RhoL",
                            steps_names_annotate = T,
                            hline_intercept = 0,
                            colors = c("Red","Blue","Green"))



# S1, MSE---------------------------------------------------------------#

mse <- matrix(NA, nrow=nrow(S1), ncol=length(Allmarginals), dimnames=list(NULL,Allmarginals))

# Extract information about mse
for(i in 1:nrow(S1)){
  mse[i,"Truth"]          <- mean(allS1[[i]]$Marginal["sq.error",,"Truth"])
  mse[i,"Unadjusted"]     <- mean(allS1[[i]]$Marginal["sq.error",,"Unadjusted"])
  mse[i,"Full"]           <- mean(allS1[[i]]$Marginal["sq.error",,"Full"])
  mse[i,"Selectedp0.157"] <- mean(allS1[[i]]$Marginal["sq.error",,"Selectedp0.157"])
  mse[i,"Selectedp0.320"] <- mean(allS1[[i]]$Marginal["sq.error",,"Selectedp0.320"])
}

S1mse <- data.frame(cbind(S1,mse))
S1mse$bL2AnL <- S1mse$bL2A*sqrt(S1mse$nL/2)

# Create nested loop plot
S1plotmse <- S1mse[,c("nobs", "nL","rhoL","bL2AnL","bAY", Allmodels)]
S1plotmse <- as.data.frame(S1plotmse)
PS1mse <- nested_loop_plot(resdf = S1plotmse,
                           grid_cols = "nL",
                           grid_rows="nobs",
                           steps = c("rhoL","bAY"),
                           x = "bL2AnL",spu_x_shift = 0.3,
                           legend_name = "Model",
                           legend_labels = Allmodels,
                           ylim=c(-0.01,0.8),
                           steps_y_height = 0.01,
                           steps_y_base = 0.05,
                           y_name = "MSE log(MRR)",
                           x_name = "RhoL1L2",
                           steps_names_annotate = T,
                           colors = c("Red","Blue","Green"))

# S1, relative MSE-------------------------------------------------------#

# Compute relative mse
relmse <- matrix(NA, nrow=nrow(S1), ncol=3, dimnames=list(NULL,Allmodels))
relmse[,"Full"] <- mse[,"Full"]/mse[,"Unadjusted"]
relmse[,"Selectedp0.157"] <- mse[,"Selectedp0.157"]/mse[,"Unadjusted"]
relmse[,"Selectedp0.320"] <- mse[,"Selectedp0.320"]/mse[,"Unadjusted"]

S1relmse <- data.frame(cbind(S1,relmse))
S1relmse$bL2AnL <- S1relmse$bL2A*sqrt(S1relmse$nL/2)

# Create nested loop plot
S1plotrelmse <- S1relmse[,c("nobs", "nL","rhoL","bL2AnL","bAY", Allmodels)]
S1plotrelmse <- as.data.frame(S1plotrelmse)
PS1relmse <- nested_loop_plot(resdf = S1plotmse,
                           grid_cols = "nL",
                           grid_rows="nobs",
                           steps = c("rhoL","bAY"),
                           x = "bL2AnL",spu_x_shift = 0.3,
                           legend_name = "Model",
                           #legend_labels = Allmodels,
                           ylim=c(-0.01,2),
                           steps_y_height = 0.01,
                           steps_y_base = 0.05,
                           y_name = "MSE log(MRR)",
                           x_name = "RhoL1L2",
                           steps_names_annotate = T,
                           colors = c("Red","Blue","Green"))

# S1, MSE check---------------------------------------------------------#

# Full
mseN <- S1mse[,"Full"] * S1mse[,"nobs"]
S1mseN <- cbind(S1,mseN)
S1mseN <- as.data.frame(S1mseN)

S1mseN$bL2AnL <- S1mseN$bL2A*sqrt(S1mseN$nL/2)
S1mseN <- S1mseN[,c("nL","bL2AnL","rhoL","bAY","nobs","mseN")]

PS1mseN <- nested_loop_plot(resdf = S1mseN,
                              grid_cols = "nL",
                              grid_rows="rhoL",
                              steps = c("bAY","bL2AnL"),
                              x = "nobs",spu_x_shift = 0.3,
                              #legend_name = "Model",
                              #legend_labels = Allmodels,
                              ylim=c(0,15),
                              steps_y_height = 0.01,
                              steps_y_base = 0.05,
                              y_name = "MSE log(MRR)",
                              x_name = "RhoL1L2",
                              steps_names_annotate = T,
                              colors = "Blue")


# S1, Coverage----------------------------------------------------------#
coverage <- matrix(NA, nrow=nrow(S1), ncol=length(Allmarginals), dimnames=list(NULL,Allmarginals))
# Extract information about coverage
for(i in 1:nrow(S1)){
  coverage[i,"Truth"]          <- mean(allS1[[i]]$Marginal["coverage",,"Truth"])
  coverage[i,"Unadjusted"]     <- mean(allS1[[i]]$Marginal["coverage",,"Unadjusted"])
  coverage[i,"Full"]           <- mean(allS1[[i]]$Marginal["coverage",,"Full"])
  coverage[i,"Selectedp0.157"] <- mean(allS1[[i]]$Marginal["coverage",,"Selectedp0.157"])
  coverage[i,"Selectedp0.320"] <- mean(allS1[[i]]$Marginal["coverage",,"Selectedp0.320"])
}

S1coverage <- data.frame(cbind(S1,coverage))
S1coverage$bL2AnL <- S1coverage$bL2A*sqrt(S1coverage$nL/2)

# Create nested loop plot
S1plotcoverage <- S1coverage[,c("nobs", "nL","rhoL","bL2AnL","bAY", Allmodels)]
S1plotcoverage <- as.data.frame(S1plotcoverage)
PS1coverage <- nested_loop_plot(resdf = S1plotcoverage,
                           grid_cols = "nL",
                           grid_rows="nobs",
                           steps = c("rhoL","bAY"),
                           x = "bL2AnL",spu_x_shift = 0.3,
                           legend_name = "Model",
                           #legend_labels = Allmodels,
                           ylim=c(-0.01,1),
                           steps_y_height = 0.01,
                           steps_y_base = 0.05,
                           y_name = "coverage log(MRR)",
                           x_name = "RhoL1L2",
                           steps_names_annotate = T,
                           colors = c("Red","Blue","Green"))



# S1, RMSD--------------------------------------------------------------#

rmsd <- matrix(NA, nrow=nrow(S1), ncol=length(Allmodels), dimnames = list(NULL, Allmodels))

# Extract rmsd
for(i in 1:nrow(S1)){
  rmsd[i,"Full"]           <- mean(allS1[[i]]$Model["RMSD bAY",,"Full"])
  rmsd[i,"Selectedp0.157"] <- mean(allS1[[i]]$Model["RMSD bAY",,"Selectedp0.157"])
  rmsd[i,"Selectedp0.320"] <- mean(allS1[[i]]$Model["RMSD bAY",,"Selectedp0.320"])
}

S1rmsd <- data.frame(cbind(S1,rmsd))
S1rmsd$bL2AnL <- S1rmsd$bL2A*sqrt(S1rmsd$nL/2)

# Create nested loop plot
S1plotrmsd <- S1rmsd[,c("nobs", "nL","rhoL","bL2AnL","bAY", Allmodels)]
S1plotrmsd <- as.data.frame(S1plotrmsd)
PS1rmsd <- nested_loop_plot(resdf = S1plotrmsd,
                           grid_cols = "nL",
                           grid_rows="nobs",
                           steps = c("rhoL","bAY"),
                           x = "bL2AnL",spu_x_shift = 0.3,
                           legend_name = "Model",
                           #legend_labels = Allmodels,
                           ylim=c(-0.01,2),
                           steps_y_height = 0.01,
                           steps_y_base = 0.05,
                           y_name = "rmsd bAY",
                           x_name = "RhoL1L2",
                           steps_names_annotate = T,
                           colors = c("Red","Blue","Green"))



#-----------------------------------------------------------------------#
# Scenario 2
#-----------------------------------------------------------------------#
allfilesS2 <- list.files(paste0(filepath,"/Simoutput/S2"), pattern = ".rds")
setwd(paste0(filepath,"/Simoutput/S2"))
allS2 <- do.call(list, lapply(1:length(allfilesS2),FUN= function(x) readRDS(allfilesS2[x])))

bias <- matrix(NA, nrow=nrow(S2), ncol=length(Allmarginals))
colnames(bias) <- Allmarginals

for(i in 1:length(allS2)){
  bias[i,"Truth"]          <- mean(allS2[[i]]$Marginal["log(error)",,"Truth"])
  bias[i,"Unadjusted"]     <- mean(allS2[[i]]$Marginal["log(error)",,"Unadjusted"])
  bias[i,"Full"]           <- mean(allS2[[i]]$Marginal["log(error)",,"Full"])
  bias[i,"Selectedp0.157"] <- mean(allS2[[i]]$Marginal["log(error)",,"Selectedp0.157"])
  bias[i,"Selectedp0.320"] <- mean(allS2[[i]]$Marginal["log(error)",,"Selectedp0.320"])
}


S2bias <- data.frame(cbind(S2,bias))
S2bias$bL2AnL <- S2bias$bL2A*sqrt(S2bias$nL/2)
S2bias$bL2YnL <- S2bias$bL2Y*sqrt(S2bias$nL/2)

S2plotbias <- S2bias[,c("nobs", "nL","rhoL","bL2AnL","bL2YnL","bAY",Allmodels)]
S2plotbias <- as.data.frame(S2plotbias)
PS2bias <- nested_loop_plot(resdf = S2plotbias,
                            grid_cols = "nL",
                            grid_rows="nobs",
                            steps = c("rhoL","bAY","bL2YnL"),
                            x = "bL2AnL",spu_x_shift = 0.1,
                            legend_name = "Model",
                            legend_labels = Allmodels,
                            ylim=c(-0.2,0.8),
                            steps_y_height = 0.1,
                            steps_y_base = 0.2,
                            y_name = "Bias log(MRR)",
                            x_name = "RhoL",
                            steps_names_annotate = T,
                            hline_intercept = 0,
                            colors = c("Red","Blue","Green"))


# S2, MSE---------------------------------------------------------------#

mse <- matrix(NA, nrow=nrow(S2), ncol=length(Allmarginals), dimnames=list(NULL,Allmarginals))

# Extract information about mse
for(i in 1:nrow(S2)){
  mse[i,"Truth"]          <- mean(allS2[[i]]$Marginal["log(sq.error)",,"Truth"])
  mse[i,"Unadjusted"]     <- mean(allS2[[i]]$Marginal["log(sq.error)",,"Unadjusted"])
  mse[i,"Full"]           <- mean(allS2[[i]]$Marginal["log(sq.error)",,"Full"])
  mse[i,"Selectedp0.157"] <- mean(allS2[[i]]$Marginal["log(sq.error)",,"Selectedp0.157"])
  mse[i,"Selectedp0.320"] <- mean(allS2[[i]]$Marginal["log(sq.error)",,"Selectedp0.320"])
}

S2mse <- data.frame(cbind(S2,mse))
S2mse$bL2AnL <- S2mse$bL2A*sqrt(S2mse$nL/2)
S2mse$bL2YnL <- S2mse$bL2Y*sqrt(S2mse$nL/2)

# Create nested loop plot
S2plotmse <- S2mse[,c("nobs", "nL","rhoL","bL2AnL","bL2YnL","bAY", Allmodels)]
S2plotmse <- as.data.frame(S2plotmse)
PS2mse <- nested_loop_plot(resdf = S2plotmse,
                           grid_cols = "nL",
                           grid_rows="nobs",
                           steps = c("rhoL","bAY","bL2YnL"),
                           x = "bL2AnL",spu_x_shift = 0.1,
                           legend_name = "Model",
                           legend_labels = Allmodels,
                           ylim=c(-0.3,0.5),
                           steps_y_height = 0.06,
                           steps_y_base = 0,
                           y_name = "MSE log(MRR)",
                           x_name = "bA-L2 * nL",
                           steps_names_annotate = T,
                           colors = c("Red","Blue","Green"))

#-----------------------------------------------------------------------#
# Scenario 3
#-----------------------------------------------------------------------#
allfilesS3 <- list.files(paste0(filepath,"/Simoutput/S3"), pattern = ".rds")
setwd(paste0(filepath,"/Simoutput/S3"))
allS3 <- do.call(list, lapply(1:length(allfilesS3),FUN= function(x) readRDS(allfilesS3[x])))

bias <- matrix(NA, nrow=nrow(S3), ncol=length(Allmarginals))
colnames(bias) <- Allmarginals

for(i in 1:length(allS3)){
  bias[i,"Truth"]          <- mean(allS3[[i]]$Marginal["log(error)",,"Truth"])
  bias[i,"Unadjusted"]     <- mean(allS3[[i]]$Marginal["log(error)",,"Unadjusted"])
  bias[i,"Full"]           <- mean(allS3[[i]]$Marginal["log(error)",,"Full"])
  bias[i,"Selectedp0.157"] <- mean(allS3[[i]]$Marginal["log(error)",,"Selectedp0.157"])
  bias[i,"Selectedp0.320"] <- mean(allS3[[i]]$Marginal["log(error)",,"Selectedp0.320"])
}


S3bias <- data.frame(cbind(S3,bias))
S3bias$bL2AnL <- S3bias$bL2A*sqrt(S3bias$nL/2)
S3bias$bL2YnL <- S3bias$bL2Y*sqrt(S3bias$nL/2)

S3plotbias <- S3bias[,c("nobs", "nL","rhoL","bL2AnL","bL2YnL","bAY",Allmodels)]
S3plotbias <- as.data.frame(S3plotbias)
PS3bias <- nested_loop_plot(resdf = S3plotbias,
                            grid_cols = "nL",
                            grid_rows="nobs",
                            steps = c("rhoL","bAY","bL2AnL"),
                            x = "bL2YnL",spu_x_shift = 0.1,
                            legend_name = "Model",
                            #legend_labels = Allmodels,
                            ylim=c(-0.2,0.8),
                            steps_y_height = 0.1,
                            steps_y_base = 0.2,
                            y_name = "Bias log(MRR)",
                            x_name = "RhoL",
                            steps_names_annotate = T,
                            hline_intercept = 0,
                            colors = c("Red","Blue","Green"))
