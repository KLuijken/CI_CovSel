####################################################
# Covariate selection project
#
# January 2020 
# Author: K Luijken
# File content: Diagnostics simulations
####################################################

# Load packages
require(loopR)

# Load simulation output
setup <- readRDS("C:/Users/Eigenaar/Desktop/LUMC/Projects/CovSel/Simulations/Versioncontrol/CI_CovSel/Simoutput/PlotsTables/Setup.rds")
for(i in 1:length(setup)){assign(names(setup)[i], unlist(setup[[i]]))}
setwd("C:/Users/Eigenaar/Desktop/LUMC/Projects/CovSel/Simulations/Versioncontrol/CI_CovSel/Simoutput/S1")
allS1 <- do.call(list, lapply(1:nrow(S1), FUN = function(x) readRDS(paste0(x,".rds"))))


# S1, MSE check---------------------------------------------------------#

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


# Full
mseNfull <- S1mse[,"Full"] * S1mse[,"nobs"]
S1mseNfull <- cbind(S1,mseNfull)
S1mseNfull <- as.data.frame(S1mseNfull)

S1mseNfull$bL2AnL <- S1mseNfull$bL2A*sqrt(S1mseNfull$nL/2)
S1mseNfull <- S1mseNfull[,c("nL","bL2AnL","rhoL","bAY","nobs","mseN")]

PS1mseNfull <- nested_loop_plot(resdf = S1mseNfull,
                            grid_cols = "nL",
                            grid_rows="rhoL",
                            steps = c("bAY","bL2AnL"),
                            x = "nobs",spu_x_shift = 0.3,
                            ylim=c(0,15),
                            steps_y_height = 0.01,
                            steps_y_base = 0.05,
                            y_name = "MSE of MRR * nobs",
                            x_name = "nobs",
                            steps_names_annotate = T,
                            colors = "Blue")

# Selected, cutoff p = 0.157
mseNsel0157 <- S1mse[,"Selectedp0.157"] * S1mse[,"nobs"]
S1mseNsel0157 <- cbind(S1,mseNsel0157)
S1mseNsel0157 <- as.data.frame(S1mseNsel0157)

S1mseNsel0157$bL2AnL <- S1mseNsel0157$bL2A*sqrt(S1mseNsel0157$nL/2)
S1mseNsel0157 <- S1mseNsel0157[,c("nL","bL2AnL","rhoL","bAY","nobs","mseN")]

PS1mseNsel0157 <- nested_loop_plot(resdf = S1mseNsel0157,
                                grid_cols = "nL",
                                grid_rows="rhoL",
                                steps = c("bAY","bL2AnL"),
                                x = "nobs",spu_x_shift = 0.3,
                                ylim=c(0,15),
                                steps_y_height = 0.01,
                                steps_y_base = 0.05,
                                y_name = "MSE of MRR * nobs",
                                x_name = "nobs",
                                steps_names_annotate = T,
                                colors = "Blue")

# Selected, cutoff p = 0.320
mseNsel0320 <- S1mse[,"Selectedp0.320"] * S1mse[,"nobs"]
S1mseNsel0320 <- cbind(S1,mseNsel0320)
S1mseNsel0320 <- as.data.frame(S1mseNsel0320)

S1mseNsel0320$bL2AnL <- S1mseNsel0320$bL2A*sqrt(S1mseNsel0320$nL/2)
S1mseNsel0320 <- S1mseNsel0320[,c("nL","bL2AnL","rhoL","bAY","nobs","mseN")]

PS1mseNsel0320 <- nested_loop_plot(resdf = S1mseNsel0320,
                                   grid_cols = "nL",
                                   grid_rows="rhoL",
                                   steps = c("bAY","bL2AnL"),
                                   x = "nobs",spu_x_shift = 0.3,
                                   ylim=c(0,15),
                                   steps_y_height = 0.01,
                                   steps_y_base = 0.05,
                                   y_name = "MSE of MRR * nobs",
                                   x_name = "nobs",
                                   steps_names_annotate = T,
                                   colors = "Blue")


pfd("Simoutput_Diagnostics_S1.pdf", width = 25, height = 15)
print(PS1mseNfull)
print(PS1mseNsel0157)
print(PS1mseNsel0320)
dev.off()