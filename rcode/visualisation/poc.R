
# Load libraries and data  ----
#------------------------------------------------------------------------------#
library(loopR)
library(ggplot2)
Comparison <- readRDS("./data/results_poc/CICovsel_POC.rds")
test <- readRDS("./data/results_poc/CICovsel_POC.rds")

# Create figure proof of concept ----
#------------------------------------------------------------------------------#
Comparison$eventrate <- Comparison$Yint
Comparison$eventrate[Comparison$Yint == 0] <- 0.5
Comparison$eventrate[Comparison$Yint == -1.65] <- 0.2


POCplot.dat <- data.frame(Comparison[,c('nobs','eventrate','bLA', 'bLY', 'Bias2_omit', 'Var')])
POCplot <- nested_loop_plot(POCplot.dat,
                 x = "bLA",
                 grid_rows = "bLY",
                 grid_cols = "eventrate",
                 steps = "nobs",
                 steps_y_height = 0.01,
                 steps_y_base = 0,
                 steps_values_annotate = T,
                 steps_names = c("Sample size"),
                 colors = c("red3","blue"),
                 ylim = c(-0.04,0.04),
                 legend_labels = c("Bias_omit^2","Var_full - Var_omit"),
                 y_name = NULL,
                 x_name = "Association covariate - exposure")



pdf("./results/figures/POC.pdf", width = 20, height = 10)
POCplot
dev.off()
