# Load libraries and data  ----
#------------------------------------------------------------------------------#
library(loopR)
library(ggplot2)
Comparison <- readRDS("./data/poc/CICovsel_POC.rds")

# Create figure proof of concept for MRR ----
#------------------------------------------------------------------------------#
Comparison$Eventrate <- Comparison$Yint
Comparison$Eventrate[Comparison$Yint == 0] <- 0.5
Comparison$Eventrate[Comparison$Yint == -1.65] <- 0.2
Comparison <- Comparison[Comparison$nobs < 300,]
names(Comparison)[names(Comparison) == "bLA"] <- "bAL"
names(Comparison)[names(Comparison) == "bLY"] <- "bYL"


POCplot.dat <- data.frame(Comparison[,c('nobs','Eventrate','bAL', 'bYL', 'Bias2_omit_MRR', 'Var_MRR')])
names(POCplot.dat)[names(POCplot.dat) == "Eventrate"] <- "Event fraction"

POCplot <- nested_loop_plot(POCplot.dat,
                            x = "bAL",
                            grid_rows = "bYL",
                            grid_cols = quote("Event fraction"),
                            steps = "nobs",
                            steps_y_height = 0.01,
                            steps_y_base = -0.01,
                            steps_values_annotate = T,
                            steps_names = c("Sample size"),
                            colors = c("red3","blue"),
                            ylim = c(-0.03,0.04),
                            legend_labels = c(expression("Bias"["omit"]^2),expression("Var"["full"] - "Var"["omit"])),
                            y_name = "Mean squared error of log(mRR)",
                            x_name = "Association exposure - covariate (bAL)",
                            post_processing = list(
                              add_custom_theme = list(
                                legend.title=element_blank(),
                                legend.position="bottom",
                                axis.text = element_text(size = 14),
                                axis.title = element_text(size = 20),
                                strip.text = element_text(size = 16),
                                legend.text = element_text(size = 16)
                              )
                            ),
                            spu_x_shift = 0.3)


pdf("./results/figures/POC.pdf", width = 10, height = 10)
POCplot
dev.off()


# Create figure proof of concept for OR ----
#------------------------------------------------------------------------------#
Comparison$Eventrate <- Comparison$Yint
Comparison$Eventrate[Comparison$Yint == 0] <- 0.5
Comparison$Eventrate[Comparison$Yint == -1.65] <- 0.2
Comparison <- Comparison[Comparison$nobs < 300,]
names(Comparison)[names(Comparison) == "bLA"] <- "bAL"
names(Comparison)[names(Comparison) == "bLY"] <- "bYL"

POCplot.dat <- data.frame(Comparison[,c('nobs','Eventrate','bAL', 'bYL', 'Bias2_omit_OR', 'Var_OR')])
POCplot <- nested_loop_plot(POCplot.dat,
                            x = "bAL",
                            grid_rows = "bYL",
                            grid_cols = "Eventrate",
                            steps = "nobs",
                            steps_y_height = 0.01,
                            steps_y_base = -0.02,
                            steps_values_annotate = T,
                            steps_names = c("Sample size"),
                            colors = c("red3","blue"),
                            ylim = c(-0.04,0.09),
                            legend_labels = c(expression("Bias"["omit"]^2),expression("Var"["full"] - "Var"["omit"])),
                            y_name = "Mean squared error of conditional OR",
                            x_name = "Association exposure - covariate (bAL)")



pdf("./results/figures/POC_OR.pdf", width = 14, height = 10)
POCplot
dev.off()
