# Load libraries and data  ----
#------------------------------------------------------------------------------#

library(looplot)
library(ggplot2)
Comparison <- readRDS("./data/poc/CICovsel_POC.rds")

# Overwrite function
label_both_custom <- function(labels, multi_line = TRUE, sep = ": ") {
  value <- ggplot2::label_value(labels, multi_line = multi_line)
  variable <- ggplot2:::label_variable(labels, multi_line = multi_line)
  variable <- Map(gsub, "_labels_", "", variable) # customization
  print(variable)
  if (multi_line) {
    out <- vector("list", length(value))
    for (i in seq_along(out)) {
      out[[i]] <- parse(text = paste(variable[[i]], value[[i]], sep = sep))
    }
    
  }
  else {
    value <- do.call("paste", list(value, sep = ", "))
    variable <- do.call("paste", list(variable, sep = ", "))
    out <- Map(paste, variable, value, sep = sep)
  }
  
  list(unname(unlist(out)))
}

assignInNamespace("label_both_custom", label_both_custom, ns="looplot")



# Create figure proof of concept for MRR ----
#------------------------------------------------------------------------------#
Comparison$Eventrate <- Comparison$Yint
Comparison$Eventrate[Comparison$Yint == 0] <- 0.5
Comparison$Eventrate[Comparison$Yint == -1.65] <- 0.2
Comparison <- Comparison[Comparison$nobs < 300,]

POCplot.dat <- data.frame(Comparison[,c('nobs','Eventrate','bLA', 'bLY', 'Bias2_omit_MRR', 'Var_MRR')])
names(POCplot.dat)[names(POCplot.dat) == "Eventrate"] <- "Event_fraction"

POCplot <- nested_loop_plot(POCplot.dat %>% rename(`beta[LY]` = bLY),
                            x = "bLA",
                            grid_rows = "beta[LY]", 
                            grid_cols = quote("Event_fraction"), 
                            steps = "nobs",
                            steps_y_height = 0.01,
                            steps_y_base = -0.01,
                            steps_values_annotate = T,
                            steps_names = c("Sample size"),
                            colors = c("black","grey"),
                            ylim = c(-0.03,0.04),
                            legend_labels = c(expression("Bias"["omit"]^2),expression("Var"["full"] - "Var"["omit"])),
                            y_name = "Mean squared error of log(mRR)",
                            x_name = expression("Association covariate - exposure" (beta["LA"])),
                            post_processing = list(
                              add_custom_theme = list(
                                legend.title=element_blank(),
                                legend.position="bottom",
                                axis.text = element_text(size = 14),
                                axis.title = element_text(size = 20),
                                strip.text = element_text(size = 16),
                                legend.text = element_text(size = 16)
                                )),
                            spu_x_shift = 0.3)


pdf("./results/figures/POC_MRR.pdf", width = 12, height = 10)
POCplot
dev.off()


# Create figure proof of concept for OR ----
#------------------------------------------------------------------------------#
Comparison$Eventrate <- Comparison$Yint
Comparison$Eventrate[Comparison$Yint == 0] <- 0.5
Comparison$Eventrate[Comparison$Yint == -1.65] <- 0.2
Comparison <- Comparison[Comparison$nobs < 300,]

POCplot.dat <- data.frame(Comparison[,c('nobs','Eventrate','bLA', 'bLY', 'Bias2_omit_OR', 'Var_OR')])
POCplot <- nested_loop_plot(POCplot.dat,
                            x = "bLA",
                            grid_rows = "bLY",
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
                            x_name = "Association covariate - exposure (bLA)")



pdf("./results/figures/POC_OR.pdf", width = 14, height = 10)
POCplot
dev.off()
