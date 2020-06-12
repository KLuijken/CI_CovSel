#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Standalone script to generate figures in supplement 1
#------------------------------------------------------------------------------#

# Load librairies ----
#------------------------------------------------------------------------------#
library(ggplot2)
library(cowplot)

# Small simulation to verify expressions ----
#------------------------------------------------------------------------------#
# Set some randomly chosen values 
sigmaA <- 1.2
sigmaL <- 1.2
sigmaY <- 1.2
alpha  <- 0.3
beta   <- 0.5
gamma  <- 0.6

# Expressions
# Full model
bias_full <- beta - beta
var_full  <- sigmaY^2/((n-1)*sigmaA^2)
mse_full  <- sigmaY^2/((n-1)*sigmaA^2)

# Reduced model
bias_red  <- (sigmaL^2*alpha*gamma)/(sigmaA^2+sigmaL^2*alpha^2)
var_red   <- (1/(n-1)) * ((sigmaA^2*sigmaY^2 + 
                             sigmaL^2*sigmaY^2*alpha^2 + 
                             sigmaA^2*sigmaL^2*gamma^2)/
                            ((sigmaA^2 + sigmaL^2*alpha^2)^2))
mse_red   <- (((n-1)*sigmaL^4*alpha^2*gamma^2) + 
                sigmaA^2*sigmaY^2 + 
                sigmaL^2*sigmaY^2*alpha^2 + 
                sigmaA^2*sigmaL^2*gamma^2)/
                             ((n - 1)*((sigmaA^2 + sigmaL^2*alpha^2)^2))

# Large sample properties ----
n      <- 1e7 

# Generate data according to model in supplement figure
L <- rnorm(n, mean = 0, sd = sigmaL)
A <- rnorm(n, mean = alpha*L, sd = sigmaA)
Y <- rnorm(n, mean = gamma*L+beta*A, sd = sigmaY)

# Fit full model
m_full      <- lm(Y~A+L)
m_bias_full <- coef(summary(m_full))["A","Estimate"] - beta
m_var_full  <- (coef(summary(m_full))["A","Std. Error"])^2

# Check values in console
c(m_bias_full, bias_full)
c(m_var_full, var_full) # not so similar?

# Fit reduced model
m_red      <- lm(Y~A)
m_bias_red <- coef(summary(m_red))["A","Estimate"] - beta
m_var_red  <- (coef(summary(m_red))["A","Std. Error"])^2

# Check values in console
c(m_bias_red, bias_red)
c(m_var_red, var_red) # not so similar?

# Finite samples ----
n      <- 1e2
rep    <- 1000000
m_coef <- 
  m_bias <- 
  m_var  <- matrix(NA, nrow=rep,ncol=2, dimnames = list(NULL,c("full","reduced")))

for(i in 1:rep){
  # Generate data according to model in supplement figure
  L <- rnorm(n, mean = 0, sd = sigmaL)
  A <- rnorm(n, mean = alpha*L, sd = sigmaA)
  Y <- rnorm(n, mean = gamma*L+beta*A, sd = sigmaY)
  
  # Fit full model
  m_full           <- lm(Y~A+L)
  m_coef[i,"full"] <- coef(summary(m_full))["A","Estimate"]
  m_bias[i,"full"] <- coef(summary(m_full))["A","Estimate"] - beta
  m_var[i,"full"]  <- (coef(summary(m_full))["A","Std. Error"])^2
  
  # Fit reduced model
  m_red               <- lm(Y~A)
  m_coef[i,"reduced"] <- coef(summary(m_red))["A","Estimate"]
  m_bias[i,"reduced"] <- coef(summary(m_red))["A","Estimate"] - beta
  m_var[i,"reduced"]  <- (coef(summary(m_red))["A","Std. Error"])^2

}

m_mse_full <- mean((m_coef[,"full"]-beta)^2)
m_mse_red <- mean((m_coef[,"reduced"]-beta)^2)

# Check values in console: mse
c(m_mse_full, mse_full)
c(m_mse_red, mse_red)

# Check values in console: bias
c(mean(m_bias[,"full"]), bias_full)
c(mean(m_bias[,"reduced"]), bias_red)

# Check values in console: variance
c(mean(m_var[,"full"]), var_full)
c(mean(m_var[,"reduced"]), var_red)

# Generate figures ----
#------------------------------------------------------------------------------#

# Helpers ----
provide_expressions <- function(sigmaA, sigmaL, sigmaY, alpha, beta, gamma, n){
  
  # Full model
  bias_full <- beta - beta
  var_full  <- sigmaY^2/((n-1)*sigmaA^2)
  mse_full  <- sigmaY^2/((n-1)*sigmaA^2)
  
  # Reduced model
  bias_red  <- (sigmaL^2*alpha*gamma)/(sigmaA^2+sigmaL^2*alpha^2)
  var_red   <- (1/(n-1)) * ((sigmaA^2*sigmaY^2 + 
                               sigmaL^2*sigmaY^2*alpha^2 + 
                               sigmaA^2*sigmaL^2*gamma^2)/
                              ((sigmaA^2 + sigmaL^2*alpha^2)^2))
  mse_red   <- (((n-1)*sigmaL^4*alpha^2*gamma^2) + 
                  sigmaA^2*sigmaY^2 + 
                  sigmaL^2*sigmaY^2*alpha^2 + 
                  sigmaA^2*sigmaL^2*gamma^2)/
    ((n - 1)*((sigmaA^2 + sigmaL^2*alpha^2)^2))
  
  return(list(bias_full=bias_full, var_full=var_full, mse_full=mse_full,
              bias_red=bias_red, var_red=var_red, mse_red=mse_red))
}

generate_plot_input <- function(sigmaA, sigmaL, sigmaY, alpha, beta, gamma, n,
                                varied_parameter,
                                varied_par_name){
  params <- expand.grid(sigmaA = sigmaA,
                        sigmaL = sigmaL,
                        sigmaY = sigmaY,
                        alpha = alpha,
                        beta = beta,
                        gamma = gamma,
                        n = n)
  
  bias <- 
    var <- 
    mse <- cbind(varied_parameter,data.frame(matrix(NA, nrow = length(varied_parameter), ncol = 2, dimnames = list(NULL, c("full","reduced")))))
  
  for(i in 1:length(varied_parameter)){
    plotinput <- provide_expressions(sigmaA = params[i,'sigmaA'],
                                     sigmaL = params[i,'sigmaL'],
                                     sigmaY = params[i,'sigmaY'],
                                     alpha = params[i,'alpha'],
                                     beta = params[i,'beta'],
                                     gamma = params[i,'gamma'],
                                     n = params[i,'n'])
    bias[i,c("full", "reduced")] <- c(plotinput$bias_full,plotinput$bias_red)
    var[i,c("full", "reduced")]  <- c(plotinput$var_full,plotinput$var_red)
    mse[i,c("full", "reduced")]  <- c(plotinput$mse_full,plotinput$mse_red)
  }
  
  plot_bias <- ggplot(bias) +
    geom_point(aes(x = varied_parameter, y = full, color = "Full")) + 
    geom_line(aes(x = varied_parameter, y = full, color = "Full")) +
    geom_point(aes(x = varied_parameter, y = reduced, color = "Reduced")) + 
    geom_line(aes(x = varied_parameter, y = reduced, color = "Reduced")) +
    scale_color_manual(name = "Model", values = c("Full" = "navyblue","Reduced"= "red")) +
    ylab(expression("Bias("~hat(beta)~")")) +
    xlab(element_text(varied_par_name))+
    ylim(0,0.5)+
    theme_classic() + theme(legend.position ="right")
  
  plot_var <- ggplot(var) +
    geom_point(aes(x = varied_parameter, y = full, color = "full")) + 
    geom_line(aes(x = varied_parameter, y = full, color = "full")) +
    geom_point(aes(x = varied_parameter, y = reduced, color = "reduced")) + 
    geom_line(aes(x = varied_parameter, y = reduced, color = "reduced")) +
    scale_color_manual(name = "Model", values = c("full" = "navyblue","reduced"= "red")) + ylab(expression("Var("~hat(beta)~")")) + 
    xlab(element_text(varied_par_name))+
    ylim(0,0.025)+
    theme_classic()
  
  plot_mse <- ggplot(mse) +
    geom_point(aes(x = varied_parameter, y = full, color = "full")) + 
    geom_line(aes(x = varied_parameter, y = full, color = "full")) +
    geom_point(aes(x = varied_parameter, y = reduced, color = "reduced")) + 
    geom_line(aes(x = varied_parameter, y = reduced, color = "reduced")) +
    scale_color_manual(name = "Model", values = c("full" = "navyblue","reduced"= "red")) +
    ylab(expression("MSE("~hat(beta)~")")) + 
    xlab(element_text(varied_par_name))+
    ylim(0,0.28)+
    theme_classic()
  
  legend <- cowplot::get_legend(plot_bias + theme(legend.position = "right"))
  p_grid <- cowplot::plot_grid(plotlist = 
                                 list(plot_bias + theme(legend.position = "none"), 
                                      plot_var + theme(legend.position = "none"),
                                      plot_mse + theme(legend.position = "none")),
                               ncol = 3)
  plot_combined <- cowplot::plot_grid(p_grid, legend, 
                                      ncol = 2,
                                      rel_widths = c(4, 0.4))
  return(plot_combined)
}



# Sample size ----

# Fix all parameters but n and generate plot
plot_samplesize <- generate_plot_input(sigmaA = 0.5,
                                       sigmaL = 0.5,
                                       sigmaY = 0.5,
                                       alpha = 0.5,
                                       beta = 0.5,
                                       gamma = 0.5,
                                       n = c(50,100,200, 300, 400, 500),
                                       varied_parameter = c(50,100,200, 300, 400, 500),
                                       varied_par_name = "Sample size")
pdf(file= paste0("./results/figures/supp_samplesize.pdf"), width = 12, height = 4)
plot_samplesize
dev.off()


# plot_samplesize <- generate_plot_input(sigmaA = 0.5,
#                                        sigmaL = 0.1,
#                                        sigmaY = 0.1,
#                                        alpha = 0.1,
#                                        beta = 0.1,
#                                        gamma = 0.1,
#                                        n = c(50,100,200, 300, 400, 500),
#                     varied_parameter = c(50,100,200, 300, 400, 500),
#                     varied_par_name = "Sample size")
# pdf(file= paste0("./results/figures/supp_samplesize_variation.pdf"), width = 12, height = 4)
# plot_samplesize
# dev.off()

# (Near)-instrument ----

# Fix all parameters but gamma and generate plot
plot_instrument <- generate_plot_input(sigmaA = 0.5,
                                       sigmaL = 0.5,
                                       sigmaY = 0.5,
                                       alpha = 0.5,
                                       beta = 0.5,
                                       gamma = c(0,0.1, 0.2,0.3,0.4,0.5),
                                       n = 200,
                                       varied_parameter = c(0,0.1, 0.2,0.3,0.4,0.5),
                                       varied_par_name = "Association Y-L")
pdf(file= paste0("./results/figures/supp_instrument.pdf"), width = 12, height = 4)
plot_instrument
dev.off()

# Additionally vary alpha with fixed gamma to study effect of instrument strength
plot_instrument <- generate_plot_input(sigmaA = 0.5,
                                       sigmaL = 0.5,
                                       sigmaY = 0.5,
                                       alpha = c(0.1, 0.2,0.3,0.4,0.5, 0.6),
                                       beta = 0.5,
                                       gamma = 0.1,
                                       n = 200,
                                       varied_parameter = c(0.1, 0.2,0.3,0.4,0.5, 0.6),
                                       varied_par_name = "Association A-L")
pdf(file= paste0("./results/figures/supp_instrument_variation.pdf"), width = 12, height = 4)
plot_instrument
dev.off()

# Predictor ----

# Fix all parameters but alpha and generate plot
plot_predictor <- generate_plot_input(sigmaA = 0.5,
                                       sigmaL = 0.5,
                                       sigmaY = 0.5,
                                       alpha = c(0,0.1, 0.2,0.3,0.4,0.5),
                                       beta = 0.5,
                                       gamma = 0.5,
                                       n = 200,
                                       varied_parameter = c(0,0.1, 0.2,0.3,0.4,0.5),
                                       varied_par_name = "Association A-L")
pdf(file= paste0("./results/figures/supp_predictor.pdf"), width = 12, height = 4)
plot_predictor
dev.off()

# Explained variance Y ----

# Fix all parameters but sigmaY and generate plot
plot_explained <- generate_plot_input(sigmaA = 0.5,
                                      sigmaL = 0.5,
                                      sigmaY = c(0,0.1, 0.2,0.3,0.4,0.5),
                                      alpha = 0.5,
                                      beta = 0.5,
                                      gamma = 0.5,
                                      n = 200,
                                      varied_parameter = c(0,0.1, 0.2,0.3,0.4,0.5),
                                      varied_par_name = "Residual variance Y")
pdf(file= paste0("./results/figures/supp_explained.pdf"), width = 12, height = 4)
plot_explained
dev.off()

