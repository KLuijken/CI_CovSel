#------------------------------------------------------------------------------#
# Causal inference covariate selection
# K Luijken
#
# Standalone script to generate figures in supplement 2
#------------------------------------------------------------------------------#

# Load librairies ----
#------------------------------------------------------------------------------#

# Helpers ----
#------------------------------------------------------------------------------#
# MSE expressions
omit <- function(n){
  ((vL/vA)*alpha*gamma)^2 + (vY - vA*beta^2 - 2*vL*alpha*beta*gamma)/vA * (1/(n-2))
}
full <- function(n){
  (vY - vA*beta^2 - vL*gamma^2 - 2*vL*alpha*beta*gamma)/(vA-vL*alpha^2) * (1/(n-3))
}

# Create plots ----
#------------------------------------------------------------------------------#

# Main text ----

gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 8

pdf("./results/figures/omit_vs_full.pdf")
par(mar = c(5,5,1,1))
curve(omit,
      xlim=c(10,200),
      ylab='Mean squared error',
      ylim=c(0,.5),
      xlab='Sample size',
      lwd=2,
      col = "blue",
      cex.lab = 2,
      cex.axis = 1.5)
curve(full,add=TRUE,col="red3",lwd=2)
legend("topright",
       c("Omit","Full"),
       col= c("blue","red3"),
       lwd=2,
       cex = 1.5)
dev.off()

# Supplementary file 2 ----

pdf("./results/figures/supp2_effects.pdf")
par(mfrow = c(3,3),
    oma=c(0,0,3,4))
gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab='MSE',ylim=c(0,.5),xlab='Sample size',lwd=2, col = "blue", main = "Reference")
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .2
alpha = .3
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Lower ", alpha)))
curve(full,add=TRUE,col="red3",lwd=2)
#mtext(expression(paste(bold("Vary parameters "),alpha,", ", beta, bold(" and "), gamma)), 
#      cex = 1.5,
#      line=4)

gamma = .2
alpha = .5
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Higher ", alpha)))
curve(full,add=TRUE,col="red3",lwd=2)
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c("Omit","Full"),col= c("blue","red3"),lwd=2)

plot(NULL)

gamma = .2
alpha = .4
beta = .2
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="MSE",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Lower ", beta)))
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .2
alpha = .4
beta = .4
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Higher ", beta)))
curve(full,add=TRUE,col="red3",lwd=2)

plot(NULL)

gamma = .1
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="MSE",ylim=c(0,.5),xlab="Sample size",lwd=2, col = "blue", main = expression(paste("Lower ", gamma)))
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .3
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="Sample size",lwd=2, col = "blue",main = expression(paste("Higher ", gamma)))
curve(full,add=TRUE,col="red3",lwd=2)

dev.off()

pdf("./results/figures/supp2_vars.pdf")
par(mfrow = c(3,3),
    oma=c(0,0,3,4))
gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab='MSE',ylim=c(0,.5),xlab='Sample size',lwd=2, col = "blue", main = "Reference")
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 1.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Lower ", sigma[A]^2)))
curve(full,add=TRUE,col="red3",lwd=2)
#mtext(expression(paste(bold("Vary parameters "),sigma[A]^2,", ", sigma[L]^2, bold(" and "), sigma[Y]^2)), 
#      cex = 1.5,
#      line=4)

gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 3.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Higher ", sigma[A]^2)))
curve(full,add=TRUE,col="red3",lwd=2)
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c("Omit","Full"),col= c("blue","red3"),lwd=2)

plot(NULL)

gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 6
curve(omit,xlim=c(10,200),ylab="MSE",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Lower ", sigma[L]^2)))
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .2
alpha = .4
beta = .3
vY <- 10
vA <- 2.5
vL <- 10
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="",lwd=2, col = "blue", main = expression(paste("Higher ", sigma[L]^2)))
curve(full,add=TRUE,col="red3",lwd=2)

plot(NULL)

gamma = .2
alpha = .4
beta = .3
vY <- 8
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="MSE",ylim=c(0,.5),xlab="Sample size",lwd=2, col = "blue", main = expression(paste("Lower ", sigma[Y]^2)))
curve(full,add=TRUE,col="red3",lwd=2)

gamma = .2
alpha = .4
beta = .3
vY <- 12
vA <- 2.5
vL <- 8
curve(omit,xlim=c(10,200),ylab="",ylim=c(0,.5),xlab="Sample size",lwd=2, col = "blue",main = expression(paste("Higher ", sigma[Y]^2)))
curve(full,add=TRUE,col="red3",lwd=2)

dev.off()
