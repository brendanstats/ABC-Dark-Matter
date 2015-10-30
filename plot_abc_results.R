compare.posterior <- function(parameters,weights,title,plot.seq,col.seq,legend.x,
                              legend.y,ylim,xlab,post.x,post.y,mle){
  par(mfrow=c(1,2), oma = c(0,0,0,0))
  n <- length(parameters)
  m <- length(plot.seq)
  legend.labels <- paste("Time Step",as.character(plot.seq))
  
  plot(density(parameters[[plot.seq[1]]], weights = weights[[1]]), ylim = ylim, xlab = xlab,
       main = "", col = col.seq[1], lwd = 2, ylab = "")
  lines(post.x,post.y, col = 4, lwd = 2)
  abline(v=mle, lty = 2, lwd = 2)
  legend(legend.x,legend.y,legend=c("ABC Posterior", "True Posterior", "MLE"),
         col = c(col.seq[1],4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)
  
  plot(density(parameters[[plot.seq[1]]], weights = weights[[1]]), ylim = ylim, xlab = xlab,
       main = "", col = 3, lwd = 2, ylab = "")
  for (ii in 2:m){
    lines(density(parameters[[plot.seq[ii]]], weights = weights[[ii]]), col = col.seq[ii], lwd = 2)
  }
  abline(v=mle, lty = 2, lwd = 2)
  legend(legend.x,legend.y,legend=c(legend.labels, "MLE"),
         col = c(col.seq,1), lty = c(rep(1,m),2), lwd=rep(1,m+1), cex = .6)
  mtext(title, outer = TRUE, line = -2)
  par(mfrow=c(1,1), oma = c(0,0,0,0))
  return()
}

##############################################################################
#Plot ABC Results
##############################################################################
setwd("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code")
load("posterior_calculations/posterior_Ec_results.R")
library(latex2exp)

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2dens3_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])
Ec.dist <- lapply(Ec.results, function(x) x[["distances"]])

pdf("figures/true_posterior_comparison.pdf")
plot(Ec.grid, Ec.posterior.reduced, col = 2, "l", xlim = c(.14, .177),
     xlab = latex2exp('$E_c$'), ylab = "Posterior Density",
     main = "True Posterior Comparison")
lines(Ec.grid, Ec.posterior.full, col = 4)
legend(.140, 240, legend = c("Precise Model", "Fast Model"),
       col = c(4, 2), lty = c(1,1))
dev.off()

pdf("figures/fast_model.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = latex2exp('$rv^2$ and 2D Density Posterior Comparison - Fast Model'),
                  plot.seq = c(26,25,20,15,10,5), col.seq = c(3, 4, 2, 5, 6, 7),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior.reduced, mle = Ec.mle.reduced$par)

plot(1:length(Ec.parameters),unlist(lapply(Ec.parameters, function(x) max(density(x)$y))),
     "l", xlab = "Time Step", ylab = "Maximum Density", main = "Fast Model")
dev.off()

pdf("figures/precise_model.pdf")
load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2densslow_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])
Ec.dist <- lapply(Ec.results, function(x) x[["distances"]])

compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = latex2exp('$rv^2$ and 2D Density Posterior Comparison - Precise Model'),
                  plot.seq = c(19, 18, 17, 15, 10, 5), col.seq = c(3, 4, 2, 5, 6, 7),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior.full, mle = Ec.mle.full$par)

plot(1:length(Ec.parameters),unlist(lapply(Ec.parameters, function(x) max(density(x)$y))),
     "l", xlab = "Time Step", ylab = "Maximum Density", main = "Precise Model")
dev.off()