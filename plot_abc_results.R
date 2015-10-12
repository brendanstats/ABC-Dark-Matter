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

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ec_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/KS_posteriors_new.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = "KS Statistic Posterior Comparison",
                  plot.seq = c(30,25,20,15,10), col.seq = c(3, 4, 2, 5, 6),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ecs_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/Smooth_KS_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = "Smoothed KS Statistic Posterior Comparison",
                  plot.seq = c(20,16,12,8,4), col.seq = c(3, 4, 2, 5, 6),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ecl2_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/L2_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = latex2exp('$L^2$ Norm Statistic Posterior Comparison'),
                  plot.seq = c(20,16,12,8,4), col.seq = c(3, 4, 2, 5, 6),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()


load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/rv2_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = latex2exp('$rv^2$ Posterior Comparison'),
                  plot.seq = c(30,25,20,15,10), col.seq = c(3, 4, 2, 5, 6),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()


load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/dens2d_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/dens2d_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = "2D Density Posterior Comparison",
                  plot.seq = c(30,25,20,15,10), col.seq = c(3, 4, 2, 5, 6),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()


load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2_output/results_add.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/rv2_add_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights, 
                  title = latex2exp('$rv^2$ Posterior Comparison'),
                  plot.seq = c(36,33,25,20,15,10), col.seq = c(3, 4, 2, 5, 6,7),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()
plot(1:36,log(unlist(lapply(Ec.parameters,function(x) max(density(x)$y)))))
plot(1:36,unlist(lapply(Ec.parameters,function(x) max(density(x)$y))))


load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2dens_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])
Ec.weights <- lapply(Ec.results, function(x) x[["weights"]])

pdf(file="figures/rv2dens_posteriors.pdf")
compare.posterior(parameters = Ec.parameters, weights = Ec.weights,
                  title = latex2exp('$rv^2$ and 2D Density Posterior Comparison'),
                  plot.seq = c(22,21,20,15,10,5), col.seq = c(3, 4, 2, 5, 6,7),
                  legend.x = .1596, legend.y = 250, ylim = c(0,250),
                  xlab = latex2exp('$E_c$'), post.x = Ec.grid,
                  post.y = Ec.posterior, mle = Ec.mle.full$par)
dev.off()

plot(1:22,log(unlist(lapply(Ec.parameters,function(x) max(density(x)$y)))))
plot(1:22,unlist(lapply(Ec.parameters,function(x) max(density(x)$y))))
hist(Ec.weights[[22]])
