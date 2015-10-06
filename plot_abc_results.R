##############################################################################
#Plot ABC Results
##############################################################################

load("posterior_calculations/posterior_Ec_results.R")
load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ec_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])

load.file <- "/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ec_output/results.R"
save.file <- "figures/KS_posteriors.pdf"

title <- "KS Statistic Posterior Comparison"
plot.seq <- c(30,25,20,10,15,10)
legend.x <- .1596
legend.y <- 250
ylim <- c(0,250)
xlab <- "Ec"
post.x <- Ec.grid
post.y <- Ec.posterior
mle <- Ec.mle.full$par

compare.posterior <- function(load.file,save.file,plot.seq,title){
  load(load.file)
  parameters <- lapply(results, function(x) x[["parameters"]])
  pdf(save.file)
  par(mfrow=c(1,2), oma = c(0,0,1,0))
  n <- length(parameters)
  legend.labels <- paste("Time Step",as.character(c(30,25,20,10,15,10)))
  
  plot(density(parameters[[n]]), ylim = ylim, xlab = xlab,
       main = "", col = 3, lwd = 2, ylab = "")
  lines(post.x,post.y, col = 4, lwd = 2)
  abline(v=mle, lty = 2, lwd = 2)
  legend(legend.x,legend.y,legend=c("ABC Posterior", "True Posterior", "MLE"),
         col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)
  
  plot(density(parameters[[n]]), ylim = ylim, xlab = xlab,
       main = "", col = 3, lwd = 2, ylab = "")
  lines(density(Ec.parameters[[25]]), col = 4, lwd = 2)
  lines(density(Ec.parameters[[20]]), col = 2, lwd = 2)
  lines(density(Ec.parameters[[15]]), col = 5, lwd = 2)
  lines(density(Ec.parameters[[10]]), col = 6, lwd = 2)
  abline(v=mle, lty = 2, lwd = 2)
  legend(legend.x,legend.y,legend=c(legend.labels, "MLE"),
         col = c(3,4,2,seq(5,m+1),1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
  mtext("KS Statistic Posterior Comparison", outer = TRUE, line = -1)
  par(mfrow=c(1,1), oma = c(0,0,0,0))
}

pdf(file="figures/KS_posteriors.pdf")
par(mfrow=c(1,2), oma = c(0,0,1,0))
plot(density(Ec.parameters[[30]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(Ec.grid,Ec.posterior, col = 4, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("ABC Posterior", "True Posterior", "MLE"),
       col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)

plot(density(Ec.parameters[[30]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(density(Ec.parameters[[25]]), col = 4, lwd = 2)
lines(density(Ec.parameters[[20]]), col = 2, lwd = 2)
lines(density(Ec.parameters[[15]]), col = 5, lwd = 2)
lines(density(Ec.parameters[[10]]), col = 6, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("Time Step 30", "Time Step 25", "Time Step 20",
                          "Time Step 15", "Time Step 10", "MLE"),
       col = c(3,4,2,5,6,1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
mtext("KS Statistic Posterior Comparison", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ecs_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])

pdf(file="figures/Smooth_KS_posteriors.pdf")
par(mfrow=c(1,2), oma = c(0,0,1,0))
plot(density(Ec.parameters[[20]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(Ec.grid,Ec.posterior, col = 4, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("ABC Posterior", "True Posterior", "MLE"),
       col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)

plot(density(Ec.parameters[[20]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(density(Ec.parameters[[16]]), col = 4, lwd = 2)
lines(density(Ec.parameters[[12]]), col = 2, lwd = 2)
lines(density(Ec.parameters[[8]]), col = 5, lwd = 2)
lines(density(Ec.parameters[[4]]), col = 6, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("Time Step 20", "Time Step 16", "Time Step 12",
                          "Time Step 8", "Time Step 4", "MLE"),
       col = c(3,4,2,5,6,1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
mtext("Smoothed KS Statistic Posterior Comparison", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/Ecl2_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])

pdf(file="figures/L2_posteriors.pdf")
par(mfrow=c(1,2), oma = c(0,0,1,0))
plot(density(Ec.parameters[[19]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(Ec.grid,Ec.posterior, col = 4, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("ABC Posterior", "True Posterior", "MLE"),
       col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)

plot(density(Ec.parameters[[19]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(density(Ec.parameters[[16]]), col = 4, lwd = 2)
lines(density(Ec.parameters[[12]]), col = 2, lwd = 2)
lines(density(Ec.parameters[[8]]), col = 5, lwd = 2)
lines(density(Ec.parameters[[4]]), col = 6, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("Time Step 19", "Time Step 16", "Time Step 12",
                          "Time Step 8", "Time Step 4", "MLE"),
       col = c(3,4,2,5,6,1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
mtext("L2 Norm Statistic Posterior Comparison", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/rv2_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])

pdf(file="figures/rv2_posteriors.pdf")
par(mfrow=c(1,2), oma = c(0,0,1,0))
plot(density(Ec.parameters[[30]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(Ec.grid,Ec.posterior, col = 4, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("ABC Posterior", "True Posterior", "MLE"),
       col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)

plot(density(Ec.parameters[[30]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(density(Ec.parameters[[25]]), col = 4, lwd = 2)
lines(density(Ec.parameters[[20]]), col = 2, lwd = 2)
lines(density(Ec.parameters[[15]]), col = 5, lwd = 2)
lines(density(Ec.parameters[[10]]), col = 6, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("Time Step 20", "Time Step 16", "Time Step 12",
                          "Time Step 8", "Time Step 4", "MLE"),
       col = c(3,4,2,5,6,1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
mtext("rv^2 Posterior Comparison", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/abc_runs/dens2d_output/results.R")
Ec.results <- results
rm(results)
Ec.parameters <- lapply(Ec.results, function(x) x[["parameters"]])

pdf(file="figures/dens2d_posteriors.pdf")
par(mfrow=c(1,2), oma = c(0,0,1,0))
plot(density(Ec.parameters[[26]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(Ec.grid,Ec.posterior, col = 4, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("ABC Posterior", "True Posterior", "MLE"),
       col = c(3,4,1), lty = c(1,1,2), lwd=c(2,2,2), cex = .6)

plot(density(Ec.parameters[[20]]), ylim = c(0,250), xlab = "Ec",
     main = "", col = 3, lwd = 2, ylab = "")
lines(density(Ec.parameters[[16]]), col = 4, lwd = 2)
lines(density(Ec.parameters[[12]]), col = 2, lwd = 2)
lines(density(Ec.parameters[[8]]), col = 5, lwd = 2)
lines(density(Ec.parameters[[4]]), col = 6, lwd = 2)
abline(v=Ec.mle.full$par, lty = 2, lwd = 2)
legend(.1596,250,legend=c("Time Step 20", "Time Step 16", "Time Step 12",
                          "Time Step 8", "Time Step 4", "MLE"),
       col = c(3,4,2,5,6,1), lty = c(1,1,1,1,1,2), lwd=c(2,2,2,2,2,2), cex = .6)
mtext("2d Density Posterior Comparison", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()