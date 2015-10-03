###############################################################################
#Test timing
###############################################################################

library(MASS)
library(rPython)
library(ggplot2)
python.load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/sampling_functions/functions_0721.py")

# time.eval <- function(param,samplesize,steps,resamplefactor){
#   start <- Sys.time()
#   python.call("sample3R",param2,steps,samplesize,resamplefactor)
#   end <- Sys.time()
#   return(end-start)
# }
# 
# #test sampling times
# param <- c(2.0,-5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5,1,3,1)
# out <- list()
# n <- 1
# for(rlim in seq(1,10,1)){
#   param[5] <- rlim
#   for(Ec in seq(.11,.21,.02)){
#     param[4] <- Ec
#     param <- as.numeric(param)
#     for(steps in 10:25){
#       times <- replicate(n = 10,time.eval(param=param,samplesize=samplesize,steps=steps,resamplefactor=resamplefactor))
#       out[[n]] <- list(rlim=rlim,Ec=Ec,steps=steps,times=times)
#       n <- n + 1
#       save(out,file = "/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/sampling_functions/sample3R_times.R")
#       cat(rlim, Ec, steps,Sys.time(),"\n")
#     }
#   }
# }

###############################################################################
#Test new distance functions
###############################################################################

#######################################
#KS Test Statistc
#######################################

ks.dist <- function(data1,data2){
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) ks.test(data1[,n],data2[,n])$statistic,data1,data2))
  ks.r <- sqrt(sum(ks.vec[1:3]^2))
  ks.v <- sqrt(sum(ks.vec[4:6]^2))
  return(c(ks.r,ks.v))
}

ks.smooth <- function(x,y){
  x.dens <- density(x)
  n <- length(x.dens$x)
  y.dens <- density(y, from = x.dens$x[1], to = x.dens$x[n], n=n)
  x.cum <- cumsum(x.dens$y)/sum(x.dens$y)
  y.cum <- cumsum(y.dens$y)/sum(y.dens$y)
  return(max(abs(x.cum-y.cum)))
}

ks.dist.sm <- function(data1,data2){
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) ks.smooth(data1[,n],data2[,n]),data1,data2))
  ks.r <- sqrt(sum(ks.vec[1:3]^2))
  ks.v <- sqrt(sum(ks.vec[4:6]^2))
  return(c(ks.r,ks.v))
}

#######################################
#L2 distance
#######################################

l2.smooth <- function(x1,x2){
  n <- length(x1)
  dens1 <- density(x1)
  dens2 <- density(x2, from = min(dens1$x), to = max(dens1$x), n = length(dens1$x))
  y1 <- cumsum(dens1$y)/sum(dens1$y)
  y2 <- cumsum(dens2$y)/sum(dens2$y)
  return(sqrt(sum((y1-y2)^2)))
}

l2.dist <- function(data1,data2){
  l2.vec <- unlist(lapply(1:6,function(n,data1,data2) l2.smooth(data1[,n],data2[,n]),data1,data2))
  l2.r <- sqrt(sum(l2.vec[1:3]^2))
  l2.v <- sqrt(sum(l2.vec[4:6]^2))
  return(c(l2.r,l2.v))
}

#######################################
#AD Test Statistc
#######################################

ad.stat <- function(x,y, dist.min, dist.max){
  nx <- length(x)
  ny <- length(y)
  
  x.sort <- sort(x)
  y.sort <- sort(y)
  
  dens.x <- density(x.sort)
  Fxx <- c(dist.min,dens.x$x,dist.max)
  Fxy <- c(0,dens.x$y/sum(sum(dens.x$y)),1)
  Fy <- approx(Fxx, Fxy, xout = y.sort)
  
  S = sum((2*1:ny-1)/ny*(log(Fy$y)+log(1-rev(Fy$y))))
  return(-ny-S)
}

ad.distance <- function(data1,data2){
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) ad.stat(data1[,n],data2[,n]),data1,data2))
  return(sqrt(sum( ks.vec^2)))
}

data1 <- data.true.E
data2 <- Ec.galaxies[[1]][[2]]

dens1 <- density(data1[,1])
dens2 <- density(data2[,1])
sqrt(sum((dens1$y - approx(dens2$x,dens2$y,xout=dens1$x)$y)^2))

for(ii in 1:2000){
  if(is.na(ad.stat(data.true[,1],Ec.galaxies[[ii]][["sample"]][,1]))) print(ii)
  if(is.na(ad.stat(data.true[,2],Ec.galaxies[[ii]][["sample"]][,2]))) print(ii)
  if(is.na(ad.stat(data.true[,3],Ec.galaxies[[ii]][["sample"]][,3]))) print(ii)
  
  if(is.na(ad.stat(data.true[,4],Ec.galaxies[[ii]][["sample"]][,4]))) print(ii)
  if(is.na(ad.stat(data.true[,5],Ec.galaxies[[ii]][["sample"]][,5]))) print(ii)
  if(is.na(ad.stat(data.true[,6],Ec.galaxies[[ii]][["sample"]][,6]))) print(ii)
}

#######################################
#Difference between maximum radius
#######################################

max.r <- function(data1,data2){
  max1 <- max(apply(data1[1:3,],1,function(x) sqrt(sum(x^2))))
  max2 <- max(apply(data2[1:3,],1,function(x) sqrt(sum(x^2))))
  return(abs(max1 - max2))
}

#######################################
#Star densities
#######################################

radius.densities <- function(quantiles,r1,r2){
  
  n <- length(quantiles)
  
  cutoffs1 <- c(0,quantile(r1, probs=quantiles))
  cutoffs2 <- c(0,quantile(r1, probs=quantiles))
  
  dens1 <- unlist(lapply(1:n+1, function(k,x,cutoffs) sum(x<=cutoffs[k] & x > cutoffs[k-1])^2/(cutoffs[k]^3 - cutoffs[k-1]^3),r1,cutoffs1))
  dens2 <- unlist(lapply(1:n+1, function(k,x,cutoffs) sum(x<=cutoffs[k] & x > cutoffs[k-1])^2/(cutoffs[k]^3 - cutoffs[k-1]^3),r2,cutoffs2))
  
  return(sqrt(sum((dens1-dens2)^2)))
}

#######################################
#Velocity Dispersion
#######################################

velocity.dispertion <- function(quantiles,v1,v2){
  
  n <- length(quantiles)
  
  cutoffs1 <- c(0,quantile(r1, probs=quantiles))
  cutoffs2 <- c(0,quantile(r1, probs=quantiles))
  
  disp1 <- unlist(lapply(1:n+1, function(k,x,cutoffs) var(x[x<=cutoffs[k] & x > cutoffs[k-1]]),v1,cutoffs1))
  disp2 <- unlist(lapply(1:n+1, function(k,x,cutoffs) var(x[x<=cutoffs[k] & x > cutoffs[k-1]]),v2,cutoffs2))
  
  return(sqrt(sum((disp1-disp2)^2)))
}


#######################################

#######################################

rv2.dist <- function(data1,data2){
  r1 <- apply(data1[,1:3],1,function(y) sqrt(sum(y^2)))
  v1 <- apply(data1[,4:6],1,function(y) sqrt(sum(y^2)))
  
  r2 <- apply(data2[,1:3],1,function(y) sqrt(sum(y^2)))
  v2 <- apply(data2[,4:6],1,function(y) sqrt(sum(y^2)))
  
  x1 <- r1*(v1^2)
  x2 <- r2*(v2^2)
  
  return(l2.smooth(x1,x2))
}

dens2d.dist <- function(data1,data2,n=25){
  r1 <- apply(data1[,1:3],1,function(y) sqrt(sum(y^2)))
  v1 <- apply(data1[,4:6],1,function(y) sqrt(sum(y^2)))
  r2 <- apply(data2[,1:3],1,function(y) sqrt(sum(y^2)))
  v2 <- apply(data2[,4:6],1,function(y) sqrt(sum(y^2)))
  
  h <- c(bandwidth.nrd(r1),bandwidth.nrd(v1))
  
  base <- kde2d(r1,v1,h=h,n=n,lims=c(range(r1),range(v1)))
  test <- kde2d(r2,v2,h=h,n=n,lims=c(range(r1),range(v1)))

  return(sum((base$z - test$z)^2))
}

###############################################################################
#Plot test stat values against samples
###############################################################################

load("posterior_calculations/posterior_Ec_results.R")

partial.true <- read.delim("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/data/simulated_galaxy_reduced_model.txt")
load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/sim_galaxies/samples_Ec.R")
Ec.list <- sample
rm(sample)

partial.r <- apply(partial.true[,1:3],1,function(y) sqrt(sum(y^2)))
partial.v <- apply(partial.true[,4:6],1,function(y) sqrt(sum(y^2)))

Ec.param <- unlist(lapply(Ec.list, function(x,item) x[[item]],"theta"))
Ec.samples <- lapply(Ec.list, function(x,item) x[[item]],"sample")
Ec.samples.r <- lapply(Ec.samples, function(x) apply(x[,1:3],1,function(y) sqrt(sum(y^2))))
Ec.samples.v <- lapply(Ec.samples, function(x) apply(x[,4:6],1,function(y) sqrt(sum(y^2))))

partial.ks.dist <- lapply(Ec.samples, function(x,true) ks.dist(true,x),partial.true)
partial.ks.dist <- do.call(rbind,partial.ks.dist)
partial.ks.smooth.dist <- lapply(Ec.samples, function(x,true) ks.dist.sm(true,x),partial.true)
partial.ks.smooth.dist <- do.call(rbind,partial.ks.smooth.dist)

pdf(file = "Figures/KS_Distance_Comparison.pdf")
par(mfrow=c(2,2), oma = c(0,0,1,0))
plot(Ec.param,partial.ks.dist[,1], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Radius KS Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.smooth.dist[,1], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Radius KS Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.dist[,2], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Velocity KS Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.smooth.dist[,2], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Velocity KS Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
mtext("KS Distance Smooth vs. Unsmoothed", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

#Put in single r and v numbers
partial.ks.r <- unlist(lapply(Ec.samples.r, function(x,true) ks.test(true,x)$statistic,partial.r))
partial.ks.v <- unlist(lapply(Ec.samples.v, function(x,true) ks.test(true,x)$statistic,partial.v))
partial.ks.smooth.r <- unlist(lapply(Ec.samples.r, function(x,true) ks.smooth(true,x),partial.r))
partial.ks.smooth.v <- unlist(lapply(Ec.samples.v, function(x,true) ks.smooth(true,x),partial.v))

pdf(file = "Figures/KS_Stat_Comparison.pdf")
par(mfrow=c(2,2), oma = c(0,0,1,0))
plot(Ec.param,partial.ks.r, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Radius KS Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.smooth.r, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Radius KS Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.v, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Velocity KS Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.ks.smooth.v, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Velocity KS Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
mtext("KS Statistic Smooth vs. Unsmoothed", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

partial.l2.dist <- lapply(Ec.samples, function(x,true) l2.dist(true,x),partial.true)
partial.l2.dist <- do.call(rbind,partial.l2.dist)
partial.l2.r <- unlist(lapply(Ec.samples.r, function(x,true) l2.smooth(true,x),partial.r))
partial.l2.v <- unlist(lapply(Ec.samples.v, function(x,true) l2.smooth(true,x),partial.v))

pdf(file = "Figures/L2_Comparison.pdf")
par(mfrow=c(2,2), oma = c(0,0,1,0))
plot(Ec.param,partial.l2.dist[,1], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Radius L2 Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.l2.r, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Radius L2 Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.l2.dist[,2], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Velocity L2 Distance")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
plot(Ec.param,partial.l2.v, pch = 19, cex = .3, xlab = "Ec",
     ylab = "Velocity L2 Stat")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
mtext("L2 Distance vs. Statistic", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))
dev.off()

partial.rv2.dist <- lapply(Ec.samples, function(x,true) rv2.dist(true,x),partial.true)
partial.rv2.dist <- do.call(rbind,partial.rv2.dist)

partial.dens2d.dist25 <- lapply(Ec.samples, function(x,true) dens2d.dist(true,x,25),partial.true)
partial.dens2d.dist25 <- do.call(rbind,partial.dens2d.dist25)

partial.dens2d.dist50 <- lapply(Ec.samples, function(x,true) dens2d.dist(true,x,50),partial.true)
partial.dens2d.dist50 <- do.call(rbind,partial.dens2d.dist50)

partial.dens2d.dist2550 <- lapply(Ec.samples, function(x,true) dens2d.dist(true,x,c(25,50)),partial.true)
partial.dens2d.dist2550 <- do.call(rbind,partial.dens2d.dist2550)

pdf("figures/rv2_Ec.pdf")
plot(Ec.param,partial.rv2.dist, pch = 19, cex = .3, xlab = "Ec",
     ylab = "rv^2", main = "rv^2 L2 distance vs. Ec")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
dev.off()

pdf("figures/2D_density_Ec.pdf")
par(mfrow=c(2,2))
plot(Ec.param,partial.dens2d.dist25, pch = 19, cex = .3, xlab = "Ec",
     ylab = "2D Density Distance", main = "n = (25,25)")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)

plot(Ec.param,partial.dens2d.dist50, pch = 19, cex = .3, xlab = "Ec",
     ylab = "2D Density Distance", main = "n = (50,50)")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)

plot(Ec.param,partial.dens2d.dist2550, pch = 19, cex = .3, xlab = "Ec",
     ylab = "2D Density Distance", main = "n = (25,50)")
abline(v=Ec.mle.full$par, lty = 2, lwd = 1)
par(mfrow=c(1,1))
dev.off()


# rdens.distances.Ec <- unlist(lapply(Ec.galaxies.r, function(r1,r2,quantiles) radius.densities(quantiles,r1,r2),true.E.r, seq(.1,1,.1)))
# vdisp.distances.Ec <- lapply(Ec.galaxies.v, function(v1,v2,quantiles) radius.densities(quantiles,v1,v2),true.E.v, seq(.1,1,.1))

# ad.distances.Ec <- unlist(lapply(Ec.galaxies.samples, function(x,true) ad.mag(true,x),data.true.E))
# rmax.distances.Ec <- lapply(Ec.galaxies.samples, max.r, data.true.E)
#rlim.galaxies <- load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/sim_galaxies/samples_rlim.R")


quantile.indicies.single <- function(y,quantiles){
  n <- length(quantiles)
  cutoffs <- quantile(y,quantiles)
  indicies <- lapply(1:n, function(n,cutoffs,y) y<cutoffs[n],cutoffs,y)
  for(ii in n:2) {
    indicies[[ii]] <- as.logical(indicies[[ii]] - indicies[[ii-1]])
  }
  return(indicies)
}

idx <- quantile.indicies.single(partial.l2.r,seq(.1,1,length.out = 10))

quantile.indicies.multiple <- function(y,quantiles){
  n <- length(quantiles)
  cutoffs <- apply(y,2,quantile,quantiles)
  indicies <- lapply(1:n, function(n,cutoffs,y) t(apply(y,1,function(y,cut) y<cut,cutoffs[n,])),cutoffs,y)
  indicies <- lapply(indicies, function(x) apply(x,1,all))
  for(ii in n:2) {
    indicies[[ii]] <- as.logical(indicies[[ii]] - indicies[[ii-1]])
  }
  return(indicies)
}

full.true <- read.delim("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/data/simulated_galaxy_full_model.txt")
load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/sim_galaxies/samples3_Ec_rlim.R")
Ec.rlim.list <- sample
rm(sample)
#2089 is weird...

full.r <- apply(full.true[,1:3],1,function(y) sqrt(sum(y^2)))
full.v <- apply(full.true[,4:6],1,function(y) sqrt(sum(y^2)))

Ec.rlim.param <- lapply(Ec.rlim.list, function(x,item) x[[item]],"theta")
Ec.rlim.param <- do.call(rbind,Ec.rlim.param)
Ec.rlim.samples <- lapply(Ec.rlim.list, function(x,item) x[[item]],"sample")
Ec.rlim.samples.r <- lapply(Ec.rlim.samples, function(x) apply(x[,1:3],1,function(y) sqrt(sum(y^2))))
Ec.rlim.samples.v <- lapply(Ec.rlim.samples, function(x) apply(x[,4:6],1,function(y) sqrt(sum(y^2))))

full.ks.dist <- lapply(Ec.rlim.samples, function(x,true) ks.dist(true,x),full.true)
full.ks.dist <- do.call(rbind,full.ks.dist)
full.ks.smooth.dist <- lapply(Ec.rlim.samples, function(x,true) ks.dist.sm(true,x),full.true)
full.ks.smooth.dist <- do.call(rbind,full.ks.smooth.dist)

pdf(file = "Figures/KS_Ecrlim_Comparison.pdf")
idx <- quantile.indicies.single(full.ks.dist[,1],seq(.1,1,length.out = 10))
title <- "KS Distance Radius Component"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)

idx <- quantile.indicies.single(full.ks.dist[,2],seq(.1,1,length.out = 10))
title <- "KS Distance Velocity Component"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)

idx <- quantile.indicies.multiple(full.ks.dist,seq(.1,1,length.out = 10))
title <- "KS Distance Radius and Velocity Components"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)
dev.off()

pdf(file = "Figures/KS_smooth_Ecrlim_Comparison.pdf")
idx <- quantile.indicies.single(full.ks.smooth.dist[,1],seq(.1,1,length.out = 10))
title <- "KS (Smoothed) Distance Radius Component"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)

idx <- quantile.indicies.single(full.ks.smooth.dist[,2],seq(.1,1,length.out = 10))
title <- "KS (Smoothed) Distance Velocity Component"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)

idx <- quantile.indicies.multiple(full.ks.smooth.dist,seq(.1,1,length.out = 10))
title <- "KS (Smoothed) Distance Radius and Velocity Components"
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)
dev.off()

full.rv2.dist <- lapply(Ec.rlim.samples, function(x,true) rv2.dist(true,x),full.true)
full.rv2.dist <- do.call(rbind,full.rv2.dist)
idx <- quantile.indicies.single(full.rv2.dist,seq(.1,1,length.out = 10))
title <- "RV^2 Distance"
pdf(file = "figures/rv2_Ecrlim.pdf")
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)
dev.off()

full.dens2d.dist <- lapply(Ec.rlim.samples, function(x,true) dens2d.dist(true,x,25),full.true)
full.dens2d.dist <- do.call(rbind,full.dens2d.dist)
idx <- quantile.indicies.single(full.dens2d.dist,seq(.1,1,length.out = 10))
title <- "2D Density Distance"

pdf(file = "figures/2D_density_Ecrlim.pdf")
plot(Ec.rlim.param[idx[[1]],1],Ec.rlim.param[idx[[1]],2], pch = 19, cex = .3,
     xlab = "Ec", ylab = "rlim", main = title)
points(Ec.rlim.param[idx[[2]],1],Ec.rlim.param[idx[[2]],2], pch = 19, cex = .3, col = 2)
points(Ec.rlim.param[idx[[3]],1],Ec.rlim.param[idx[[3]],2], pch = 19, cex = .3, col = 3)
points(Ec.rlim.param[idx[[4]],1],Ec.rlim.param[idx[[4]],2], pch = 19, cex = .3, col = 4)
points(Ec.rlim.param[idx[[5]],1],Ec.rlim.param[idx[[5]],2], pch = 19, cex = .3, col = 5)
points(Ec.rlim.param[idx[[6]],1],Ec.rlim.param[idx[[6]],2], pch = 19, cex = .3, col = 6)
points(Ec.rlim.param[idx[[7]],1],Ec.rlim.param[idx[[7]],2], pch = 19, cex = .3, col = 7)
points(Ec.rlim.param[idx[[8]],1],Ec.rlim.param[idx[[8]],2], pch = 19, cex = .3, col = 8)
legend(.70,12, legend = c("Smallest 10%", "Smallest 20%", "Smallest 30%","Smallest 40%" ,
                          "Smallest 50%", "Smallest 60%", "Smallest 70%", "Smallest 80%"),
       pch = rep(19,8),col = 1:8)
dev.off()

par(mfrow=c(2,2), oma = c(0,0,1,0))
plot(Ec.rlim.param[,1],full.ks.dist[,1], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Radius KS Distance")
plot(Ec.rlim.param[,1],full.ks.smooth.dist[,1], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Radius KS Distance")
plot(Ec.rlim.param[,1],full.ks.dist[,2], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Velocity KS Distance")
plot(Ec.rlim.param[,1],full.ks.smooth.dist[,2], pch = 19, cex = .3, xlab = "Ec",
     ylab = "Smoothed Velocity KS Distance")
mtext("KS Distance Smooth vs. Unsmoothed", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))


par(mfrow=c(2,2), oma = c(0,0,1,0))
plot(Ec.rlim.param[,2],full.ks.dist[,1], pch = 19, cex = .3, xlab = "rlim",
     ylab = "Radius KS Distance")
plot(Ec.rlim.param[,2],full.ks.smooth.dist[,1], pch = 19, cex = .3, xlab = "rlim",
     ylab = "Smoothed Radius KS Distance")
plot(Ec.rlim.param[,2],full.ks.dist[,2], pch = 19, cex = .3, xlab = "rlim",
     ylab = "Velocity KS Distance")
plot(Ec.rlim.param[,2],full.ks.smooth.dist[,2], pch = 19, cex = .3, xlab = "rlim",
     ylab = "Smoothed Velocity KS Distance")
mtext("KS Distance Smooth vs. Unsmoothed", outer = TRUE, line = -1)
par(mfrow=c(1,1), oma = c(0,0,0,0))

ks.distances.rlim.Ec <- unlist(lapply(rlim.Ec.galaxies.samples, function(x,true) ks.mag(true,x),data.true))
#ad.distances.rlim.Ec <- unlist(lapply(rlim.Ec.galaxies.samples, function(x,true) ad.mag(true,x),data.true))
#rmax.distances.rlim.Ec <- lapply(rlim.Ec.galaxies.samples, max.r, data.true)
rdens.distances.rlim.Ec <- lapply(rlim.Ec.galaxies.r, function(r1,r2,quantiles) radius.densities(quantiles,r1,r2),true.E.r, seq(.1,1,.1))
vdisp.distances.rlim.Ec <- lapply(rlim.Ec.galaxies.v, function(v1,v2,quantiles) radius.densities(quantiles,v1,v2),true.E.v, seq(.1,1,.1))

pdf(file="../figures/Ecrlim_ks_scatter.pdf")
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
plot(rlim.Ec.galaxies.parameters[,1],ks.distances.rlim.Ec, cex = .3, xlab = "Ec",
     ylab = "KS Stat")
plot(rlim.Ec.galaxies.parameters[,2],ks.distances.rlim.Ec, cex = .3, xlab = "rlim",
     ylab = "KS Stat")
mtext("KS Stat vs. Single Parameters", outer = TRUE, line = 0)
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
dev.off()

df <- data.frame(Ec=rlim.Ec.galaxies.parameters[,1],
      rlim=rlim.Ec.galaxies.parameters[,2],
      ks=ks.distances.rlim.Ec)
df$group <- 1
df$group[df$ks>quantile(df$ks,.1)] <- 2
df$group[df$ks>quantile(df$ks,.2)] <- 3
df$group[df$ks>quantile(df$ks,.3)] <- 4
df$group[df$ks>quantile(df$ks,.4)] <- 5
df$group[df$ks>quantile(df$ks,.5)] <- 6
df$group[df$ks>quantile(df$ks,.6)] <- 7
df$group[df$ks>quantile(df$ks,.7)] <- 8
df$group[df$ks>quantile(df$ks,.8)] <- 9
df$group[df$ks>quantile(df$ks,.9)] <- 10


pdf(file="../figures/Ecrlim_ks_continuous.pdf")
ggplot(df, aes(x=Ec,y=rlim)) + geom_point(aes(col = ks))
dev.off()
#ggplot(df, aes(x=Ec,y=rlim)) + geom_point(aes(col = group))

cutoffs <- as.character(round(unlist(lapply(seq(.1,.9,.1), function(x) quantile(df$ks,x))),2))

pdf(file="../figures/Ecrlim_ks_quantile.pdf")
plot(df$Ec[df$group==1],df$rlim[df$group==1],pch=19,cex = .3,
     ylab = "rlim", xlab = "Ec")
points(df$Ec[df$group==2],df$rlim[df$group==2],pch=19,cex = .3,col=2)
points(df$Ec[df$group==3],df$rlim[df$group==3],pch=19,cex = .3,col=3)
points(df$Ec[df$group==4],df$rlim[df$group==4],pch=19,cex = .3,col=4)
points(df$Ec[df$group==5],df$rlim[df$group==5],pch=19,cex = .3,col=5)
points(df$Ec[df$group==6],df$rlim[df$group==6],pch=19,cex = .3,col=6)
points(df$Ec[df$group==7],df$rlim[df$group==7],pch=19,cex = .3,col=7)
points(df$Ec[df$group==8],df$rlim[df$group==8],pch=19,cex = .3,col=8)
legend(.8,12,legend=c(paste("ks stat <", cutoffs[1:8])),
       col = 1:8, pch = rep(19,8), cex = .8)
dev.off()


###############################################################################
#Accepted and rejected parameters
###############################################################################

load("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/2015_09_10/Ec_all_15_add_results.R")

step <- 7
sampled.values <- results[[step]][["all.tested"]][,1]
tested.values <- results[[step]][["all.tested"]][,2]
accepted.flag <- results[[step]][["all.tested"]][,5]

plot(sampled.values[!accepted.flag],tested.values[!accepted.flag], cex = .4, pch = 16, col = "red")
points(sampled.values[accepted.flag],tested.values[accepted.flag], cex = .4, pch = 16)

plot(density(sampled.values))

plot(density(tested.values))

plot(density(tested.values[accepted.flag]))
lines(density(tested.values[!accepted.flag]))
