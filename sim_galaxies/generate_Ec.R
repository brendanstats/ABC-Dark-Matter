#####################
#Advanced Data Analysis
#Brendan McVeigh
#July 28, 2015
#####################

setwd("~/Google Drive/2015_S1_Spring/ADA/2015_07_28")
library(rPython)
library(parallel)
library(doParallel) 
library(foreach)


python.load("functions_0624.py")
source("ABC_Code/priors.R")

##########################################
#Generate Datasets
##########################################

generate.sample <- function(python.file,sample.fun,param,ct,prior){
  
  python.load(python.file)
  
	proposal <- prior(1,"s")
	param[4] <- proposal
  param <- as.numeric(param)
  
	sample <- do.call(cbind,python.call(sample.fun,param,ct))
	out.pair <- list(theta=proposal,sample=sample)
	return(out.pair)
}

file.name <- "Ec_samples.R"
n <- 2000
ct <- 5000
param <- c(2.0,-5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)

ncore <- detectCores()
cl <- makeCluster(ncore)
registerDoParallel(cl)
sample <- foreach(i=1:n, .packages="rPython") %dopar%
  generate.sample(python.file = "functions_0624.py", 
  sample.fun = "sampleR",
  param = param,ct = ct, prior = prior.Ec)
stopCluster(cl)

save(sample,file = file.name)

##########################################
#Distance Functions to test
##########################################

mean.x <- function(data.sim){
  return(mean(data.sim[,1]))
}

mean.y <- function(data.sim){
  return(mean(data.sim[,2]))
}

mean.r <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  return(mean(r))
}

mean.vz <- function(data.sim){
  return(mean(data.sim[,6]))
}

m2.x <- function(data.sim){
  return(mean(data.sim[,1]*data.sim[,1]))
}

m2.y <- function(data.sim){
  return(mean(data.sim[,2]*data.sim[,2]))
}

m2.r <- function(data.sim){
  r2 <- data.sim[,1]^2+data.sim[,2]^2
  return(mean(r2))
}

m2.vz <- function(data.sim){
  return(mean(data.sim[,6]*data.sim[,6]))
}

m3.x <- function(data.sim){
  return(mean(data.sim[,1]*data.sim[,1]*data.sim[,1]))
}

m3.y <- function(data.sim){
  return(mean(data.sim[,2]*data.sim[,2]*data.sim[,2]))
}

m3.r <- function(data.sim){
  r3 <- (data.sim[,1]^2+data.sim[,2]^2)^(3/2)
  return(mean(r3))
}

m3.vz <- function(data.sim){
  return(mean(data.sim[,6]*data.sim[,6]*data.sim[,6]))
}

vdisp.25 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(0))
  max.r <- quantile(r,probs = c(.25))
  keep <- min.r <= r & r < max.r
  return(var(data.sim[keep,6]))
}

vdisp.5 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.25))
  max.r <- quantile(r,probs = c(.5))
  keep <- min.r <= r & r < max.r
  return(var(data.sim[keep,6]))
}

vdisp.75 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.5))
  max.r <- quantile(r,probs = c(.75))
  keep <- min.r <= r & r < max.r
  return(var(data.sim[keep,6]))
}

vdisp1 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.75))
  max.r <- quantile(r,probs = c(1))
  keep <- min.r <= r & r < max.r
  return(var(data.sim[keep,6]))
}
dens.25 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(0))
  max.r <- quantile(r,probs = c(.25))
  num <- dim(data.sim)[1]/4
  return(num/(max.r^2-min.r^2))	
}

dens.5 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.25))
  max.r <- quantile(r,probs = c(.5))
  num <- dim(data.sim)[1]/4
  return(num/(max.r^2-min.r^2))	
}

dens.75 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.5))
  max.r <- quantile(r,probs = c(.75))
  num <- dim(data.sim)[1]/4
  return(num/(max.r^2-min.r^2))	
}

dens1 <- function(data.sim){
  r <- sqrt(data.sim[,1]^2+data.sim[,2]^2)
  min.r <- quantile(r,probs = c(.75))
  max.r <- quantile(r,probs = c(1))
  num <- dim(data.sim)[1]/4
  return(num/(max.r^2-min.r^2))	
}

load("distance_samples_test.R")

thetas <- unlist(lapply(sample, function(x) x[["theta"]]))

mu.x <- unlist(lapply(sample,function(x) mean.x(x$sample)))
plot(thetas, mu.x, xlab = "Ec", ylab = "Mean of x coordinate")

mu.y <- unlist(lapply(sample,function(x) mean.y(x$sample)))
plot(thetas, mu.y, xlab = "Ec", ylab = "Mean of y coordinate")

mu.r <- unlist(lapply(sample,function(x) mean.r(x$sample)))
plot(thetas,mu.r, xlab = "Ec", ylab = "Mean of radius")

mu.vz <- unlist(lapply(sample,function(x) mean.vz(x$sample)))
plot(thetas,mu.vz, xlab = "Ec", ylab = "Mean of z velocity")

mom2.x <- unlist(lapply(sample,function(x) m2.x(x$sample)))
plot(thetas,mom2.x, xlab = "Ec", ylab = "Second Moment of x coordinate")

mom2.y <- unlist(lapply(sample,function(x) m2.y(x$sample)))
plot(thetas,mom2.y, xlab = "Ec", ylab = "Second Moment of y coordinate")

mom2.r <- unlist(lapply(sample,function(x) m2.r(x$sample)))
plot(thetas,mom2.r, xlab = "Ec", ylab = "Second Moment of radius")

mom2.vz <- unlist(lapply(sample,function(x) m2.vz(x$sample)))
plot(thetas,mom2.vz, xlab = "Ec", ylab = "Second Moment of z velocity")

mom3.x <- unlist(lapply(sample,function(x) m3.x(x$sample)))
plot(thetas,mom3.x, xlab = "Ec", ylab = "Third Moment of x coordinate")

mom3.y <- unlist(lapply(sample,function(x) m3.y(x$sample)))
plot(thetas,mom3.y, xlab = "Ec", ylab = "Third Moment of y coordinate")

mom3.r <- unlist(lapply(sample,function(x) m3.r(x$sample)))
plot(thetas,mom3.r, xlab = "Ec", ylab = "Third Moment of radius")

mom3.vz <- unlist(lapply(sample,function(x) m3.vz(x$sample)))
plot(thetas,mom3.vz, xlab = "Ec", ylab = "Third Moment of z velocity")

vd.25 <- unlist(lapply(sample,function(x) vdisp.25(x$sample)))
plot(thetas,vd.25, xlab = "Ec", ylab = "Sample Variance of z velocity",
     main = "Sample Variance of Velocity \n 1st quartile of radius")

vd.5 <- unlist(lapply(sample,function(x) vdisp.5(x$sample)))
plot(thetas,vd.5, xlab = "Ec", ylab = "Sample Variance of z velocity",
     main = "Sample Variance of Velocity \n 2nd quartile of radius")

vd.75 <- unlist(lapply(sample,function(x) vdisp.75(x$sample)))
plot(thetas,vd.75, xlab = "Ec", ylab = "Sample Variance of z velocity",
     main = "Sample Variance of Velocity \n 3rd quartile of radius")

vd1 <- unlist(lapply(sample,function(x) vdisp1(x$sample)))
plot(thetas,vd1, xlab = "Ec", ylab = "Sample Variance of z velocity",
     main = "Sample Variance of Velocity \n 4th quartile of radius")

rden.25 <- unlist(lapply(sample,function(x) dens.25(x$sample)))
plot(thetas,rden.25, xlab = "Ec", ylab = "Stars per unit area",
     main = "Sample Density \n 1st quartile of radius")

rden.5 <- unlist(lapply(sample,function(x) dens.5(x$sample)))
plot(thetas,rden.5, xlab = "Ec", ylab = "Stars per unit area",
     main = "Sample Density \n 2nd quartile of radius")

rden.75 <- unlist(lapply(sample,function(x) dens.75(x$sample)))
plot(thetas,rden.75, xlab = "Ec", ylab = "Stars per unit area",
     main = "Sample Density \n 3rd quartile of radius")

rden1 <- unlist(lapply(sample,function(x) dens1(x$sample)))
plot(thetas,rden1, xlab = "Ec", ylab = "Stars per unit area",
     main = "Sample Density \n 4th quartile of radius")

start <- Sys.time()
rden1 <- unlist(lapply(sample,function(x) dens1(x$sample)))
end <- Sys.time()
end-start

start <- Sys.time()
rden2 <- unlist(mclapply(sample,function(x) dens1(x$sample)))
end <- Sys.time()
end-start