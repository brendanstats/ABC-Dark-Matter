###################
#Advanced Data Analysis
#Brendan McVeigh
#May 25, 2015
###################

#######################################
#Define Distance Functions
#######################################

###################
#helper functions
###################
magnitude <- function(x){
  return(sqrt(sum(x^2)))
}

calc.m2 <- function(x){
  sum(x^2)/length(x)
}

#######################################
#functions from python
#######################################
getprojdensity <- function(x){
  n <- dim(x)[1]
  
  r <- apply(x,1,magnitude)
  r.sort <- sort(r)
  rhoR <- n:1/n
  return(cbind(r.sort,rhoR))
}

sd.mean <- function(bin,r,vz){
  rows <- 1:length(r)
  rows <- which(rows > bin[1] & rows <= bin[2])
  r.mean <- mean(r[rows])
  vz.sd <- sqrt(((length(rows)-1)/length(rows))*var(vz[rows]))
  return(c(r.mean,vz.sd))
}

getlosvdisp <- function(x,bins.num=30){
  r <- apply(x[,1:2],1,magnitude)
  ord <- order(r)
  r.sort <- r[ord]
  vz.sort <- x[ord,3]
  
  lower <- seq(0,length(r),length.out=bins.num)
  upper <- lower+lower[2]
  bins <- cbind(lower,upper)
  
  disp <- t(apply(bins,1,sd.mean,r.sort,vz.sort))
  return(disp)
}

#######################################
#Single Distance Functions
#######################################

###################
#Gaussian Distance, assumes sd of 2, sample size 1000
###################

mean.dist <- function(theta, suff.stats){
  sim.data <- rnorm(1000,theta,2)
  return(abs(mean(sim.data-suff.stats)))
}

###################
#Energy Distance
###################

energy.dist <- function(theta, suff.stats) {
  param <- as.numeric(c(theta,c(2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)))  
  
  n <- dim(suff.stats)[1]
  sim.data <- do.call(cbind,python.call("sample",param,n-1))
  pool <- rbind(suff.stats,sim.data[,c(1,2,6)])  
  pool.dist <- dist(pool)
  epsilon <- as.numeric(edist(pool.dist,c(n,n),distance=TRUE))
  return(epsilon)
}

###################
#Chi-Squared Galaxy Distance
###################

chi.sq.dist <- function(theta, suff.stats){
  param.1 <- as.numeric(c(theta[1:2],c(2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)))
  param.2 <- as.numeric(c(theta[3:4],c(1.1, 0.17, 1.5, 0, 8.2, 0.086, 21.0, 3.0)))
  
  n <- dim(data)[1]
  sim.data.1 <- do.call(cbind,python.call("sample",param.1,n-1))
  sim.data.2 <- do.call(cbind,python.call("sample",param.2,n-1))
  
  #proj.den.1 <- getprojdensity(sim.data.1[,c(1,2)])
  vdisp.1 <- getlosvdisp(sim.data.1[,c(1,2,6)])
  num.bins <- dim(vdisp.1)[1]
  if(is.na(vdisp.1[num.bins,1])){
    vdisp.1 <- vdisp.1[-num.bins,]
  }
  #proj.den.2 <- getprojdensity(sim.data.1[,c(1,2)])
  vdisp.2 <- getlosvdisp(sim.data.2[,c(1,2,6)])
  if(is.na(vdisp.2[num.bins,1])){
    vdisp.2 <- vdisp.2[-num.bins,]
  }
  
  chi.total <- sum(vdisp.1,vdisp.2)
  return(abs(suff.stats - chi.total))  
}

###################
#Gamma Galaxy Distance
###################

gamma.dist <- function(theta, suff.stats){
  param.other.1 <- c(2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)
  param.other.2 <- c(1.1, 0.17, 1.5, 0, 8.2, 0.086, 21.0, 3.0)
  
  param.1 <- as.numeric(c(theta[1:2], param.other.1))
  param.2 <- as.numeric(c(theta[3:4], param.other.2))
  
  n <- dim(data)[1]
  sim.data.1 <- do.call(cbind,python.call("sample",param.1,n-1))
  sim.data.2 <- do.call(cbind,python.call("sample",param.2,n-1))
  
  var.1 <- (n-1)/n*var(sim.data.1[,6])
  var.2 <- (n-1)/n*var(sim.data.2[,6])
  
  r.1 <- median(apply(sim.data.1[,1:2],1,magnitude))
  r.2 <- median(apply(sim.data.2[,1:2],1,magnitude))
  
  gamma.sim <- 1+log10(var.2/var.1)/log10(r.2/r.1)
  return(abs(gamma.sim-suff.stats))
}

#######################################
#Vector Distance Functions
#######################################

###################
#Second Moment Galaxy Distance
###################

m1.m2.dist <- function(theta, suff.stats){

  #simulate data
  sim.data <- rnorm(1000,theta[1],theta[2])
  
  #calculate first and second moments
  sim.stats <- c(mean(sim.data),calc.m2(sim.data))
  
  #return difference from summary statistics
  return(matrix(abs(sim.stats-suff.stats),nrow=1))
}


galaxy.m2.a.dist <- function(theta, suff.stats) {

#assume parameters
 param <- as.numeric(c(theta,-5.3,2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5))  
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.d.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0, theta,2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5))
  
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.e.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,theta,0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.Ec.dist <- function(theta, suff.stats) {
  
#assume parameters
  cat("dist 1", "\n")
  cat(theta, "\n")
  param <- as.numeric(c(2.0,-5.3,2.5,theta,1.5, -9.0, 6.9, 0.086, 21.0, 1.5))  
  cat(param, "\n")
#assume 5,000 stars
  n <- 5000
  cat("dist 2", "\n")
  
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sampleR",param,n))[,c(1,2,6)]
  cat("dist 3", "\n")
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  cat("dist 4", "\n")
#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  cat("dist 5", "\n")
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  cat("dist 6", "\n")
  return(epsilon)
}

galaxy.m1m2.Ec.dist <- function(theta, suff.stats) {
  
  #assume parameters
  cat("dist 1", "\n")
  cat(theta, "\n")
  param <- as.numeric(c(2.0,-5.3,2.5,theta,1.5, -9.0, 6.9, 0.086, 21.0, 1.5))  
  cat(param, "\n")
  #assume 5,000 stars
  n <- 5000
  cat("dist 2", "\n")
  
  #call python function to simulate data
  sim.data <- do.call(cbind,python.call("sampleR",param,n))[,c(1,2,6)]
  cat("dist 3", "\n")
  #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  cat("dist 4", "\n")
  #calculate second moments
  moments <- c(mean(sim.data[,1]),calc.m2(sim.data[,2]))
  cat("dist 5", "\n")
  #compute distance
  epsilon <- matrix(abs(moments-suff.stats),nrow=1)
  cat("dist 6", "\n")
  return(epsilon)
}

galaxy.m1m2m2.Ec.dist <- function(theta, suff.stats) {
  
  param <- suff.stats$param
  stars <- suff.stats$stars
  value <- suff.stats$value
  file.name <- suff.stats$file.name
  
  #assume parameters
  cat("dist 1", "\n")
  cat(theta, "\n")
  param <- as.numeric(c(param[1:3],theta,param[5:10]))
  cat(param, "\n")
  
  #set number of stars
  n <- stars
  cat("dist 2", "\n")
  
  #call python function to simulate data
  sim.data <- do.call(cbind,python.call("sampleR",param,n))
  means <- c(as.numeric(apply(sim.data,2,mean)), Sys.time())
  cat(means,"\n",
      sep = "\t", file = file.name, append = TRUE)
  sim.data <- sim.data[,c(1,2,6)]
  cat("dist 3", "\n")
  #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  cat("dist 4", "\n")
  #calculate second moments
  moments <- c(mean(sim.data[,1]),calc.m2(sim.data[,1]),calc.m2(sim.data[,2]))
  cat("dist 5", "\n")
  #compute distance
  epsilon <- matrix(abs(moments-value),nrow=1)
  cat("dist 6", "\n")
  return(epsilon)
}

galaxy.q1.vdisp.dens.Ec.dist <- function(theta, suff.stats) {
  
  param <- suff.stats$param
  stars <- suff.stats$stars
  value <- suff.stats$value
  file.name <- suff.stats$file.name
  
  #assume parameters
  cat("dist 1", "\n")
  cat(theta, "\n")
  param <- as.numeric(c(param[1:3],theta,param[5:10]))
  cat(param, "\n")
  
  #set number of stars
  n <- stars
  cat("dist 2", "\n")
  
  #call python function to simulate data
  sim.data <- do.call(cbind,python.call("sampleR",param,n))
  r <- apply(sim.data[,1:2],1,magnitude)
  min.r <- quantile(r,probs = c(0))
  max.r <- quantile(r,probs = c(.25))
  
  keep <- min.r <= r & r < max.r
  
  dens <- sum(keep) / (max.r^2-min.r^2)
  samp.var <- var(sim.data[keep,6])
  
  #compute distance
  epsilon <- matrix(abs(c(dens,samp.var)-value),nrow=1)
  cat("dist 6", "\n")
  return(epsilon)
}

Ec.full.dist <- function(theta, suff.stats) {
  
  param <- suff.stats$param
  stars <- suff.stats$stars
  value <- suff.stats$value
  file.name <- suff.stats$file.name
  
  #assume parameters
  cat("dist 1", "\n")
  cat(theta, "\n")
  param <- as.numeric(c(param[1:3],theta,param[5:10]))
  cat(param, "\n")
  
  #set number of stars
  n <- stars
  cat("dist 2", "\n")
  
  #call python function to simulate data
  sim.data <- do.call(cbind,python.call("sampleR",param,n))
  r <- apply(sim.data[,1:3],1,magnitude)
  v <- apply(sim.data[,4:6],1,magnitude)
  
  min.r <- quantile(r,probs = c(0))
  max.r <- quantile(r,probs = c(.25))
  
  keep <- min.r <= r & r < max.r
  
  dens <- sum(keep) / (max.r^2-min.r^2)
  samp.var <- var(v[keep])
  
  #compute distance
  epsilon <- matrix(abs(c(dens,samp.var)-value),nrow=1)
  cat("dist 6", "\n")
  return(epsilon)
}

galaxy.m2.rlim.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16,theta,-9.0, 6.9, 0.086, 21.0, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.b.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16, 1.5,theta,6.9, 0.086, 21.0, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.q.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16, 1.5, -9.0, theta,0.086, 21.0, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.Jb.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16, 1.5, -9.0, 6.9, theta, 21.0, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.vmax.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16, 1.5, -9.0, 6.9, 0.086, theta, 1.5))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.rmax.dist <- function(theta, suff.stats) {

#assume parameters
  param <- as.numeric(c(2.0,-5.3,2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, theta))  
 
#assume 2,000 stars
  n <- 2000
 
#call python function to simulate data
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
 #convert to two dimentions
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])

#calculate second moments
  second.moments <- apply(sim.data,2,calc.m2)
  
#compute distance
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.ad.dist <- function(theta, suff.stats) {
  
  param <- as.numeric(c(theta,2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5))  
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.ab.dist <- function(theta, suff.stats) {
  
  param <- as.numeric(c(theta[1],-5.3,2.5, 0.16, 1.5, theta[2], 6.9, 0.086, 21.0, 1.5))  
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

galaxy.m2.db.dist <- function(theta, suff.stats) {
  
  param <- as.numeric(c(2.0,theta[1],2.5, 0.16, 1.5, theta[2], 6.9, 0.086, 21.0, 1.5))  
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",param,n-1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

###################
#Second Moment King Model distance functions
###################

king.sigma.m2.dist <- function(theta, suff.stats) {
  sigma <- as.numeric(theta)
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",sigma,3.0,n+1))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,3.0,20,n+1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

king.ratio.m2.dist <- function(theta, suff.stats) {
  ratio <- as.numeric(theta)
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",5.0,ratio,n+1))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",5.0,ratio,20,n+1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}

king.all.m2.dist <- function(theta, suff.stats) {
  
  sigma <- as.numeric(theta[1])
  ratio <- as.numeric(theta[2])
  
  n <- 2000
  
  sim.data <- do.call(cbind,python.call("sample",sigma,ratio,n+1))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,ratio,20,n+1))[,c(1,2,6)]
  
  sim.data <- cbind(apply(sim.data[,1:2],1,magnitude),sim.data[,3])
  second.moments <- apply(sim.data,2,calc.m2)
  
  epsilon <- matrix(abs(second.moments-suff.stats),nrow=1)
  
  return(epsilon)
}


king.sigma.maxr.dist <- function(theta, suff.stats) {
  sigma <- as.numeric(theta)
  n <- 5000
  
  sim.data <- do.call(cbind,python.call("sample",sigma,3.0,n))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,3.0,20,n+1))[,c(1,2,6)]
  sim.data <- apply(sim.data[,1:2],1,magnitude)
  
  max.r <- max(sim.data)
  epsilon <- matrix(max.r-suff.stats)
  
  return(epsilon)
}

king.ratio.maxr.dist <- function(theta, suff.stats) {
  cat("king 1", "\n")
  ratio <- as.numeric(theta)
  n <- 5000
  
  cat("king 2", "\n")
  cat("ratio:", ratio, "stars:", n+1, "\n")
  sim.data <- do.call(cbind,python.call("sample",5.0,ratio,n))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",5.0,ratio,20,n+1))[,c(1,2,6)]
  cat("king 3", "\n")
  sim.data <- apply(sim.data[,1:2],1,magnitude)
  
  cat("king 4", "\n")
  max.r <- max(sim.data)
  cat("king 5", "\n")
  epsilon <- abs(matrix(max.r-suff.stats))
  
  return(epsilon)
}

king.all.maxr.dist <- function(theta, suff.stats) {
  
  cat("king 1", "\n")
  sigma <- as.numeric(theta[1])
  ratio <- as.numeric(theta[2])
  
  n <- 5000
  
  cat("king 2", "\n")
  cat("sigma:", sigma, "ratio:", ratio, "stars:", n, "\n")
  sim.data <- do.call(cbind,python.call("sample",sigma,ratio,n+1))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,ratio,20,n+1))[,c(1,2,6)]
  cat("king 3", "\n")
  sim.data <- apply(sim.data[,1:2],1,magnitude)
  cat("king 4", "\n")
  max.r <- max(sim.data)
  cat("king 5", "\n")
  epsilon <- abs(matrix(max.r-suff.stats))
  cat("king 6", "\n")
  return(epsilon)
}

king.all.maxrm2v.dist <- function(theta, suff.stats) {
  
  time.1 <- Sys.time()
  cat("king 1", time.1, "\n")
  sigma <- as.numeric(theta[1])
  ratio <- as.numeric(theta[2])
  
  n <- 5000
  
  time.2 <- Sys.time()
  cat("king 2", time.2, time.2-time.1, "\n")
  
  cat("sigma:", sigma, "ratio:", ratio, "stars:", n+1, "\n")
  sim.data <- do.call(cbind,python.call("sample",sigma,ratio,n))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,ratio,20,n+1))[,c(1,2,6)]
  time.3 <- Sys.time()
  cat("king 3", time.3, time.3-time.2, "\n")
  
  radius <- apply(sim.data[,1:2],1,magnitude)
  v.disp <- calc.m2(sim.data[,3])
  time.4 <- Sys.time()
  cat("king 4", time.4, time.4-time.3, "\n")
  
  max.r <- max(radius)
  time.5 <- Sys.time()
  cat("king 5", time.5, time.5-time.4, "\n")
  
  epsilon <- abs(c(max.r-suff.stats[1], v.disp-suff.stats[2]))
  time.6 <- Sys.time()
  cat("king 6", time.6, time.6-time.5, "\n")
  
  return(epsilon)
}

king.all.maxrm2v.dist.times <- function(theta, suff.stats) {
  
  time.1 <- Sys.time()
  cat("king 1", time.1, "\n")
  sigma <- as.numeric(theta[1])
  ratio <- as.numeric(theta[2])
  
  n <- 5000
  
  time.2 <- Sys.time()
  cat("king 2", time.2, time.2-time.1, "\n")
  
  cat("sigma:", sigma, "ratio:", ratio, "stars:", n+1, "\n")
  sim.data <- do.call(cbind,python.call("sample",sigma,ratio,n))[,c(1,2,6)]
  #sim.data <- do.call(cbind,python.call("sample2",sigma,ratio,20,n+1))[,c(1,2,6)]
  time.3 <- Sys.time()
  cat("king 3", time.3, time.3-time.2, "\n")
  
  radius <- apply(sim.data[,1:2],1,magnitude)
  v.disp <- calc.m2(sim.data[,3])
  time.4 <- Sys.time()
  cat("king 4", time.4, time.4-time.3, "\n")
  
  max.r <- max(radius)
  time.5 <- Sys.time()
  cat("king 5", time.5, time.5-time.4, "\n")
  
  epsilon <- abs(c(max.r-suff.stats[1], v.disp-suff.stats[2]))
  time.6 <- Sys.time()
  cat("king 6", time.6, time.6-time.5, "\n")
  
  return(list(epsilon=epsilon,times=c(time.1,time.2,time.3,time.4,time.5,time.6)))
}


ks.distance <-  function(theta, suff.stats){
  python.fun <- suff.stats[["python.fun"]]
  param <- suff.stats[["param"]]
  steps <- suff.stats[["steps"]]
  samplesize <- suff.stats[["samplesize"]]
  resamplefactor <- suff.stats[["resamplefactor"]]
  
  data.true <- suff.stats[["data.true"]]
  indicies <- suff.stats[["indicies"]]
  param[indicies] <- theta
  param <- as.numeric(param)
  
  data.sim <- do.call(cbind,python.call(python.fun,param,steps,samplesize,resamplefactor))
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) ks.test(data1[,n],data2[,n])$statistic,data.true,data.sim))
  ks.r <- sqrt(sum(ks.vec[1:3]^2))
  ks.v <- sqrt(sum(ks.vec[4:6]^2))
  return(c(ks.r,ks.v))
}

ks.distance.smooth <-  function(theta, suff.stats){
  ks.smooth <- function(x,y){
    x.dens <- density(x)
    n <- length(x.dens$x)
    y.dens <- density(y, from = x.dens$x[1], to = x.dens$x[n], n=n)
    x.cum <- cumsum(x.dens$y)/sum(x.dens$y)
    y.cum <- cumsum(y.dens$y)/sum(y.dens$y)
    return(max(abs(x.cum-y.cum)))
  }
  python.fun <- suff.stats[["python.fun"]]
  param <- suff.stats[["param"]]
  steps <- suff.stats[["steps"]]
  samplesize <- suff.stats[["samplesize"]]
  resamplefactor <- suff.stats[["resamplefactor"]]
  
  data.true <- suff.stats[["data.true"]]
  indicies <- suff.stats[["indicies"]]
  param[indicies] <- theta
  param <- as.numeric(param)
  
  data.sim <- do.call(cbind,python.call(python.fun,param,steps,samplesize,resamplefactor))
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) ks.smooth(data1[,n],data2[,n]),data.true,data.sim))
  ks.r <- sqrt(sum(ks.vec[1:3]^2))
  ks.v <- sqrt(sum(ks.vec[4:6]^2))
  return(c(ks.r,ks.v))
}

ls.distance <-  function(theta, suff.stats){
  l2.smooth <- function(x1,x2){
    n <- length(x1)
    dens1 <- density(x1)
    dens2 <- density(x2, from = min(dens1$x), to = max(dens1$x), n = length(dens1$x))
    y1 <- cumsum(dens1$y)/sum(dens1$y)
    y2 <- cumsum(dens2$y)/sum(dens2$y)
    return(sqrt(sum((y1-y2)^2)))
  }
  python.fun <- suff.stats[["python.fun"]]
  param <- suff.stats[["param"]]
  steps <- suff.stats[["steps"]]
  samplesize <- suff.stats[["samplesize"]]
  resamplefactor <- suff.stats[["resamplefactor"]]
  
  data.true <- suff.stats[["data.true"]]
  indicies <- suff.stats[["indicies"]]
  param[indicies] <- theta
  param <- as.numeric(param)
  
  data.sim <- do.call(cbind,python.call(python.fun,param,steps,samplesize,resamplefactor))
  ks.vec <- unlist(lapply(1:6,function(n,data1,data2) l2.smooth(data1[,n],data2[,n]),data.true,data.sim))
  ks.r <- sqrt(sum(ks.vec[1:3]^2))
  ks.v <- sqrt(sum(ks.vec[4:6]^2))
  return(c(ks.r,ks.v))
}

rv2.distance <-  function(theta, suff.stats){
  l2.smooth <- function(x1,x2){
    n <- length(x1)
    dens1 <- density(x1)
    dens2 <- density(x2, from = min(dens1$x), to = max(dens1$x), n = length(dens1$x))
    y1 <- cumsum(dens1$y)/sum(dens1$y)
    y2 <- cumsum(dens2$y)/sum(dens2$y)
    return(sqrt(sum((y1-y2)^2)))
  }
  rv2.dist <- function(data1,data2){
    r1 <- apply(data1[,1:3],1,function(y) sqrt(sum(y^2)))
    v1 <- apply(data1[,4:6],1,function(y) sqrt(sum(y^2)))
    
    r2 <- apply(data2[,1:3],1,function(y) sqrt(sum(y^2)))
    v2 <- apply(data2[,4:6],1,function(y) sqrt(sum(y^2)))
    
    x1 <- r1*(v1^2)
    x2 <- r2*(v2^2)
    
    return(l2.smooth(x1,x2))
  }
  python.fun <- suff.stats[["python.fun"]]
  param <- suff.stats[["param"]]
  steps <- suff.stats[["steps"]]
  samplesize <- suff.stats[["samplesize"]]
  resamplefactor <- suff.stats[["resamplefactor"]]
  
  data.true <- suff.stats[["data.true"]]
  indicies <- suff.stats[["indicies"]]
  param[indicies] <- theta
  param <- as.numeric(param)
  
  data.sim <- do.call(cbind,python.call(python.fun,param,steps,samplesize,resamplefactor))
  return(rv2.dist(data.true,data.sim))
}

dens2d.distance <-  function(theta, suff.stats){
  library(MASS)
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
  python.fun <- suff.stats[["python.fun"]]
  param <- suff.stats[["param"]]
  steps <- suff.stats[["steps"]]
  samplesize <- suff.stats[["samplesize"]]
  resamplefactor <- suff.stats[["resamplefactor"]]
  
  data.true <- suff.stats[["data.true"]]
  indicies <- suff.stats[["indicies"]]
  param[indicies] <- theta
  param <- as.numeric(param)
  
  data.sim <- do.call(cbind,python.call(python.fun,param,steps,samplesize,resamplefactor))
  return(dens2d.dist(data.true,data.sim))
}