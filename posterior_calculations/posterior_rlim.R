###################
#Advanced Data Analysis
#Brendan McVeigh
#July 28, 2015
###################

setwd("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/posterior_calculations")

#########################################################
#Define Functions
#########################################################

vesc <- function(x, r.lim = 1.5, vmax = 21, rmax = 1.5){
  rs = rmax/2.16
  Ps = (vmax/0.465)^2
  
  xlim = r.lim/rs  	 #turn rlim in unit of rs
  Pr   = Ps * ( 1 - (log(1+x))/x )
  Plim = Ps * ( 1 - (log(1+xlim))/xlim ) #lms 0.45*Ps
  
  vesc = (2 * (Plim - Pr))^0.5
  return(vesc)
}

hprob <- function(v, x, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                  b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  
  phi.s <- (vmax/.465)^2
  r.s <- rmax/2.16  
  phi.r <- phi.s * (1-log(1+x)/x)
  
  E <- v^2/2 + phi.r
  Ec <- Ec * phi.s
  xlim <- r.lim/r.s
  phi.lim <- phi.s * (1-log(1+xlim)/xlim)
  
  if(E<phi.lim){
    hE <- E^a*(E^q+Ec^q)^(d/q)*(phi.lim-E)^e
  } else hE = 0
  return(hE*v^2*x^2)
}

hprob.vec.v <- Vectorize(hprob, vectorize.args = "v")

intv.dens <- function(x, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                      b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  out <- integrate(hprob.vec.v, lower = 0, upper = vesc(x), x, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                   q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  return(out$value)
}

intv.dens.vec.x <- Vectorize(intv.dens, vectorize.args = "x")

norm <- function(a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                 b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  N <- integrate(intv.dens.vec.x, lower = 0, upper = r.lim/r.s, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
          q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)$value
  return(N)
}

hprob.norm <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                       b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  x <- xv[1]
  v <- xv[2]
  N <- norm(a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)
  prob <- hprob(v, x, a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)/N
  return(prob)
}

hprob.xv <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                  b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  x <- xv[1]
  v <- xv[2]
  phi.s <- (vmax/.465)^2
  r.s <- rmax/2.16  
  phi.r <- phi.s * (1-log(1+x)/x)
  
  E <- v^2/2 + phi.r
  Ec <- Ec * phi.s
  xlim <- r.lim/r.s
  phi.lim <- phi.s * (1-log(1+xlim)/xlim)
  
  if(E<phi.lim){
    hE <- E^a*(E^q+Ec^q)^(d/q)*(phi.lim-E)^e
  } else hE = 0
  return(hE*v^2*x^2)
}

nll.norm <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                     b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  
  liklihood <- apply(xv,1,hprob.xv, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                     q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax) / norm(a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)
  ll <- sum(log(liklihood))
  return(-ll)
}

nll.norm.rlim <- Vectorize(function(r.lim,xv){nll.norm(xv,r.lim=r.lim)}, vectorize.args = "r.lim")

density.adj <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                     b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5, adj = 0){
  nll <- nll.norm(xv=xv, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                 q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  
  return(exp((adj-nll)))
}

density.adj.rlim <- Vectorize(function(r.lim,xv,adj){density.adj(xv,r.lim=r.lim,adj=adj)}, vectorize.args = "r.lim")

######################################
#Set parameters
######################################

a <- 2.0
d <- -5.3
e <- 2.5
Ec <- .16
r.lim <- 1.5
b <- -9.0
q <- 6.9
Jb <- 8.6*10^(-2)
vmax <- 21
rmax <- 1.5

phi.s <- (vmax/.465)^2
r.s <- rmax/2.16
xlim <- r.lim / r.s

######################################
#Norm Calculation
######################################

#Ec.grid <- seq(.01,1,.001)
rlim.grid <- unique(c(seq(1,5,.1),seq(1,2,.05)))
rlim.grid <- sort(rlim.grid)


#########################################################
#Calculate Posterior
#########################################################

data.rlim <- read.delim("../data/simulated_galaxy_reduced_model.txt")

r <- sqrt(data.rlim[,1]^2+data.rlim[,2]^2+data.rlim[,3]^2)
x <- r / r.s
v <- sqrt(data.rlim[,4]^2+data.rlim[,5]^2+data.rlim[,6]^2)
xv <- cbind(x,v)

#MLE Adjustment
mle.full <- optim(1.5,nll.norm.rlim,xv=xv, method = 'L-BFGS-B', lower = 1, upper = 12)
denom.int <- integrate(density.adj.rlim,lower=.1, upper=.3,xv=xv, adj = mle.full$value)$value

posterior <- density.adj.rlim(r.lim=rlim.grid, xv=xv, adj = mle.full$value)/denom.int

# #Plot posterior density
# pdf(file="Analytic_Posterior.pdf")
# plot(Ec.grid,posterior,"l", ylab = "Posterior Density", xlab = "Ec", main="Posterior Density")
# abline(v=.16, lty=2,lwd=2)
# dev.off()
# 
# pdf(file="Analytic_Posterior_Zoom.pdf")
# plot(Ec.grid,posterior,"l", xlim = c(.14,.18), ylab = "Posterior Density", xlab = "Ec", main="Posterior Density")
# abline(v=.16, lty=2,lwd=2)
# dev.off()

save(data.E,mle.full,denom.int,Ec.grid,posterior,file="Posterior_Results.R")

###################
#Compare with run #1
###################
load("Posterior_Results.R")
load("Ec_all_15_add_results.R")
n <- length(results)
parameters <- lapply(results,function(x) x[["parameters"]])
# weights <-do.call(cbind,lapply(results,function(x) x[["weights"]]))
# distances <- do.call(cbind,lapply(results,function(x) x[["distances"]]))
# samples <- do.call(cbind,lapply(results,function(x) x[["total.samples"]]))
# cutoff <- do.call(rbind,lapply(results,function(x) x[["epsilon"]]))

pdf(file="Posterior_Comparison.pdf")
plot(density(results[[16]][["parameters"]]),lwd = 2, ylim=c(0,250),
     main="Density of Parameter \"Ec\"", xlim=c(.13,.18),lty=2)
lines(Ec.grid,posterior,lwd = 2)
#abline(v=.16, lty=3,lwd=2)
abline(v=mle.full$par, lty=3,lwd=2)
legend(x=.13,y=230,
       legend=c("Analytic Posterior", "ABC Posterior", "MLE"),
       lty=1:3, lwd = rep(2,3))
dev.off()

pdf(file="ABC_Densities.pdf")
plot(density(parameters[[1]]),lwd = 2, ylim=c(0,150),
     xlim = c(.1,.22),
     main=expression(paste("Density of Ec")),
     xlab = expression(paste("Ec")))
for(i in 2:5){
  lines(density(parameters[[seq(0,16,4)[i]]]),lwd = 2,col=i)
}
abline(v=.16,lwd=2,col=1,lty=2)
abline(v=mle.full$par,lwd=2,col=1,lty=3)
legend(x=.18,y=150,
       legend=c(paste("Time Step",c(1,seq(4,16,4))),"True Value", "MLE"),
       col = c(1:5,1,1), lty=c(rep(1,5),2,3), lwd = rep(2,5+2))
dev.off()

pdf(file="ABC_Densities_tail.pdf")
start <- 12
series <- n - start + 1
plot(density(parameters[[start]]),lwd = 2, ylim=c(0,150),
     main=expression(paste("Density of Ec")),
     xlab = expression(paste("Ec")), col = 1)
for(i in (start+1):n){
  lines(density(parameters[[i]]),lwd = 2,col=i-start+1)
}
abline(v=.16,lwd=2,col=1,lty=2)
abline(v=mle.full$par,lwd=2,col=1,lty=3)
legend(x=.162,y=150,
       legend=c(paste("Time Step",as.character(start:n)),"True Value", "MLE"),
       col = c(1:series,1,1), lty=c(rep(1,series),2,3), lwd = rep(2,series+2))
dev.off()

plot(2:16,apply(samples,2,sum)[-1])